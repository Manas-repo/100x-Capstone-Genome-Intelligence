"""Model step (T11) and verdict logic, plus the end-to-end pipeline.

Model step: for each lab's free-text reasoning, ask whether that reasoning carries over to THIS
person (an apparently healthy, unrelated South Asian adult with no phenotype or family data).
A rule cannot do this: the reasoning is unstructured text, and whether "seen in an affected
family" or "found in trans with a pathogenic variant" applies depends on reading it.

Verdict: plain code over the rule flags and the model's per-lab answers, so it can be explained.
"""
import json
import os
from pathlib import Path

from dotenv import load_dotenv

from evidence import assemble, flags
from vcf_clinvar import match_clinvar, read_vcf

ROOT = Path(__file__).resolve().parent.parent
load_dotenv(ROOT / ".env")
CACHE = ROOT / "data" / "llm_cache.json"
MAX_SUPPORT, MAX_DISPUTE, MAX_CHARS = 8, 4, 1200
SHOWN = 5  # findings shown as cards; the rest are listed as declined or held back

PERSON = (
    "An apparently healthy adult of South Asian (Indian) ancestry, sequenced as part of a population "
    "study. No symptoms, no clinical phenotype and no family members were tested. Genotype at this "
    "variant: {zygosity}."
)

SYSTEM = """You help a clinical variant curator check ClinVar evidence. Labs classified a variant.
For each lab's written reasoning, decide whether that reasoning carries over to the person described.
Evidence tied to the lab's own patient (their phenotype, a de novo finding in their trio, segregation
in their family, a second variant in trans in their patient) does not carry over to an unrelated
healthy person. Evidence about the variant itself (functional studies, population rarity, many
independent affected carriers, an expert panel) can carry over. Be strict and brief. Quote the
lab's words in reasons. Never invent evidence the text does not state. Reply in JSON only."""

USER = """Variant: {name} in {gene} ({consequence}). Gene disease mechanism: {mechanism}.
Condition(s): {conditions}
Person: {person}

Lab submissions (id, call, year, method, reasoning):
{subs}

Return JSON:
{{"labs": [{{"id": <id>, "evidence_types": [<from: population_rarity, case_observations, segregation,
 de_novo, in_trans, functional_study, computational, null_variant, same_codon, founder_effect,
 expert_review, none_stated>], "carries_over": "yes" | "partly" | "no" | "cant_tell",
 "reason": "<one sentence quoting the lab>"}}],
 "strongest_transferable_evidence": "<one sentence, or 'none'>",
 "cheapest_check": "<the single quickest thing the curator should check to confirm or reject this, one sentence>"}}"""


def _cache() -> dict:
    try:
        return json.loads(CACHE.read_text())
    except (FileNotFoundError, ValueError):
        return {}


def model_review(ev: dict) -> dict | None:
    """Per-lab carry-over judgments from the Groq model. None if no key is set."""
    key = os.getenv("GROQ_API_KEY")
    if not key:
        return None
    cache = _cache()
    ck = str(ev["variation_id"])
    if ck in cache:
        return cache[ck]

    from groq import Groq

    subs = ev["support"][:MAX_SUPPORT] + ev["dispute"][:MAX_DISPUTE]
    lines = []
    for i, s in enumerate(subs):
        text = (s["comment"] or "(no reasoning given)")[:MAX_CHARS].replace("\n", " ")
        lines.append(f"[{i}] {s['call']} | {s['year']} | {s['method']} | {s['submitter']}: {text}")
    prompt = USER.format(
        name=ev["name"], gene=ev["gene"], consequence=ev["consequence"], mechanism=ev["mechanism_note"],
        conditions=ev["conditions"], person=PERSON.format(zygosity=ev["zygosity"]), subs="\n".join(lines),
    )
    resp = Groq(api_key=key).chat.completions.create(
        model=os.getenv("GROQ_MODEL", "openai/gpt-oss-120b"),
        response_format={"type": "json_object"},
        temperature=0,
        messages=[{"role": "system", "content": SYSTEM}, {"role": "user", "content": prompt}],
    )
    out = json.loads(resp.choices[0].message.content)
    for lab in out.get("labs", []):  # attach submitter names back by index
        i = lab.get("id")
        if isinstance(i, int) and 0 <= i < len(subs):
            lab["submitter"], lab["call"] = subs[i]["submitter"], subs[i]["call"]
    cache[ck] = out
    tmp = CACHE.with_suffix(".tmp")
    tmp.write_text(json.dumps(cache, indent=1))
    tmp.replace(CACHE)
    return out


def verdict(ev: dict, fl: list[dict], review: dict | None) -> tuple[str, list[str]]:
    """Defensible / Contested / Declined, with the reasons that decided it."""
    declines = [f["reason"] for f in fl if f["effect"] == "decline"]
    if declines:
        return "Declined", declines
    contests = [f["reason"] for f in fl if f["effect"] == "contest"]
    expert = any(f["code"] == "expert_panel" for f in fl)

    carried = None
    if review:
        sup_labs = [l for l in review.get("labs", []) if l.get("call") in ("P", "LP")]
        carried = sum(l.get("carries_over") in ("yes", "partly") for l in sup_labs)
        if sup_labs and carried == 0:
            contests.append("None of the pathogenic labs' stated reasoning carries over to this person.")

    if expert and not any(f["code"] in ("conflict", "null_not_lof", "sas_enriched") for f in fl):
        return "Defensible", ["Expert panel review with no conflicting labs."] + contests
    if not contests:
        return "Defensible", ["Several labs agree, the evidence is recent, and it carries over to this person."]
    return "Contested", contests


def analyse(vcf_path) -> dict:
    """End to end: VCF -> candidates -> evidence -> flags -> model -> verdicts."""
    hits = match_clinvar(read_vcf(vcf_path))
    # Candidates: anything at least one lab calls P/LP (headline label alone would hide conflicts).
    cand = hits[hits["ClinicalSignificance"].str.contains("athogenic", na=False)
                & ~hits["ClinicalSignificance"].str.contains("^Benign|^Likely benign", na=False)]
    results = []
    for _, row in cand.iterrows():
        ev = assemble(row)
        if not ev["support"]:
            continue
        fl = flags(ev)
        review = None
        if not any(f["effect"] == "decline" for f in fl):
            review = model_review(ev)  # only spend model calls on variants still in play
        v, why = verdict(ev, fl, review)
        results.append({"evidence": ev, "flags": fl, "review": review, "verdict": v, "why": why})

    order = {"Defensible": 0, "Contested": 1, "Declined": 2}
    results.sort(key=lambda r: (order[r["verdict"]], -len(r["evidence"]["support"])))
    live = [r for r in results if r["verdict"] != "Declined"]
    return {
        "findings": live[:SHOWN],
        "held_back": live[SHOWN:],
        "declined": [r for r in results if r["verdict"] == "Declined"],
        "n_variants": len(hits),
        "model_used": any(r["review"] for r in results),
    }
