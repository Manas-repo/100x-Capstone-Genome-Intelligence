"""Evidence assembly (T9) and rule checks (T10) for one candidate variant.

The rule checks encode how a practising curator (interview 2026-10-03) says he
weighs a ClinVar record: count of submitters, whether each submitter's criteria apply, evidence
age, the null-variant-in-a-non-LoF-gene trap, plus inheritance and population fit for this person.
No model is used here; every flag is a plain rule with a plain-language reason.
"""
import re
from datetime import date

import pandas as pd

from gnomad import lookup as gnomad_lookup
from vcf_clinvar import panel, submissions_for

TODAY = date(2026, 10, 5)
STALE_YEARS = 5

CALL = {
    "pathogenic": "P", "likely pathogenic": "LP", "pathogenic/likely pathogenic": "P",
    "uncertain significance": "VUS", "likely benign": "LB", "benign": "B", "benign/likely benign": "B",
}
SUPPORT, DISPUTE = {"P", "LP"}, {"VUS", "LB", "B"}


def _call(label: str) -> str | None:
    return CALL.get(str(label).strip().lower())


def _year(d: str) -> int | None:
    m = re.search(r"(\d{4})", str(d))
    return int(m.group(1)) if m else None


def consequence(name: str) -> str:
    """Rough variant type from the ClinVar HGVS name."""
    n = str(name)
    if re.search(r"p\.\([^)]*(Ter|\*)\)|p\.[A-Za-z]{3}\d+(Ter|\*)", n) or "fs" in n:
        return "null (nonsense/frameshift)"
    if re.search(r"c\.[\d*-]+[+-][12][ACGTdi]", n):
        return "null (canonical splice)"
    if re.search(r"p\.\(?[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)?", n):
        return "missense"
    if "=" in n:
        return "synonymous"
    return "other"


def assemble(row: pd.Series, tag: str = "", with_population: bool = True) -> dict:
    """row: one sample variant matched to ClinVar (from vcf_clinvar.match_clinvar).
    tag: ClinVar release to read submissions from ('' = current, '2020-01' for the backtest)."""
    vid = int(row["VariationID"])
    subs = submissions_for(vid, tag)
    subs["call"] = subs["ClinicalSignificance"].map(_call)
    subs["year"] = subs["DateLastEvaluated"].map(_year)
    subs = subs.sort_values("year", ascending=False, na_position="last")

    gene = str(row["GeneSymbol"]).split(";")[0]
    pinfo = panel().set_index("gene")
    g = pinfo.loc[gene].to_dict() if gene in pinfo.index else {}

    pop = gnomad_lookup(str(row["chrom"]), int(row["pos"]), row["ref"], row["alt"]) if with_population else {}

    def _sub(s):
        return {
            "submitter": s["Submitter"], "call": s["call"] or s["ClinicalSignificance"],
            "year": s["year"], "method": s["CollectionMethod"], "review": s["ReviewStatus"],
            "comment": None if s["Description"] in (None, "-", "") else s["Description"],
            "condition": s["ReportedPhenotypeInfo"], "scv": s["SCV"],
        }

    return {
        "variation_id": vid,
        "gene": gene,
        "name": row["Name"],
        "consequence": consequence(row["Name"]),
        "zygosity": row.get("zygosity"),
        "call_filter": row.get("filter"),
        "clinvar_label": row["ClinicalSignificance"],
        "review_status": row["ReviewStatus"],
        "conditions": g.get("conditions", row.get("PhenotypeList")),
        "inheritance": g.get("inheritance"),
        "panel": g.get("panel", "acmg_sf"),
        "finding_type": (
            "carrier (one copy of a recessive variant: matters for children, not for this person's health)"
            if g.get("inheritance") == "AR" and row.get("zygosity") == "heterozygous"
            else "health risk for this person"
        ),
        "lof_mechanism": g.get("lof_mechanism"),
        "mechanism_note": g.get("mechanism_note"),
        "acmg_restriction": g.get("acmg_restriction") or "",
        "support": [_sub(s) for _, s in subs[subs["call"].isin(SUPPORT)].iterrows()],
        "dispute": [_sub(s) for _, s in subs[subs["call"].isin(DISPUTE)].iterrows()],
        "population": {**pop, "af_1kg_sas": row.get("af_1kg_sas"), "af_1kg_all": row.get("af_1kg_all")},
        "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/",
    }


def flags(ev: dict, as_of_year: int = TODAY.year, person: bool = True) -> list[dict]:
    """Rule checks. Each flag: code, effect (decline / contest / support), plain reason."""
    out = []

    def add(code, effect, reason):
        out.append({"code": code, "effect": effect, "reason": reason})

    sup, dis = ev["support"], ev["dispute"]
    labs = {s["submitter"] for s in sup}

    # What kind of finding is this for this person? (skipped in the variant-only backtest)
    # One copy in a recessive gene = carrier finding (reproductive relevance), not a health risk.
    # Two copies in a healthy adult = a claim the person should be affected; if they are not,
    # that itself argues against the classification (ACMG BS2), a low-penetrance condition, or a bad call.
    if person and ev["inheritance"] == "AR" and ev["zygosity"] == "homozygous":
        add("healthy_homozygote", "contest",
            f"This apparently healthy adult has two copies. If {ev['gene']} disease were fully penetrant "
            "and early-onset, they would be expected to be affected. Check the condition's age of onset "
            "and penetrance before trusting the 'pathogenic' label.")
    if ev["acmg_restriction"] == "truncating variants only" and not ev["consequence"].startswith("null"):
        add("acmg_restriction", "decline", f"ACMG only reports truncating {ev['gene']} variants; this one is not truncating.")
    if ev["gene"] == "HFE" and not ("Cys282Tyr" in ev["name"] and ev["zygosity"] == "homozygous"):
        add("acmg_restriction", "decline", "ACMG only reports HFE p.C282Y when homozygous.")
    if ev["call_filter"] not in (None, "PASS", "."):
        add("low_quality_call", "decline", f"The sequencing call did not pass quality filters ({ev['call_filter']}).")

    # Population fit.
    p = ev["population"]
    if not person:
        pass
    elif p.get("af_sas", 0) >= 0.01 and ev["inheritance"] != "AR":
        add("too_common_sas", "decline",
            f"Found in {p['af_sas']:.1%} of South Asian alleles in gnomAD, too common to cause a rare dominant disease.")
    elif p.get("af_sas", 0) >= 0.05:
        add("too_common_sas", "decline",
            f"Found in {p['af_sas']:.1%} of South Asian alleles in gnomAD. Above 5% is treated as benign "
            "even for recessive disease (ACMG BA1).")
    elif p.get("af_sas", 0) >= 0.001 and p.get("af_sas", 0) > 5 * max(p.get("af_all", 0), 1e-6):
        add("sas_enriched", "contest",
            f"{p['af_sas'] / max(p['af_all'], 1e-6):.0f}x more common in South Asians than overall "
            f"({p['af_sas']:.2%} vs {p['af_all']:.3%}). Evidence measured elsewhere may not fit.")
    if person and p.get("found") is None:
        add("no_sas_data", "contest", "gnomAD could not be reached, so South Asian frequency is unknown.")
    elif person and (not p.get("found") or p.get("an_sas", 0) < 1000):
        add("no_sas_data", "contest", "No usable South Asian frequency data in gnomAD for this variant.")

    # The curator's pattern: null variant called LP in a gene where disease is not loss of function.
    if ev["consequence"].startswith("null") and ev["lof_mechanism"] in ("no", "mixed"):
        add("null_not_lof", "contest",
            f"Null variant, but {ev['gene']} disease works by {ev['mechanism_note']}. "
            "'Null so pathogenic' (PVS1) may not apply.")

    # Strength and provenance of the ClinVar evidence itself.
    if "expert panel" in str(ev["review_status"]) or "practice guideline" in str(ev["review_status"]):
        add("expert_panel", "support", "Reviewed by a ClinGen expert panel.")
    if len(labs) == 1:
        add("single_lab", "contest", f"Only one lab ({next(iter(labs))}) calls it pathogenic.")
    if len(labs) == 0:
        add("no_pathogenic_submitter", "decline", "No submitter currently calls it pathogenic.")
    if dis:
        n_vus = sum(d["call"] == "VUS" for d in dis)
        n_ben = sum(d["call"] in ("LB", "B") for d in dis)
        add("conflict", "contest", f"{len(labs)} lab(s) say pathogenic; {n_vus} say uncertain, {n_ben} say benign.")
    years = [s["year"] for s in sup if s["year"]]
    if years and as_of_year - max(years) > STALE_YEARS:
        add("stale", "contest", f"Newest pathogenic assessment is from {max(years)}, over {STALE_YEARS} years old.")
    if sup and all(s["method"] in ("literature only", "not provided") for s in sup):
        add("literature_only", "contest", "Pathogenic calls come only from literature, not clinical testing.")
    if sup and all(s["comment"] is None for s in sup):
        add("no_reasoning", "contest", "No pathogenic submitter explains its reasoning.")
    return out
