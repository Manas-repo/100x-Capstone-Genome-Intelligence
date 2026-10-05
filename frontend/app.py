"""Gradio frontend (T12c). Talks to the FastAPI backend at BACKEND_URL.

Run locally (backend running on :8000):  python frontend/app.py
"""
import html
import os
import uuid

import gradio as gr
import requests

BACKEND = os.getenv("BACKEND_URL", "http://localhost:8000").rstrip("/")
SLOTS = 5  # matches SHOWN in core/adjudicate.py
TIMEOUT = 300  # first call can wait on a sleeping Render instance plus model calls

NOTE = (
    "**Not a diagnosis.** This tool helps a variant curator decide which ClinVar 'pathogenic' calls "
    "deserve trust for one person. It is a research prototype built on public 1000 Genomes data. "
    "Nothing here should be used for medical decisions."
)
COLORS = {"Defensible": "#1f7a3a", "Contested": "#a15c00", "Declined": "#666"}


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def _lab_line(s: dict, model: dict | None) -> str:
    comment = s.get("comment") or "no reasoning given"
    if len(comment) > 260:
        comment = comment[:260] + "..."
    carry = ""
    if model:
        tag = {"yes": "carries over", "partly": "partly carries over", "no": "does not carry over",
               "cant_tell": "can't tell"}.get(model.get("carries_over"), "")
        if tag:
            carry = f" <b>[{esc(tag)}: {esc(model.get('reason', ''))}]</b>"
    return (f"<li><b>{esc(s['call'])}</b>, {esc(s['year'] or 'no date')}, {esc(s['submitter'])} "
            f"({esc(s['method'])}): <i>{esc(comment)}</i>{carry}</li>")


def card_html(r: dict) -> str:
    e, v = r["evidence"], r["verdict"]
    review = r.get("review") or {}
    by_sub = {l.get("submitter"): l for l in review.get("labs", []) if l.get("submitter")}
    pop = e["population"]
    years = [s["year"] for s in e["support"] + e["dispute"] if s["year"]]
    age = f"{min(years)} to {max(years)}" if years else "no dates given"
    if pop.get("found"):
        popline = (f"gnomAD South Asian {pop['af_sas']:.4%} (from {pop['an_sas']:,} alleles) "
                   f"vs all populations {pop['af_all']:.4%}.")
    else:
        popline = "No South Asian frequency data in gnomAD."
    sup = "".join(_lab_line(s, by_sub.get(s["submitter"])) for s in e["support"][:6]) or "<li>none</li>"
    dis = "".join(_lab_line(s, by_sub.get(s["submitter"])) for s in e["dispute"][:4]) or "<li>none</li>"
    why = "".join(f"<li>{esc(w)}</li>" for w in r["why"])
    check = review.get("cheapest_check")
    return f"""
<div style="border:1px solid #ccc;border-left:6px solid {COLORS[v]};border-radius:8px;padding:12px 16px">
  <div style="font-size:1.1em"><span style="background:{COLORS[v]};color:#fff;padding:2px 8px;border-radius:4px">{v}</span>
    &nbsp;<b>{esc(e['gene'])}</b> {esc(e['name'])}</div>
  <div style="margin:6px 0">{esc(e['consequence'])}, {esc(e['zygosity'])}. {esc(e['finding_type'])}.<br>
    Condition: {esc(e['conditions'])}. ClinVar headline: <i>{esc(e['clinvar_label'])}</i> ({esc(e['review_status'])}).</div>
  <b>Why this verdict</b><ul>{why}</ul>
  <b>Supporting (pathogenic) calls</b><ul>{sup}</ul>
  <b>Disputing calls</b><ul>{dis}</ul>
  <b>Evidence age:</b> {age}. &nbsp; <b>Population match:</b> {popline}<br>
  {f'<b>Quickest check:</b> {esc(check)}<br>' if check else ''}
  <a href="{e['clinvar_url']}" target="_blank">Open the ClinVar record</a>
</div>"""


def short(r: dict) -> str:
    e = r["evidence"]
    return f"{e['gene']} {e['name'][:70]} | {e['variation_id']}"


def load(choice, upload):
    try:
        if upload is not None:
            with open(upload, "rb") as f:
                resp = requests.post(f"{BACKEND}/analyse", files={"file": (os.path.basename(upload), f)},
                                     timeout=TIMEOUT)
        else:
            sample = "heldout" if choice.startswith("Held") else "dev"
            resp = requests.post(f"{BACKEND}/analyse/{sample}", timeout=TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as ex:
        msg = getattr(getattr(ex, "response", None), "text", "") or str(ex)
        return [None, f"**Error:** {msg}"] + [gr.update(visible=False)] * SLOTS * 2 + ["", gr.update(choices=[]), [],
                                                                                      gr.update(choices=[])]

    sysr = data["system"]
    cards = sysr["findings"]
    summary = (f"**{data['sample']}**: {sysr['n_variants']:,} variants matched ClinVar in the panel genes. "
               f"The naive list reports **{len(data['baseline'])}** findings. This system shows "
               f"**{len(cards)}** findings, holds back {len(sysr['held_back'])} more, and declines "
               f"**{len(sysr['declined'])}**. Model step: {'on' if sysr['model_used'] else 'off (rules only)'}.")
    if not cards:
        summary += "\n\nNo finding survived. For a healthy person that is the most common correct answer."

    htmls, groups = [], []
    for i in range(SLOTS):
        if i < len(cards):
            htmls.append(gr.update(value=card_html(cards[i]), visible=True))
            groups.append(gr.update(visible=True))
        else:
            htmls.append(gr.update(value="", visible=False))
            groups.append(gr.update(visible=False))

    rest = sysr["held_back"] + sysr["declined"]
    declined_md = "\n".join(
        f"- **{r['verdict']}**: {short(r)}. " + " ".join(r["why"][:2]) for r in rest
    ) or "Nothing declined."
    base_rows = [[b["GeneSymbol"], b["Name"], b["label"], b["stars"], b["zygosity"], b["VariationID"]]
                 for b in data["baseline"]]
    base_choices = [f"{b['GeneSymbol']} {b['Name'][:70]} | {b['VariationID']}" for b in data["baseline"]]
    return ([data, summary] + htmls + groups
            + [declined_md, gr.update(choices=[short(r) for r in rest], value=None), base_rows,
               gr.update(choices=base_choices, value=None)])


def send(data, session, source, variation_id, gene, verdict, answer, reason):
    if not data:
        return "Load a sample first."
    if answer not in ("agree", "disagree"):
        return "Pick agree or disagree."
    try:
        r = requests.post(f"{BACKEND}/feedback", timeout=60, json={
            "session_id": session, "sample": data["sample"], "variation_id": int(variation_id), "gene": gene,
            "source": source, "shown_verdict": verdict, "answer": answer, "reason": reason or None})
        r.raise_for_status()
    except requests.RequestException as ex:
        return f"Not saved: {ex}"
    return f"Saved ({answer})."


def send_card(i):
    def fn(data, session, answer, reason):
        if not data or i >= len(data["system"]["findings"]):
            return "Load a sample first."
        r = data["system"]["findings"][i]
        return send(data, session, "system", r["evidence"]["variation_id"], r["evidence"]["gene"], r["verdict"],
                    answer, reason)
    return fn


def send_pick(source):
    def fn(data, session, pick, answer, reason):
        if not pick:
            return "Pick a finding from the list."
        vid = int(pick.rsplit("|", 1)[1])
        if source == "baseline":
            b = next(b for b in data["baseline"] if b["VariationID"] == vid)
            return send(data, session, "baseline", vid, b["GeneSymbol"], b["label"], answer, reason)
        r = next(r for r in data["system"]["held_back"] + data["system"]["declined"]
                 if r["evidence"]["variation_id"] == vid)
        return send(data, session, "system", vid, r["evidence"]["gene"], r["verdict"], answer, reason)
    return fn


with gr.Blocks(title="Genome Intelligence") as demo:
    data = gr.State(None)
    session = gr.State(lambda: str(uuid.uuid4()))
    gr.Markdown("# Challengeable findings\nA second opinion on ClinVar 'pathogenic' calls, for one real genome. "
                "Each finding shows the evidence for and against it. Tell us where you disagree.")
    gr.Markdown(NOTE)
    with gr.Row():
        choice = gr.Radio(["Development sample (HG03774)", "Held-out sample"], value="Development sample (HG03774)",
                          label="Bundled 1000 Genomes sample")
        upload = gr.File(label="Or upload a VCF (GRCh38, single sample)", file_types=[".vcf", ".gz"], type="filepath")
    go = gr.Button("Analyse", variant="primary")
    summary = gr.Markdown()

    with gr.Tabs():
        with gr.Tab("Findings (this system)"):
            card_html_out, card_groups = [], []
            for i in range(SLOTS):
                h = gr.HTML(visible=False)
                with gr.Row(visible=False) as g:
                    ans = gr.Radio(["agree", "disagree"], label="Do you agree with this verdict?")
                    why = gr.Textbox(label="Why (optional)", lines=1)
                    btn = gr.Button("Save answer")
                    status = gr.Markdown()
                btn.click(send_card(i), [data, session, ans, why], status)
                card_html_out.append(h)
                card_groups.append(g)
        with gr.Tab("Declined and held back"):
            declined_md = gr.Markdown()
            d_pick = gr.Dropdown(label="Disagree with one of these?", choices=[])
            with gr.Row():
                d_ans = gr.Radio(["agree", "disagree"], label="Agree with the verdict?")
                d_why = gr.Textbox(label="Why (optional)", lines=1)
                d_btn = gr.Button("Save answer")
            d_status = gr.Markdown()
            d_btn.click(send_pick("system"), [data, session, d_pick, d_ans, d_why], d_status)
        with gr.Tab("Naive list (baseline)"):
            gr.Markdown("What an obvious tool reports: every variant whose ClinVar headline says pathogenic, "
                        "ranked by label and review stars. No weighing of evidence.")
            base_table = gr.Dataframe(headers=["Gene", "Variant", "ClinVar label", "Stars", "Zygosity", "ClinVar ID"],
                                      interactive=False, wrap=True)
            b_pick = gr.Dropdown(label="Finding", choices=[])
            with gr.Row():
                b_ans = gr.Radio(["agree", "disagree"], label="Would you report this as pathogenic for this person?")
                b_why = gr.Textbox(label="Why (optional)", lines=1)
                b_btn = gr.Button("Save answer")
            b_status = gr.Markdown()
            b_btn.click(send_pick("baseline"), [data, session, b_pick, b_ans, b_why], b_status)

    go.click(load, [choice, upload],
             [data, summary] + card_html_out + card_groups + [declined_md, d_pick, base_table, b_pick])

if __name__ == "__main__":
    demo.launch()
