"""Gradio frontend (T12c).

If BACKEND_URL is set, it talks to the FastAPI backend over HTTP (separate deploy). If not, it
calls the same backend functions in-process (one Hugging Face Space, or a plain local run).

Layout: Step 1 is the naive list (what other tools show); Step 2 is this tool's second opinion,
split into Reported findings and Not reported. Every item gets the same answer row: Agree/Disagree,
an optional reason and a Save button. The button is orange while there is something to save, and
gray "Saved" once saved. Changing the answer turns it orange again ("Save change"); a re-save adds
a new row, and the latest row per finding is the answer that counts.

Run locally:  python app.py   (repo root)
"""
import html
import os
import uuid
from pathlib import Path

import gradio as gr
import requests

BACKEND = os.getenv("BACKEND_URL", "").rstrip("/")
CARD_SLOTS = 5  # matches SHOWN in core/adjudicate.py
LIST_SLOTS = 20  # naive list and Not reported tab; anything beyond is summarised in one line
TIMEOUT = 300  # first call can wait on a sleeping backend plus model calls

SAMPLES = {"New sample": "heldout", "Development sample (HG03774)": "dev"}
COLORS = {"Defensible": "#1f7a3a", "Contested": "#a15c00", "Declined": "#666"}
PILL = {"P": "#b42318", "LP": "#b42318", "B": "#1f7a3a", "LB": "#1f7a3a"}  # anything else (VUS...) gray
CARRY = {"yes": ("applies to this person", "#1f7a3a"), "partly": ("partly applies", "#a15c00"),
         "no": ("does not apply to this person", "#b42318")}  # the model's carries_over answer, in plain words
LAB_CAP = 10  # labs listed per side; ClinVar has the rest
CLINVAR = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}/"

INTRO = """# Challengeable findings
A second opinion on ClinVar "pathogenic" calls for one real genome. Instead of a ranked list to accept,
each finding shows the evidence for and against it, so you can judge it.

**Not a diagnosis.** A research prototype on public 1000 Genomes data (no patient data). Not for medical decisions.

### How to use it
1. Press **Start**. It loads a new genome: a South Asian person from the public 1000 Genomes project.
   The first load can take up to a minute.
2. **Step 1: What other tools show.** The usual ranked list. Would you report each one as pathogenic for this person?
3. **Step 2: Our second opinion.** The same genome, after this tool weighed the evidence. It has two parts,
   **Reported findings** and **Not reported**. Tell us where you agree and where you don't.
4. For each item, pick an answer, add a short reason if you like, and press the orange **Save** button.
   It turns gray and says **Saved**. Changed your mind? Pick the other answer and press **Save change**.
   The latest answer is the one that counts.
"""

STEP1_LABEL = '<span class="step-label s1">Typical tool output</span>'
STEP2_LABEL = '<span class="step-label s2">This tool</span>'

TAB_NAIVE = """**What this is:** the kind of list tools like Promethease or a simple ClinVar filter give you. Every
variant in this genome whose ClinVar headline says pathogenic, ranked by label and review stars, with no
weighing of the evidence. We ask you about it first so we can compare it with Step 2.

**What to do:** for each one, would you report it as pathogenic for this person? Go with your first instinct."""

TAB_STEP2 = """**What this is:** the same genome, after this tool checked the evidence behind each ClinVar call. It splits
them into two groups. **Reported findings** are the ones it thinks deserve your attention, each with the evidence
for and against. **Not reported** are the ones it declined. Please look at both."""

TAB_FINDINGS = """**What to do:** agree or disagree with each verdict, and say why if you can. The legend below
explains the colored labels on the cards."""

TAB_REST = """**What to do:** tell us if any of these should have been reported. <b>Declined</b> means the evidence does
not hold up (for example, the variant is too common in South Asians to cause a rare disease). <b>Held back</b>
means it passed, but only the top 5 are shown as cards."""

NOT_LOADED = "_Press **Start** above to load a genome._"

# Save button states are set by class so Gradio's focus styling can't leave a saved button orange
CSS = """
.item-box { border: 1px solid #ccc !important; border-radius: 10px !important; padding: 6px !important;
            margin-bottom: 14px !important; }
.item-box > div > div:first-child > div { border: none !important; }
.save-btn { transition: none !important; }
.save-btn.idle, .save-btn.saved { background: #9ca3af !important; border-color: #9ca3af !important;
                                  color: #fff !important; opacity: .8; }
.save-btn.unsaved { background: #ea580c !important; border-color: #ea580c !important; color: #fff !important; }

.loading-banner { display: flex; align-items: center; gap: 12px; border: 2px solid #ea580c; border-radius: 8px;
                  background: rgba(234, 88, 12, .10); padding: 12px 16px; font-size: 15px; font-weight: 600; }
.spinner { width: 20px; height: 20px; flex: none; border: 3px solid rgba(234, 88, 12, .3);
           border-top-color: #ea580c; border-radius: 50%; animation: gi-spin .8s linear infinite; }
@keyframes gi-spin { to { transform: rotate(360deg); } }
.summary-box { border: 1px solid var(--border-color-primary); border-radius: 8px; padding: 10px 14px; }
.summary-box.error { border-color: #b42318; }

.step-label { display: inline-block; padding: 2px 10px; border-radius: 12px; font-size: 12px; font-weight: 600;
              color: #fff; }
.step-label.s1 { background: #6b7280; }
.step-label.s2 { background: #ea580c; }

.fc { font-size: 14px; line-height: 1.5; padding: 4px 6px; }
.fc-head { font-size: 1.08em; }
.fc-badge { color: #fff; padding: 2px 8px; border-radius: 4px; font-weight: 600; }
.fc-sub { font-size: 13px; margin: 4px 0 8px; }
.fc-chips { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 10px; }
.fc-chip { background: var(--background-fill-secondary); border: 1px solid var(--border-color-primary);
           border-radius: 6px; padding: 3px 8px; font-size: 12.5px; }
.fc-why { background: var(--background-fill-secondary); border-left: 4px solid; padding: 8px 12px;
          margin-bottom: 12px; }
.fc-why ul { margin: 4px 0 0 18px !important; padding: 0 !important; }
.fc-why li { margin: 2px 0 !important; }
.fc-cols { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; }
.fc-col-title { font-weight: 600; font-size: 13px; margin-bottom: 2px; }
.fc-lab { border-top: 1px solid var(--border-color-primary); padding: 6px 0; font-size: 13px; }
.fc-pill { display: inline-block; min-width: 30px; text-align: center; border-radius: 4px; font-size: 11.5px;
           font-weight: 600; padding: 0 5px; margin-right: 6px; color: #fff; }
.fc-meta { color: #4b5563; font-size: 12px; }
.dark .fc-meta { color: #b0b6c0; }
.fc-tag { font-size: 11.5px; font-weight: 600; border: 1px solid; border-radius: 10px; padding: 0 7px;
          margin-left: 4px; white-space: nowrap; }
.fc-judge { font-size: 12.5px; margin-top: 2px; }
.fc-lab details { margin-top: 2px; font-size: 12.5px; }
.fc-lab summary { cursor: pointer; list-style: none; }
.fc-lab summary::-webkit-details-marker { display: none; }
.fc-lab details[open] .teaser { display: none; }
.fc-lab .more { color: #ea580c; font-weight: 600; white-space: nowrap; }
.fc-lab .more::after { content: "show more"; }
.fc-lab details[open] .more::after { content: "show less"; }
.fc-foot { margin-top: 10px; font-size: 13px; }

.legend { border: 1px solid var(--border-color-primary); border-radius: 8px; padding: 10px 14px; margin: 6px 0 14px;
          background: var(--background-fill-secondary); font-size: 13px; }
.legend-title { font-weight: 700; font-size: 14px; margin-bottom: 6px; }
.legend-group { display: grid; gap: 4px; margin-bottom: 8px; }
.legend-sub { font-weight: 600; font-size: 12.5px; color: #4b5563; }
.dark .legend-sub { color: #b0b6c0; }
.legend-row { display: grid; grid-template-columns: 210px 1fr; gap: 10px; align-items: baseline; }
.legend-key .fc-tag { margin-left: 0; }
@media (max-width: 600px) { .legend-row { grid-template-columns: 1fr; gap: 2px; } }

/* Step 2 sub-tabs as clear buttons so Reported findings / Not reported is easy to spot */
.subtabs-hint { font-weight: 600; margin: 8px 0 2px; }
.subtabs > div:first-child { gap: 8px !important; border-bottom: none !important; margin-bottom: 8px; }
.subtabs button[role="tab"] { border: 2px solid #ea580c !important; border-radius: 8px !important;
                              padding: 8px 18px !important; font-size: 15px !important; font-weight: 600 !important;
                              color: #ea580c !important; background: transparent !important; }
.subtabs button[role="tab"][aria-selected="true"] { background: #ea580c !important; color: #fff !important; }
.subtabs button[role="tab"]::after { display: none !important; }
"""


def legend_html() -> str:
    """The key to every colored label on a card, built from the same colors the cards use."""
    def badge(v):
        return f'<span class="fc-badge" style="background:{COLORS[v]}">{v}</span>'

    def tag(key):
        text, color = CARRY[key]
        return f'<span class="fc-tag" style="color:{color};border-color:{color}">{text}</span>'

    def pill(call, color):
        return f'<span class="fc-pill" style="background:{color}">{call}</span>'

    rows = [
        ("Verdict", [
            (badge("Defensible"), "the evidence holds up for this person."),
            (badge("Contested"), "labs disagree, or the evidence is weak or does not apply to this person."),
        ]),
        ("Does the lab's reasoning apply to this person?", [
            (tag("yes"), "it is about the variant itself (experiments, rarity, predictions), so it holds for anyone."),
            (tag("partly"), "some of it is about the variant, some about the lab's own patient."),
            (tag("no"), "it is about the lab's own patient (their symptoms, their family), not this person."),
        ]),
        ("Each lab's call", [
            (pill("P", "#b42318") + pill("LP", "#b42318"), "pathogenic, likely pathogenic."),
            (pill("B", "#1f7a3a") + pill("LB", "#1f7a3a"), "benign, likely benign."),
            (pill("VUS", "#6b7280"), "uncertain significance (or another call)."),
        ]),
    ]
    out = ['<div class="legend"><div class="legend-title">Legend</div>']
    for title, items in rows:
        out.append(f'<div class="legend-group"><div class="legend-sub">{title}</div>')
        out += [f'<div class="legend-row"><span class="legend-key">{k}</span><span>{d}</span></div>'
                for k, d in items]
        out.append('</div>')
    out.append('<div class="fc-meta">"This person" is an unrelated, healthy adult. An AI model reads each lab&#39;s '
               'written reasoning and judges whether it applies.</div></div>')
    return "".join(out)


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def pct(x: float) -> str:
    return f"{x:.2%}" if x >= 0.0001 else f"{x:.4%}"


def _lab_row(s: dict, model: dict | None) -> str:
    """One lab as a compact row: call pill, year, lab, the model's judgement in view, the lab's own words folded."""
    call = s["call"] or "?"
    pill = f'<span class="fc-pill" style="background:{PILL.get(call, "#6b7280")}">{esc(call)}</span>'
    head = (f'{pill}<b>{esc(s["year"] or "no date")}</b> · {esc(s["submitter"])} '
            f'<span class="fc-meta">({esc(s["method"])})</span>')
    tag_text, tag_color = CARRY.get((model or {}).get("carries_over"), ("", ""))
    judge = ""
    if tag_text:
        head += f' <span class="fc-tag" style="color:{tag_color};border-color:{tag_color}">{tag_text}</span>'
        if (model or {}).get("reason"):
            judge = f'<div class="fc-judge"><b>Why:</b> {esc(model["reason"])}</div>'
    comment = s.get("comment")
    if not comment:
        return f'<div class="fc-lab">{head}<div class="fc-meta">No reasoning given.</div></div>'
    if len(comment) > 700:
        comment = comment[:700] + "..."
    if len(comment) <= 140:
        words = f'<div class="fc-meta"><b>Lab says:</b> <i>{esc(comment)}</i></div>'
    else:
        teaser = comment[:140].rsplit(" ", 1)[0] + "..."
        words = (f'<details class="fc-meta"><summary><b>Lab says:</b> <i class="teaser">{esc(teaser)}</i> '
                 f'<span class="more"></span></summary><i>{esc(comment)}</i></details>')
    return f'<div class="fc-lab">{head}{judge}{words}</div>'


def card_html(r: dict) -> str:
    e, v = r["evidence"], r["verdict"]
    review = r.get("review") or {}
    by_sub = {l.get("submitter"): l for l in review.get("labs", []) if l.get("submitter")}
    pop = e["population"]
    years = [s["year"] for s in e["support"] + e["dispute"] if s["year"]]
    age = f"{min(years)} to {max(years)}" if years else "no dates given"
    if pop.get("found"):
        fold = pop["af_sas"] / pop["af_all"] if pop["af_all"] else 0
        ratio = f" ({fold:.0f}x)" if fold >= 2 else ""
        popline = (f'South Asian <b>{pct(pop["af_sas"])}</b> vs all <b>{pct(pop["af_all"])}</b>{ratio}'
                   f' <span class="fc-meta">gnomAD, {pop["an_sas"]:,} South Asian alleles</span>')
    else:
        popline = "No South Asian frequency in gnomAD"
    cond = e["conditions"] or "not listed"
    if len(cond) > 120:
        cond = cond[:120] + "..."
    def labs(rows: list, cap: int) -> str:
        out = "".join(_lab_row(s, by_sub.get(s["submitter"])) for s in rows[:cap])
        if len(rows) > cap:
            out += (f'<div class="fc-lab fc-meta">Showing {cap} of {len(rows)} labs. '
                    f'The other {len(rows) - cap} are on the ClinVar record.</div>')
        return out or '<div class="fc-lab fc-meta">None</div>'

    sup, dis = labs(e["support"], LAB_CAP), labs(e["dispute"], LAB_CAP)
    why = "".join(f"<li>{esc(w)}</li>" for w in r["why"])
    check = review.get("cheapest_check")
    if check and ("gnomad" in check.lower() or "frequency" in check.lower()):
        check = None  # the card already shows population frequency
    check_html = f"<b>Quickest check:</b> {esc(check)} &nbsp;·&nbsp; " if check else ""
    return f"""
<div class="fc">
  <div class="fc-head"><span class="fc-badge" style="background:{COLORS[v]}">{v}</span>
    &nbsp;<b>{esc(e['gene'])}</b> {esc(e['name'])}</div>
  <div class="fc-sub">{esc(e['consequence'])}, {esc(e['zygosity'])}. {esc(e['finding_type'])}.</div>
  <div class="fc-chips">
    <span class="fc-chip">Condition: <b>{esc(cond)}</b></span>
    <span class="fc-chip">ClinVar: <b>{esc(e['clinvar_label'])}</b> <span class="fc-meta">({esc(e['review_status'])})</span></span>
    <span class="fc-chip">Evidence: <b>{age}</b></span>
    <span class="fc-chip">{popline}</span>
  </div>
  <div class="fc-why" style="border-left-color:{COLORS[v]}"><b>Why {v}</b><ul>{why}</ul></div>
  <div class="fc-cols">
    <div><div class="fc-col-title">For pathogenic ({len(e['support'])})</div>{sup}</div>
    <div><div class="fc-col-title">Against ({len(e['dispute'])})</div>{dis}</div>
  </div>
  <div class="fc-foot">{check_html}<a href="{e['clinvar_url']}" target="_blank">Open the ClinVar record</a></div>
</div>"""


def naive_html(i: int, b: dict) -> str:
    stars = "&#9733;" * b["stars"] + "&#9734;" * (4 - b["stars"])
    names = (b.get("PhenotypeList") or "").replace(";", "|").split("|")
    cond = ", ".join(dict.fromkeys(n.strip() for n in names if n.strip() not in ("", "not provided", "not specified")))
    if len(cond) > 160:
        cond = cond[:160] + "..."
    return f"""
<div style="border-radius:8px;padding:10px 14px">
  <b>#{i + 1}&nbsp; {esc(b['GeneSymbol'])}</b> {esc(b['Name'])}<br>
  ClinVar: <b>{esc(b['label'])}</b> &nbsp;<span title="{esc(b['ReviewStatus'])}">{stars}</span>
  &nbsp;| {esc(b['zygosity'])} &nbsp;| {esc(cond) or 'no condition listed'}
  &nbsp;<a href="{CLINVAR.format(b['VariationID'])}" target="_blank">ClinVar record</a>
</div>"""


def rest_html(r: dict, held: bool) -> str:
    e = r["evidence"]
    label, color = ("Held back", COLORS["Contested"]) if held else ("Declined", COLORS["Declined"])
    why = " ".join(esc(w) for w in r["why"][:2])
    return f"""
<div style="border-left:6px solid {color};border-radius:8px;padding:10px 14px">
  <span style="background:{color};color:#fff;padding:1px 7px;border-radius:4px">{label}</span>
  &nbsp;<b>{esc(e['gene'])}</b> {esc(e['name'])}<br>
  {f'Verdict: {esc(r["verdict"])}. ' if held else ''}{why}
  &nbsp;<a href="{e['clinvar_url']}" target="_blank">ClinVar record</a>
</div>"""


# ---- what each slot points at -------------------------------------------------------------

def _rest(data: dict) -> list[tuple[dict, bool]]:
    return ([(r, True) for r in data["system"]["held_back"]]
            + [(r, False) for r in data["system"]["declined"]])


def item_for(data, kind: str, i: int):
    """(source, variation_id, gene, shown_verdict) for slot i of a tab, or None."""
    if not data:
        return None
    if kind == "naive":
        rows = data["baseline"]
        if i < len(rows):
            b = rows[i]
            return "baseline", b["VariationID"], b["GeneSymbol"], b["label"]
        return None
    rows = data["system"]["findings"] if kind == "card" else [r for r, _ in _rest(data)]
    if i < len(rows):
        r = rows[i]
        return "system", r["evidence"]["variation_id"], r["evidence"]["gene"], r["verdict"]
    return None


# ---- loading a genome ---------------------------------------------------------------------

def fetch(sample: str | None, upload: str | None) -> dict:
    if not BACKEND:
        from api.main import analyse_file, analyse_named

        if upload:
            return analyse_file(Path(upload), os.path.basename(upload))
        return analyse_named(sample)
    if upload:
        with open(upload, "rb") as f:
            resp = requests.post(f"{BACKEND}/analyse", files={"file": (os.path.basename(upload), f)},
                                 timeout=TIMEOUT)
    else:
        resp = requests.post(f"{BACKEND}/analyse/{sample}", timeout=TIMEOUT)
    resp.raise_for_status()
    return resp.json()


def _slot(html_value: str | None):
    """Updates for one slot: html, box, answer, reason, save button, saved flag."""
    show = html_value is not None
    return [gr.update(value=html_value or ""), gr.update(visible=show),
            gr.update(value=None), gr.update(value=""), btn_idle(), False]


def btn(value: str, state: str):
    return gr.update(value=value, interactive=state == "unsaved", elem_classes=["save-btn", state])


def btn_idle():
    return btn("Save", "idle")


def will_load(label: str) -> str:
    name = "new sample" if SAMPLES[label] == "heldout" else "development sample (HG03774)"
    return f"Start loads the **{name}**. Uploading your own VCF is optional (single sample, GRCh38)."


def loading():
    """Shown the moment Start is pressed: spinner banner, Start and Upload locked until the load ends."""
    banner = ('<div class="loading-banner"><span class="spinner"></span>'
              'Loading the genome. This can take up to a minute the first time.</div>')
    return banner, gr.update(value="Loading...", interactive=False), gr.update(interactive=False)


def ready():
    return [gr.update(value="Start", interactive=True), gr.update(interactive=True)]


def load(sample_label: str, upload: str | None = None):
    try:
        data = fetch(SAMPLES[sample_label], upload)
    except Exception as ex:  # bad upload, backend down: show it instead of crashing the page
        msg = getattr(getattr(ex, "response", None), "text", "") or f"{type(ex).__name__}: {ex}"
        out = [None, f'<div class="summary-box error"><b>Could not load this genome.</b> {esc(msg)}</div>',
               NOT_LOADED, NOT_LOADED, NOT_LOADED]
        for _ in range(CARD_SLOTS + 2 * LIST_SLOTS):
            out += _slot(None)
        return out + ready()

    sysr = data["system"]
    cards, rest, base = sysr["findings"], _rest(data), data["baseline"]
    name = {"heldout": "new sample", "dev": "development sample"}.get(data["sample"], data["sample"])
    summary = (f"<b>Loaded: {esc(name)}.</b> {sysr['n_variants']:,} of this person's variants are in ClinVar, "
               f"in the disease genes this tool checks. Other tools would report <b>{len(base)}</b>. "
               f"This tool reports <b>{len(cards)}</b> and does not report <b>{len(rest)}</b>. "
               f"Begin with <b>Step 1</b> below.")
    if not sysr["model_used"]:
        summary += "<br><i>Note: the model step is off, so this ran on rules only.</i>"
    summary = f'<div class="summary-box">{summary}</div>'

    def count(n: int, noun: str, empty: str) -> str:
        if not n:
            return empty
        extra = f" Showing the first {LIST_SLOTS}." if n > LIST_SLOTS else ""
        return f"_{n} {noun}{'' if n == 1 else 's'} below.{extra}_"

    out = [data, summary,
           count(len(base), "item", "_The naive list is empty for this genome._"),
           count(len(cards), "finding", "_No finding survived. For a healthy person that is often the "
                                        "correct answer._"),
           count(len(rest), "item", "_Nothing was declined or held back._")]
    for i in range(LIST_SLOTS):
        out += _slot(naive_html(i, base[i]) if i < len(base) else None)
    for i in range(CARD_SLOTS):
        out += _slot(card_html(cards[i]) if i < len(cards) else None)
    for i in range(LIST_SLOTS):
        out += _slot(rest_html(*rest[i]) if i < len(rest) else None)
    return out + ready()


# ---- saving answers -----------------------------------------------------------------------

def on_edit(answer, saved):
    """Orange while there is an unsaved answer; 'Save change' once a previous answer was saved."""
    if answer not in ("agree", "disagree"):
        return btn_idle()
    return btn("Save change" if saved else "Save", "unsaved")


def saver(kind: str, i: int):
    def fn(data, session, answer, reason):
        item = item_for(data, kind, i)
        if not item:
            raise gr.Error("Load a genome first.")
        if answer not in ("agree", "disagree"):
            raise gr.Error("Pick an answer first.")
        source, vid, gene, shown = item
        row = {"session_id": session, "sample": data["sample"], "variation_id": int(vid), "gene": gene,
               "source": source, "shown_verdict": shown, "answer": answer, "reason": (reason or "").strip() or None}
        try:
            if BACKEND:
                requests.post(f"{BACKEND}/feedback", timeout=60, json=row).raise_for_status()
            else:
                from api.main import save_feedback

                save_feedback(row)
        except Exception as ex:  # button stays orange so the answer can be saved again
            raise gr.Error(f"Not saved, please try again. ({type(ex).__name__})")
        return btn("Saved ✓", "saved"), True
    return fn


def slots(kind: str, n: int, question: str, choices: list[tuple[str, str]], session, data) -> list:
    comps = []
    for i in range(n):
        with gr.Column(visible=False, elem_classes="item-box") as box:
            h = gr.HTML()
            with gr.Row(equal_height=True):
                ans = gr.Radio(choices, label=question, scale=2)
                why = gr.Textbox(label="Why (optional)", lines=1, scale=3)
                save = gr.Button("Save", interactive=False, scale=0, min_width=130,
                                 elem_classes=["save-btn", "idle"])
        saved = gr.State(False)
        ans.input(on_edit, [ans, saved], save)
        why.input(on_edit, [ans, saved], save)
        save.click(saver(kind, i), [data, session, ans, why], [save, saved])
        comps += [h, box, ans, why, save, saved]
    return comps


with gr.Blocks(title="Genome Intelligence") as demo:
    data = gr.State(None)
    session = gr.State(lambda: str(uuid.uuid4()))
    gr.Markdown(INTRO)

    with gr.Row(equal_height=True):
        start = gr.Button("Start", variant="primary", size="md", scale=0, min_width=140)
        upload = gr.UploadButton("Upload your own VCF", file_types=[".vcf", ".gz"], size="md",
                                 variant="secondary", scale=0, min_width=200)
        target = gr.Markdown(will_load(next(iter(SAMPLES))))
    with gr.Accordion("Builder options", open=False):
        sample = gr.Radio(list(SAMPLES), value=next(iter(SAMPLES)), label="Sample that Start loads")
    sample.change(will_load, sample, target)
    summary = gr.HTML()

    with gr.Tabs():
        with gr.Tab("Step 1: What other tools show", render_children=True):
            gr.HTML(STEP1_LABEL)
            gr.Markdown(TAB_NAIVE)
            naive_count = gr.Markdown(NOT_LOADED)
            naive = slots("naive", LIST_SLOTS, "Report as pathogenic for this person?",
                          [("Yes", "agree"), ("No", "disagree")], session, data)
        with gr.Tab("Step 2: Our second opinion", render_children=True):
            gr.HTML(STEP2_LABEL)
            gr.Markdown(TAB_STEP2)
            gr.HTML('<div class="subtabs-hint">Switch between the two groups here:</div>')
            with gr.Tabs(elem_classes="subtabs"):
                with gr.Tab("Reported findings", render_children=True):
                    gr.Markdown(TAB_FINDINGS)
                    gr.HTML(legend_html())
                    card_count = gr.Markdown(NOT_LOADED)
                    cards = slots("card", CARD_SLOTS, "Do you agree with this verdict?",
                                  [("Agree", "agree"), ("Disagree", "disagree")], session, data)
                with gr.Tab("Not reported", render_children=True):
                    gr.Markdown(TAB_REST)
                    rest_count = gr.Markdown(NOT_LOADED)
                    rest = slots("rest", LIST_SLOTS, "Agree it should not be reported?",
                                 [("Agree", "agree"), ("Disagree", "disagree")], session, data)

    outputs = ([data, summary, naive_count, card_count, rest_count] + naive + cards + rest
               + [start, upload])
    start.click(loading, None, [summary, start, upload]).then(load, sample, outputs, show_progress="hidden")
    upload.upload(loading, None, [summary, start, upload]).then(load, [sample, upload], outputs,
                                                                show_progress="hidden")

if __name__ == "__main__":
    demo.launch(css=CSS)
