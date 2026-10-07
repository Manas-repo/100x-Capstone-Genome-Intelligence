---
title: Genome Intelligence
emoji: 🧬
colorFrom: green
colorTo: gray
sdk: gradio
sdk_version: 6.29.1
app_file: app.py
python_version: "3.12"
short_description: Challengeable findings from ClinVar P/LP calls
startup_duration_timeout: 30m
pinned: false
---

# Genome Intelligence: a second opinion on ClinVar "pathogenic" calls

Takes a real variant file (VCF), finds the variants that ClinVar calls pathogenic or likely pathogenic, and checks whether the reasoning behind each label holds up **for this person**. It shows a handful of challengeable findings and openly declines to call the rest.

Each finding shows what supports it, what disputes it, how old that evidence is, and which population it was measured in. The verdict is one of:

- **Defensible**: the evidence is strong and carries over to this person.
- **Contested**: something in the evidence should make a curator stop and check (named on the card).
- **Declined**: the tool will not call it, and says why.

Sense-making, not diagnosis. Not a substitute for a clinician or genetic counselor.

Built for the 100X Engineers C7 Capstone (Genome Intelligence brief).

- **Live app:** https://manasnapte-genome-intelligence-100xcapstone.hf.space
- **Test user:** a Genome Analyst at a Bangalore genomics company (germline exome and carrier screening, ACMG/AMP classification).

## Why it exists

ClinVar collects labs' classifications, it does not judge between them. Most entries rest on one lab, many conflict, and the population evidence is mostly European. A tool that reads the headline label and ranks it (like Promethease, or the "naive list" in this app) has quietly picked a side.

The checks in this tool come from an interview with a practising curator. He named the error he sees most: rare null variants (nonsense, frameshift, splice) called "likely pathogenic" by labs and tools without checking whether loss of function is actually how that gene causes disease. He also explained how he weighs a ClinVar conflict: how many labs, which criteria they used, whether those criteria apply to his patient, and how old the evidence is.

## How it works

```
VCF  ->  match to ClinVar (593 genes)  ->  evidence per variant  ->  rule checks  ->  model step  ->  verdict
```

1. **Gene panel (593 genes).** The 81 ACMG secondary findings genes (v3.2, health risks) plus 512 recessive carrier-screening genes (ClinGen autosomal recessive, Definitive or Strong). Each gene carries its disease mechanism (loss of function or not).
2. **Candidates.** Every variant the person carries where at least one lab says pathogenic or likely pathogenic.
3. **Evidence.** Every lab's call, year, method, review status and written reasoning from ClinVar, plus South Asian and overall frequency from gnomAD v4.
4. **Rule checks (no model)** in [core/evidence.py](core/evidence.py): one lab only, labs in conflict, newest evidence over 5 years old, literature-only, no reasoning given, null variant in a gene whose disease is not loss of function, too common in South Asians to cause the disease, much more common in South Asians than overall, a healthy adult with two copies of a recessive variant, ACMG reporting restrictions, expert panel review.
5. **Model step** in [core/adjudicate.py](core/adjudicate.py) (Groq, `openai/gpt-oss-120b`, temperature 0). For each lab's free-text reasoning, the model says whether it carries over to this person (an unrelated healthy adult with no phenotype or family data). "Found in trans in our patient" or "segregates in this family" does not carry over. A functional study usually does. This is the part a rule cannot do: the reasoning is unstructured text.
6. **Verdict** is plain code over the flags and the model's answers, so every verdict can be traced to named reasons.

The app shows two tabs side by side: **Step 1, what other tools show** (the naive list) and **Step 2, our second opinion**. The tester marks agree or disagree on each item, with a reason, and the answers are saved to Supabase.

## Data

| Source | Use |
|---|---|
| [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) `variant_summary` and `submission_summary`, Oct 2026 and Jan 2020 | Labels, per-lab calls and reasoning; the 2020 copy for the backtest |
| [1000 Genomes high-coverage GRCh38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/) | Two real South Asian genomes, picked with random seed `20261007`. Dev: HG03774 (Indian Telugu). Held-out: not opened by the builder before the user test |
| [gnomAD v4](https://gnomad.broadinstitute.org/) API | South Asian vs overall allele frequency, cached |
| [ACMG SF v3.2](https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/), [ClinGen gene curation and gene-disease validity](https://search.clinicalgenome.org/) | Gene panel, inheritance, mechanism |

Large files (ClinVar SQLite, the two sample VCFs) are in the GitHub release [`data-v1`](https://github.com/Manas-repo/100x-Capstone-Genome-Intelligence/releases/tag/data-v1) and downloaded on first start. No patient data is used anywhere.

## Evaluation

- **Baseline first.** [core/baseline.py](core/baseline.py) is the naive approach: ClinVar's own pathogenic filter, ranked, no weighing. Output on the dev genome: [results/baseline_dev.csv](results/baseline_dev.csv).
- **Backtest against ClinVar's own history** ([core/backtest.py](core/backtest.py)). Take every Jan 2020 P/LP call in the 593 genes (45,999 variants), run the rule checks on 2020 data only, and see which calls ClinVar later walked back. 7.2% were walked back. A call we marked Contested in 2020 was about 4x as likely to be walked back as one we did not (9.5% vs 2.2%). Expert-panel calls were almost never walked back (0.1%). Results in [results/](results/).
- **User test** on the held-out genome with a practising curator, run end to end on the live app. Scheduled for the week of Oct 7, 2026; results will be added here when it is done. Report built with [eval/user_test_report.py](eval/user_test_report.py).
- **Logic freeze.** The verdict logic, rules and model prompt were frozen at git tag `logic-freeze` before the user test. Check with `git diff logic-freeze -- core/`.

## Run locally

```bash
python -m venv .venv
.venv/Scripts/activate        # Windows; use .venv/bin/activate on macOS/Linux
pip install -r requirements.txt
python app.py                 # downloads data on first run, then open http://localhost:7860
```

Optional `.env` (never committed): `GROQ_API_KEY` (model step; without it, verdicts use the rules only), `SUPABASE_URL` and `SUPABASE_SERVICE_ROLE_KEY` (saving answers; without them, answers go to `results/feedback_local.jsonl`). The Supabase table is in [api/supabase_feedback.sql](api/supabase_feedback.sql).

The same functions are also exposed over HTTP with FastAPI: `uvicorn api.main:app`.

## Stack

Gradio UI and backend in one Hugging Face Space (FastAPI kept for a separate deploy), Groq for the model step, Supabase (Postgres) for the tester's answers.

## Limits

- One gene set (germline monogenic), one ancestry group, two genomes, one tester.
- gnomAD's South Asian group is not India-specific, and founder variants in specific Indian communities can look "too common" or "enriched" for the wrong reason.
- The 1000 Genomes phased callset drops variants seen in only one person, so some rare pathogenic variants can be missed.
- The gene mechanism list was built from public sources and has not yet been checked by a curator.
- ClinVar is a snapshot (Oct 2026). If a label changes later, findings already shown are not updated.
