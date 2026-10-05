"""Backtest (T14): would our rules have warned about pathogenic calls that ClinVar later walked back?

Setup
- Universe: panel variants whose ClinVar headline label in January 2020 was P, LP or P/LP.
  This is exactly what the naive baseline reports, with no warning attached.
- Our system, run on the January 2020 data only: variant-level rule checks (single lab, conflict,
  stale, literature only, no reasoning given, null variant in a non-LoF gene, expert panel).
  Person-specific checks (zygosity, population) and the model step are left out, because there is
  no person here and model calls over thousands of variants are not needed to test the rules.
- Outcome, from today's ClinVar (Oct 2026): is the headline label still P/LP, or was it walked back
  (now VUS, conflicting, or benign)?

Question: among calls later walked back, how many did we flag as Contested in 2020 (catch rate),
and among calls that held up, how many did we flag anyway (false alarm rate)?
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adjudicate import verdict  # noqa: E402
from evidence import assemble, flags  # noqa: E402
from vcf_clinvar import clinvar_variants  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"
PLP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}


def outcome(label: str) -> str:
    first = str(label).split(";")[0].strip()
    if first in PLP:
        return "held up"
    if first.startswith("Conflicting") or first == "Uncertain significance":
        return "walked back (now uncertain/conflicting)"
    if "enign" in first:
        return "walked back (now benign)"
    return "other"


def main() -> None:
    old = clinvar_variants("2020-01")
    now = clinvar_variants().set_index("VariationID")["ClinicalSignificance"]
    old = old[old["ClinicalSignificance"].str.split(";").str[0].str.strip().isin(PLP)]
    old = old[old["VariationID"].isin(now.index)]
    print(f"variants the baseline would report as pathogenic in Jan 2020: {len(old):,}")

    rows = []
    for i, (_, r) in enumerate(old.iterrows()):
        ev = assemble(r, tag="2020-01", with_population=False)
        fl = [f for f in flags(ev, as_of_year=2020, person=False) if f["effect"] != "decline"]
        v, why = verdict(ev, fl, None)
        rows.append({
            "VariationID": r["VariationID"], "gene": ev["gene"], "name": ev["name"],
            "label_2020": r["ClinicalSignificance"], "label_now": now[r["VariationID"]],
            "outcome": outcome(now[r["VariationID"]]), "system_2020": v,
            "flags": ", ".join(f["code"] for f in fl),
        })
        if i % 2000 == 0:
            print(f"  {i:,} done")
    df = pd.DataFrame(rows)
    df = df[df["outcome"] != "other"]
    df["walked_back"] = df["outcome"].str.startswith("walked back")
    df["flagged"] = df["system_2020"] == "Contested"
    RESULTS.mkdir(exist_ok=True)
    df.to_csv(RESULTS / "backtest_rows.csv", index=False)

    wb, held = df[df["walked_back"]], df[~df["walked_back"]]
    summary = pd.DataFrame([
        {"group": "walked back since 2020", "n": len(wb), "flagged_contested_2020": wb["flagged"].mean()},
        {"group": "held up", "n": len(held), "flagged_contested_2020": held["flagged"].mean()},
    ])
    # Which single flags carry the signal?
    per_flag = []
    for code in ["single_lab", "conflict", "stale", "literature_only", "no_reasoning", "null_not_lof", "expert_panel"]:
        has = df["flags"].str.contains(code)
        if has.sum():
            per_flag.append({
                "flag": code, "n_flagged": int(has.sum()),
                "walked_back_rate_if_flag": df.loc[has, "walked_back"].mean(),
                "walked_back_rate_if_not": df.loc[~has, "walked_back"].mean(),
            })
    per_flag = pd.DataFrame(per_flag)
    summary.to_csv(RESULTS / "backtest_summary.csv", index=False)
    per_flag.to_csv(RESULTS / "backtest_per_flag.csv", index=False)
    pd.set_option("display.width", 200)
    print("\nBaseline flags 0% of these: it reports every one as pathogenic.\n")
    print(summary.to_string(index=False))
    print()
    print(per_flag.to_string(index=False))
    print(f"\noverall walked-back rate: {df['walked_back'].mean():.1%} of {len(df):,}")


if __name__ == "__main__":
    main()
