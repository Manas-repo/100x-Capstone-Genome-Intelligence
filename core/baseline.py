"""Naive baseline: what an obvious tool does (Promethease-style).

VCF in -> look up each variant in ClinVar -> keep the ones whose headline label says pathogenic
-> rank by label and review stars. No weighing of evidence, no zygosity or inheritance logic,
no population check, no reading of submitter reasoning.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vcf_clinvar import DATA, match_clinvar, read_vcf  # noqa: E402

LABEL_RANK = {"Pathogenic": 0, "Pathogenic/Likely pathogenic": 1, "Likely pathogenic": 2}
STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
}


def run_baseline(vcf_path) -> pd.DataFrame:
    hits = match_clinvar(read_vcf(vcf_path))
    label = hits["ClinicalSignificance"].str.split(";").str[0].str.strip()
    hits = hits[label.isin(LABEL_RANK)].copy()
    hits["label"] = label[hits.index]
    hits["stars"] = hits["ReviewStatus"].map(STARS).fillna(0).astype(int)
    hits["rank_key"] = hits["label"].map(LABEL_RANK)
    hits = hits.sort_values(["rank_key", "stars"], ascending=[True, False])
    cols = ["GeneSymbol", "Name", "label", "stars", "ReviewStatus", "NumberSubmitters", "zygosity",
            "PhenotypeList", "VariationID", "chrom", "pos", "ref", "alt"]
    return hits[cols].reset_index(drop=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "dev"
    out = run_baseline(DATA / "samples" / f"{which}_sample.vcf.gz")
    results = DATA.parent / "results"
    results.mkdir(exist_ok=True)
    out.to_csv(results / f"baseline_{which}.csv", index=False)
    pd.set_option("display.width", 220, "display.max_colwidth", 60)
    print(f"{len(out)} findings reported by the baseline")
    print(out[["GeneSymbol", "Name", "label", "stars", "zygosity", "NumberSubmitters"]].to_string())
