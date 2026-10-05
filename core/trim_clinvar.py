"""Trim ClinVar tab-delimited files down to panel genes (GRCh38).

Inputs in data/raw/ (from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/):
  variant_summary[_TAG].txt.gz     one row per variant, headline classification, VCF-style alleles
  submission_summary[_TAG].txt.gz  one row per lab submission (SCV): call, date, method, comment

Outputs in data/:
  clinvar_variants[_TAG].parquet
  clinvar_submissions[_TAG].parquet

Usage: python core/trim_clinvar.py            (current release)
       python core/trim_clinvar.py 2020-01    (archived release, for the backtest)
"""
import gzip
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
DATA = ROOT / "data"

VARIANT_COLS = [
    "VariationID", "AlleleID", "Type", "Name", "GeneSymbol", "ClinicalSignificance", "ClinSigSimple",
    "LastEvaluated", "RS# (dbSNP)", "PhenotypeList", "Assembly", "Chromosome", "ReviewStatus",
    "NumberSubmitters", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
]


def in_panel(df: pd.DataFrame, panel: pd.DataFrame) -> pd.Series:
    keep = pd.Series(False, index=df.index)
    chrom = "chr" + df["Chromosome"].astype(str)
    for p in panel.itertuples():
        keep |= (chrom == p.chrom) & df["PositionVCF"].between(p.start, p.end)
    return keep


def main(tag: str = "") -> None:
    suffix = f"_{tag}" if tag else ""
    panel = pd.read_csv(DATA / "panel_genes.csv")

    vpath = DATA / f"clinvar_variants{suffix}.parquet"
    if vpath.exists():
        variants = pd.read_parquet(vpath)
        print(f"variants in panel (cached): {len(variants):,}")
    else:
        variants = trim_variants(suffix, panel)
        variants.to_parquet(vpath, index=False)
        print(f"variants in panel: {len(variants):,}")
    trim_submissions(suffix, set(variants["VariationID"]))


def trim_variants(suffix: str, panel: pd.DataFrame) -> pd.DataFrame:
    if suffix:
        return trim_archived_variants(suffix)
    parts = []
    for chunk in pd.read_csv(
        RAW / f"variant_summary{suffix}.txt.gz", sep="\t", chunksize=500_000, low_memory=False,
        usecols=lambda c: c.lstrip("#") in VARIANT_COLS,
    ):
        chunk.columns = [c.lstrip("#") for c in chunk.columns]
        chunk = chunk[chunk["Assembly"] == "GRCh38"]
        chunk = chunk[pd.to_numeric(chunk["PositionVCF"], errors="coerce").fillna(-1) > 0]
        chunk["PositionVCF"] = chunk["PositionVCF"].astype(int)
        parts.append(chunk[in_panel(chunk, panel)])
    return pd.concat(parts).drop_duplicates("VariationID")


def trim_archived_variants(suffix: str) -> pd.DataFrame:
    """Older releases have no VCF-style position columns, so keep variants by VariationID
    using the current release's panel set (run the current release first)."""
    current_ids = set(pd.read_parquet(DATA / "clinvar_variants.parquet", columns=["VariationID"])["VariationID"])
    parts = []
    for chunk in pd.read_csv(RAW / f"variant_summary{suffix}.txt.gz", sep="\t", chunksize=500_000, low_memory=False):
        chunk.columns = [c.lstrip("#") for c in chunk.columns]
        chunk = chunk[(chunk["Assembly"] == "GRCh38") & chunk["VariationID"].isin(current_ids)]
        parts.append(chunk[[c for c in VARIANT_COLS if c in chunk.columns]])
    return pd.concat(parts).drop_duplicates("VariationID")


def trim_submissions(suffix: str, ids: set) -> None:
    # The submission file has comment lines before the header; find the header line.
    path = RAW / f"submission_summary{suffix}.txt.gz"
    with gzip.open(path, "rt", encoding="utf8", errors="replace") as f:
        skip = 0
        for line in f:
            if line.startswith("#VariationID\t"):
                break
            skip += 1
    parts = []
    for chunk in pd.read_csv(
        path, sep="\t", skiprows=skip, chunksize=500_000, low_memory=False, dtype=str,
        on_bad_lines="skip", quoting=3,
    ):
        chunk.columns = [c.lstrip("#") for c in chunk.columns]
        chunk = chunk[chunk["VariationID"].astype(int).isin(ids)]
        parts.append(chunk)
    subs = pd.concat(parts)
    subs.to_parquet(DATA / f"clinvar_submissions{suffix}.parquet", index=False)
    print(f"submissions for those variants: {len(subs):,}")
    print("submission columns:", list(subs.columns))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
