"""Read a single-sample VCF and match its variants to ClinVar (panel genes, GRCh38)."""
import gzip
from functools import lru_cache
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"


def read_vcf(path) -> pd.DataFrame:
    """One row per (site, alt allele) carried by the sample, with zygosity."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alts, qual, filt, info, fmt, sample = line.rstrip("\n").split("\t")[:10]
            gt = sample.split(":")[fmt.split(":").index("GT")] if "GT" in fmt else sample.split(":")[0]
            alleles = gt.replace("|", "/").split("/")
            info_d = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            for i, alt in enumerate(alts.split(","), start=1):
                n = alleles.count(str(i))
                if n == 0 or alt.startswith("<"):  # skip symbolic structural variants
                    continue
                rows.append({
                    "chrom": chrom.removeprefix("chr"), "pos": int(pos), "ref": ref, "alt": alt,
                    "zygosity": "homozygous" if n == 2 else ("hemizygous" if len(alleles) == 1 else "heterozygous"),
                    "filter": filt,
                    "af_1kg_all": _f(info_d.get("AF")),
                    "af_1kg_sas": _f(info_d.get("AF_SAS")),
                })
    return pd.DataFrame(rows)


def _f(x):
    try:
        return float(x.split(",")[0])
    except (AttributeError, ValueError):
        return None


@lru_cache(maxsize=2)
def clinvar_variants(tag: str = "") -> pd.DataFrame:
    suffix = f"_{tag}" if tag else ""
    df = pd.read_parquet(DATA / f"clinvar_variants{suffix}.parquet")
    if "PositionVCF" in df:
        df = df.rename(columns={
            "Chromosome": "chrom", "PositionVCF": "pos", "ReferenceAlleleVCF": "ref", "AlternateAlleleVCF": "alt",
        })
        df["chrom"] = df["chrom"].astype(str)
    return df


@lru_cache(maxsize=2)
def clinvar_submissions(tag: str = "") -> pd.DataFrame:
    suffix = f"_{tag}" if tag else ""
    df = pd.read_parquet(DATA / f"clinvar_submissions{suffix}.parquet")
    df["VariationID"] = df["VariationID"].astype(int)
    return df


def match_clinvar(sample: pd.DataFrame) -> pd.DataFrame:
    """Sample variants that have a ClinVar record, with the ClinVar headline fields."""
    cv = clinvar_variants()
    return sample.merge(cv, on=["chrom", "pos", "ref", "alt"], how="inner")


def panel() -> pd.DataFrame:
    return pd.read_csv(DATA / "panel_genes.csv")
