"""Read a single-sample VCF and match its variants to ClinVar (panel genes, GRCh38)."""
import gzip
import sqlite3
from functools import lru_cache
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
# The app reads current ClinVar from an indexed SQLite file (built by build_app_db below) so it fits
# in a small server's memory. The parquet files are only needed to build it and for the backtest.
APP_DB = DATA / "clinvar_app.sqlite"
VAR_COLS = ["chrom", "pos", "ref", "alt", "VariationID", "Name", "GeneSymbol", "ClinicalSignificance",
            "ReviewStatus", "NumberSubmitters", "PhenotypeList"]
SUB_COLS = ["VariationID", "ClinicalSignificance", "DateLastEvaluated", "Description", "ReportedPhenotypeInfo",
            "ReviewStatus", "CollectionMethod", "Submitter", "SCV"]


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


def submissions_for(vid: int, tag: str = "") -> pd.DataFrame:
    """All lab submissions for one ClinVar variant."""
    if not tag and APP_DB.exists():
        with sqlite3.connect(APP_DB) as con:
            return pd.read_sql("SELECT * FROM submissions WHERE VariationID = ?", con, params=(int(vid),))
    subs = clinvar_submissions(tag)
    return subs[subs["VariationID"] == int(vid)].copy()


def match_clinvar(sample: pd.DataFrame) -> pd.DataFrame:
    """Sample variants that have a ClinVar record, with the ClinVar headline fields."""
    if APP_DB.exists():
        with sqlite3.connect(APP_DB) as con:
            sample.to_sql("sample", con, if_exists="replace", index=False)
            cv = pd.read_sql(
                "SELECT v.* FROM variants v JOIN (SELECT DISTINCT chrom, pos, ref, alt FROM sample) s "
                "ON v.chrom = s.chrom AND v.pos = s.pos AND v.ref = s.ref AND v.alt = s.alt", con)
            con.execute("DROP TABLE sample")
    else:
        cv = clinvar_variants()
    return sample.merge(cv, on=["chrom", "pos", "ref", "alt"], how="inner")


def build_app_db() -> None:
    """Write the app's SQLite file: every panel variant, and submissions only for variants that at
    least one lab calls pathogenic or likely pathogenic (nothing else can become a finding)."""
    v = clinvar_variants()[VAR_COLS].copy()
    s = clinvar_submissions()
    plp = set(s.loc[s["ClinicalSignificance"].isin(
        ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]), "VariationID"])
    s = s.loc[s["VariationID"].isin(plp), SUB_COLS]
    APP_DB.unlink(missing_ok=True)
    with sqlite3.connect(APP_DB) as con:
        v.to_sql("variants", con, index=False)
        s.to_sql("submissions", con, index=False)
        con.execute("CREATE INDEX v_pos ON variants (chrom, pos)")
        con.execute("CREATE INDEX s_vid ON submissions (VariationID)")
        con.execute("VACUUM")
    print(f"{APP_DB.name}: {len(v):,} variants, {len(s):,} submissions, {APP_DB.stat().st_size / 1e6:.0f} MB")


if __name__ == "__main__":
    build_app_db()


def panel() -> pd.DataFrame:
    return pd.read_csv(DATA / "panel_genes.csv")
