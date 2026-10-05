"""Add recessive carrier-screening genes to data/panel_genes.csv.

Selection (data-driven, reproducible):
1. ClinGen gene-disease validity: autosomal recessive, classification Definitive or Strong.
2. Keep genes where at least one variant seen in 1000 Genomes carries a ClinVar pathogenic or
   likely pathogenic annotation (from the 1000G "annotated.clinical.txt" files). These are the
   genes where a healthy person can plausibly turn out to be a carrier.

Sources:
- https://search.clinicalgenome.org/kb/gene-validity/download
- http://ftp.1000genomes.ebi.ac.uk/.../20201028_3202_raw_GT_with_annot/*chrN*.annotated.clinical.txt
"""
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PANEL = ROOT / "data" / "panel_genes.csv"
CLIN_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
    "20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}"
    ".recalibrated_variants.annotated.clinical.txt"
)


def ar_genes() -> pd.DataFrame:
    v = pd.read_csv(RAW / "clingen_validity.csv", skiprows=[0, 1, 2, 3, 5])
    v = v[(v["MOI"] == "AR") & v["CLASSIFICATION"].isin(["Definitive", "Strong"])]
    return v.groupby("GENE SYMBOL")["DISEASE LABEL"].apply(lambda d: "; ".join(sorted(set(d)))).reset_index()


def genes_with_1kg_pathogenic() -> set[str]:
    out_dir = RAW / "1kg_clinical"
    out_dir.mkdir(exist_ok=True)
    genes = set()
    for c in [str(i) for i in range(1, 23)] + ["X"]:
        p = out_dir / f"chr{c}.txt"
        if not p.exists():
            r = requests.get(CLIN_URL.format(chrom=f"chr{c}"), timeout=(30, 300))
            r.raise_for_status()
            p.write_bytes(r.content)
        df = pd.read_csv(p, sep="\t", dtype=str, usecols=["GENE", "CLNSIG", "FILTER"])
        # Old ClinVar numeric codes: 5 = pathogenic, 4 = likely pathogenic.
        sig = df["CLNSIG"].fillna("").str.split(r"[|,]")
        hit = sig.map(lambda xs: any(x.strip() in ("4", "5") for x in xs)) & (df["FILTER"] == "PASS")
        genes |= set(df.loc[hit, "GENE"].dropna())
    return genes


def coords(genes: list[str]) -> pd.DataFrame:
    """GRCh38 coordinates via the Ensembl REST batch lookup."""
    rows = []
    for i in range(0, len(genes), 500):
        r = requests.post(
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            json={"symbols": genes[i : i + 500]},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=120,
        )
        r.raise_for_status()
        for sym, j in r.json().items():
            if j and j.get("seq_region_name") in [str(x) for x in range(1, 23)] + ["X"]:
                rows.append({"gene": sym, "chrom": f"chr{j['seq_region_name']}", "start": j["start"], "end": j["end"]})
    return pd.DataFrame(rows)


def main() -> None:
    panel = pd.read_csv(PANEL)
    panel = panel[panel.get("panel", pd.Series("acmg_sf", index=panel.index)).fillna("acmg_sf") == "acmg_sf"].copy()
    panel["panel"] = "acmg_sf"

    ar = ar_genes()
    seen = genes_with_1kg_pathogenic()
    pick = ar[ar["GENE SYMBOL"].isin(seen) & ~ar["GENE SYMBOL"].isin(panel["gene"])]
    print(f"ClinGen AR definitive/strong genes: {len(ar)}; with a ClinVar P/LP allele in 1000G: {len(pick)}")

    c = coords(sorted(pick["GENE SYMBOL"]))
    carrier = pick.rename(columns={"GENE SYMBOL": "gene", "DISEASE LABEL": "conditions"}).merge(c, on="gene")
    carrier["hi_score"] = 30.0
    carrier["inheritance"] = carrier["chrom"].map(lambda ch: "XL" if ch == "chrX" else "AR")
    carrier["lof_mechanism"] = "yes"  # recessive disease is overwhelmingly loss of function
    carrier["mechanism_note"] = "loss of function (recessive)"
    carrier["acmg_restriction"] = ""
    carrier["panel"] = "carrier"

    out = pd.concat([panel, carrier[panel.columns]], ignore_index=True).sort_values(["panel", "gene"])
    out.to_csv(PANEL, index=False)
    span = (carrier["end"] - carrier["start"]).sum() / 1e6
    print(f"carrier genes added: {len(carrier)} (total span {span:.1f} Mb); panel now {len(out)} genes")


if __name__ == "__main__":
    main()
