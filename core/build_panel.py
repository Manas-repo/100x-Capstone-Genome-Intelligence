"""Build data/panel_genes.csv: ACMG SF v3.2 genes with GRCh38 coordinates and disease mechanism.

Sources (downloaded into data/raw/ by hand, see README):
- NCBI ClinVar ACMG page (gene list and conditions): https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
- ClinGen gene curation list, GRCh38 (coordinates, haploinsufficiency score):
  https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
"""
import html
import re
from pathlib import Path

import pandas as pd
import requests

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
OUT = Path(__file__).resolve().parent.parent / "data" / "panel_genes.csv"

# Inheritance that changes what counts as a finding: one hit in a recessive gene is carrier status.
RECESSIVE = {"ATP7B", "BTD", "CASQ2", "GAA", "HFE", "MUTYH", "RPE65", "TRDN"}
X_LINKED = {"GLA", "OTC"}

# Genes where disease is NOT (only) caused by loss of function, so "null variant, therefore LP"
# (PVS1) does not hold. This is the error pattern the curator named in the 2026-10-03 interview.
NOT_LOF = {
    "RET": "gain of function",
    "TTR": "gain of function (amyloid)",
    "RYR1": "gain of function (malignant hyperthermia)",
    "CACNA1S": "gain of function (malignant hyperthermia)",
    "MYH7": "dominant negative missense",
    "MYL2": "missense",
    "MYL3": "missense",
    "ACTC1": "missense",
    "TNNT2": "missense",
    "TNNI3": "missense",
    "TNNC1": "missense",
    "TPM1": "missense",
    "PRKAG2": "gain of function",
    "ACTA2": "missense (dominant negative)",
    "MYH11": "missense",
    "CALM1": "missense",
    "CALM2": "missense",
    "CALM3": "missense",
    "RYR2": "gain of function",
    "PCSK9": "gain of function",
    "HFE": "specific variant only (p.C282Y homozygous)",
    "TMEM43": "single missense founder variant (p.S358L)",
    "RBM20": "missense hotspot (RS domain)",
    "TGFBR2": "missense (kinase domain)",
    "TGFBR1": "missense (kinase domain)",
}

# Genes where both mechanisms cause disease, so a null variant may or may not fit the condition.
MIXED = {
    "SCN5A": "loss of function (Brugada) and gain of function (long QT 3)",
    "LMNA": "loss of function and missense",
    "SMAD3": "loss of function and missense",
}

# ACMG reporting restrictions written on the NCBI page.
RESTRICTIONS = {
    "TTN": "truncating variants only",
    "HFE": "p.C282Y homozygotes only",
}


def acmg_genes() -> pd.DataFrame:
    s = (RAW / "acmg_page.html").read_text(encoding="utf8", errors="ignore")
    t = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", s)))
    t = t.replace("MedGenC", "MedGen C")  # typo on the page that hides RPE65
    body = t[t.find("Gene via GTR Variations that may be pathogenic") :]
    body = body[len("Gene via GTR Variations that may be pathogenic ") :]
    rows = re.findall(r"(.+?) MedGen [A-Z0-9, ]+? ([A-Z0-9]+) \(MIM \d+\) ClinVar ", body)
    conditions: dict[str, list[str]] = {}
    for cond, gene in rows:
        cond = re.sub(r"\s*\( MIM [\d ,MI]+\)", "", cond).strip()
        conditions.setdefault(gene, [])
        if cond not in conditions[gene]:
            conditions[gene].append(cond)
    return pd.DataFrame({"gene": list(conditions), "conditions": ["; ".join(c) for c in conditions.values()]})


def clingen() -> pd.DataFrame:
    df = pd.read_csv(RAW / "clingen_curation_grch38.tsv", sep="\t", skiprows=5)
    df = df.rename(columns={"#Gene Symbol": "gene", "Genomic Location": "loc", "Haploinsufficiency Score": "hi_score"})
    loc = df["loc"].str.extract(r"(chr[\dXY]+):(\d+)-(\d+)")
    df["chrom"], df["start"], df["end"] = loc[0], loc[1].astype("Int64"), loc[2].astype("Int64")
    return df[["gene", "chrom", "start", "end", "hi_score"]]


def fill_missing_coords(panel: pd.DataFrame) -> None:
    """ClinGen lacks a few genes; take GRCh38 coordinates from the Ensembl REST API."""
    for i, row in panel[panel["chrom"].isna()].iterrows():
        r = requests.get(
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{row['gene']}",
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        r.raise_for_status()
        j = r.json()
        panel.loc[i, ["chrom", "start", "end"]] = [f"chr{j['seq_region_name']}", j["start"], j["end"]]


def main() -> None:
    panel = acmg_genes().merge(clingen(), on="gene", how="left")
    panel["inheritance"] = panel["gene"].map(
        lambda g: "AR" if g in RECESSIVE else ("XL" if g in X_LINKED else "AD")
    )
    panel["lof_mechanism"] = panel["gene"].map(
        lambda g: "mixed" if g in MIXED else ("no" if g in NOT_LOF else "yes")
    )
    panel["mechanism_note"] = panel["gene"].map({**NOT_LOF, **MIXED}).fillna("loss of function")
    fill_missing_coords(panel)
    panel["acmg_restriction"] = panel["gene"].map(RESTRICTIONS).fillna("")
    missing = panel[panel["chrom"].isna()]["gene"].tolist()
    if missing:
        print("WARNING no coordinates for:", missing)
    panel.sort_values("gene").to_csv(OUT, index=False)
    print(f"{len(panel)} genes written to {OUT}")


if __name__ == "__main__":
    main()
