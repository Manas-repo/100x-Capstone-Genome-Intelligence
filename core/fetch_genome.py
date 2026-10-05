"""Pull two real South Asian genomes from the 1000 Genomes high-coverage (NYGC, GRCh38) callset,
restricted to the ACMG SF panel genes, and write one single-sample VCF per person.

- dev sample:      used while building
- held-out sample: not opened until the user test (only its file is written)

Samples are picked with a fixed random seed from the South Asian (SAS) superpopulation, so
neither file is hand-picked. The remote VCFs are bgzipped and tabix-indexed; this script reads
the .tbi index and fetches only the byte ranges covering panel genes (pure Python, no bcftools).
"""
import gzip
import io
import random
import json
import struct
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SAMPLES_DIR = DATA / "samples"

BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
VCF = BASE + "/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PED = BASE + "/20130606_g1k_3202_samples_ped_population.txt"
SEED = 20261007  # demo day; recorded so the pick is reproducible
FLANK = 1000  # bp around each gene, to keep splice sites near the ends


# ---------- minimal remote tabix ----------

def _get(url: str, **kw) -> requests.Response:
    """GET with retries; the EBI server has short outages."""
    for attempt in range(6):
        try:
            r = requests.get(url, timeout=(30, 600), **kw)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            if attempt == 5:
                raise
            print(f"    retry {attempt + 1} after error: {type(e).__name__}")
            time.sleep(10 * (attempt + 1))


def _bgzf_blocks(raw: bytes):
    """Yield decompressed BGZF blocks from a byte string (stops at a truncated block)."""
    pos = 0
    while pos + 18 <= len(raw):
        bsize = struct.unpack_from("<H", raw, pos + 16)[0] + 1
        if pos + bsize > len(raw):
            return
        yield zlib.decompress(raw[pos + 18 : pos + bsize - 8], -15)
        pos += bsize


def _read_tbi(url: str) -> dict:
    raw = _get(url + ".tbi").content
    data = gzip.decompress(raw)
    buf = io.BytesIO(data)
    magic, n_ref, *_ = struct.unpack("<4s7i", buf.read(32))
    assert magic == b"TBI\x01"
    names = buf.read(struct.unpack("<i", buf.read(4))[0]).split(b"\x00")[:-1]
    index = {}
    for name in names:
        bins = {}
        for _ in range(struct.unpack("<i", buf.read(4))[0]):
            b, n_chunk = struct.unpack("<Ii", buf.read(8))
            bins[b] = [struct.unpack("<QQ", buf.read(16)) for _ in range(n_chunk)]
        n_intv = struct.unpack("<i", buf.read(4))[0]
        linear = struct.unpack(f"<{n_intv}Q", buf.read(8 * n_intv))
        index[name.decode()] = (bins, linear)
    return index


def _reg2bins(beg: int, end: int) -> list[int]:
    end -= 1
    bins = [0]
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        bins.extend(range(offset + (beg >> shift), offset + (end >> shift) + 1))
    return bins


def _fetch_region(url: str, index: dict, chrom: str, beg: int, end: int) -> list[str]:
    bins, linear = index[chrom]
    min_off = linear[beg >> 14] if (beg >> 14) < len(linear) else 0
    chunks = [c for b in _reg2bins(beg, end) for c in bins.get(b, []) if c[1] > min_off]
    if not chunks:
        return []
    # Chunks can sit far apart in the file (large structural variants live in coarse bins),
    # so fetch each group of nearby chunks on its own instead of one min-to-max span.
    chunks.sort()
    groups = [list(chunks[0])]
    for b, e in chunks[1:]:
        if (b >> 16) <= (groups[-1][1] >> 16) + 0x10000:
            groups[-1][1] = max(groups[-1][1], e)
        else:
            groups.append([b, e])
    found = {}
    for start_v, end_v in groups:
        c_start, u_start = start_v >> 16, start_v & 0xFFFF
        c_end = (end_v >> 16) + 0x10000  # include the whole last block
        raw = _get(url, headers={"Range": f"bytes={c_start}-{c_end}"}).content
        text = b"".join(_bgzf_blocks(raw))[u_start:].decode()
        for line in text.split("\n")[:-1]:  # last line may be cut off
            f = line.split("\t", 5)
            if len(f) < 6 or f[0] != chrom:
                continue
            pos = int(f[1])
            if pos > end:
                break
            if pos >= beg:
                found[(pos, f[3], f[4])] = line
    return [found[k] for k in sorted(found)]


def _header(url: str) -> list[str]:
    raw = _get(url, headers={"Range": "bytes=0-4000000"}).content
    text = b"".join(_bgzf_blocks(raw)).decode(errors="replace")
    lines = []
    for line in text.split("\n"):
        lines.append(line)
        if line.startswith("#CHROM"):
            return lines
    raise RuntimeError("no #CHROM line in the first 4 MB")


# ---------- sample choice and writing ----------

def pick_samples() -> tuple[str, str, pd.DataFrame]:
    ped = pd.read_csv(PED, sep=r"\s+")
    sas = ped[ped["Superpopulation"] == "SAS"].sort_values("SampleID")
    dev, held = random.Random(SEED).sample(list(sas["SampleID"]), 2)
    return dev, held, ped.set_index("SampleID")


def regions(panel: pd.DataFrame) -> pd.DataFrame:
    """What to download.
    - ACMG SF genes: the whole gene plus flanks.
    - Carrier genes: only windows around positions where at least one ClinVar lab calls a variant
      pathogenic or likely pathogenic (anything else could never become a finding), merged when
      closer than 2 kb. Keeps the download small.
    """
    acmg = panel[panel["panel"] == "acmg_sf"]
    rows = [("acmg", g.chrom, int(g.start) - FLANK, int(g.end) + FLANK) for g in acmg.itertuples()]

    car = panel[panel["panel"] == "carrier"]
    v = pd.read_parquet(DATA / "clinvar_variants.parquet", columns=["VariationID", "Chromosome", "PositionVCF"])
    s = pd.read_parquet(DATA / "clinvar_submissions.parquet", columns=["VariationID", "ClinicalSignificance"])
    plp = set(s.loc[s["ClinicalSignificance"].isin(
        ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]), "VariationID"].astype(int))
    v = v[v["VariationID"].isin(plp)]
    v_chrom = "chr" + v["Chromosome"].astype(str)
    for ch, genes in car.groupby("chrom"):
        pos = v.loc[v_chrom == ch, "PositionVCF"].to_numpy()
        inside = np.zeros(len(pos), dtype=bool)
        for g in genes.itertuples():
            inside |= (pos >= g.start) & (pos <= g.end)
        pos = np.sort(pos[inside])
        if not len(pos):
            continue
        start = end = int(pos[0])
        for x in pos[1:]:
            if x - end > 2000:
                rows.append(("carrier", ch, start - 50, end + 50))
                start = int(x)
            end = int(x)
        rows.append(("carrier", ch, start - 50, end + 50))
    return pd.DataFrame(rows, columns=["kind", "chrom", "start", "end"])


def main() -> None:
    SAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    cache_dir = DATA / "raw" / "1kg_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    dev, held, ped = pick_samples()
    picks = {"dev": dev, "heldout": held}
    print(f"seed {SEED}: dev={dev} ({ped.loc[dev, 'Population']}), held-out sample chosen (population hidden)")

    panel = pd.read_csv(DATA / "panel_genes.csv")
    head_path = cache_dir / "header.json"
    if not head_path.exists():
        head_path.write_text(json.dumps(_header(VCF.format(chrom="chr1"))))
    head = json.loads(head_path.read_text())
    cols = head[-1].split("\t")
    col_of = {k: cols.index(s) for k, s in picks.items()}

    regs = regions(panel)
    print(f"regions: {len(regs)} ({(regs['kind'] == 'carrier').sum()} carrier windows)")

    # Per (kind, chromosome): keep only the two picked samples' non-reference sites; cache so reruns resume.
    for (kind, chrom), rs in regs.groupby(["kind", "chrom"], sort=False):
        out = cache_dir / (f"{chrom}.json" if kind == "acmg" else f"carrier_{chrom}.json")
        if out.exists():
            continue
        url = VCF.format(chrom=chrom)
        if chrom == "chrX":  # chrX was re-released as v2 on the server
            url = url.replace("phased_panel.vcf.gz", "phased_panel.v2.vcf.gz")
        index = _read_tbi(url)

        def one(r):
            kept = {k: [] for k in picks}
            for line in _fetch_region(url, index, chrom, int(r.start), int(r.end)):
                f = line.split("\t")
                for k, col in col_of.items():
                    if f[col].split(":")[0].replace("|", "/") not in ("0/0", "./.", "0", "."):
                        kept[k].append("\t".join(f[:9] + [f[col]]))
            return kept

        with ThreadPoolExecutor(max_workers=8) as pool:
            results = list(pool.map(one, rs.itertuples()))
        merged = {
            k: sorted({l for r in results for l in r[k]}, key=lambda l: int(l.split("\t")[1])) for k in picks
        }
        out.write_text(json.dumps(merged))
        print(f"  {kind} {chrom}: done ({len(rs)} regions)")

    counts = {}
    order = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    for k, s in picks.items():
        with gzip.open(SAMPLES_DIR / f"{k}_sample.vcf.gz", "wt") as w:
            w.writelines(l + "\n" for l in head[:-1] if l.startswith("##"))
            w.write(f"##source_sample=1000 Genomes high coverage {s}; ACMG SF v3.2 genes and "
                    "ClinVar-pathogenic positions in recessive carrier genes only\n")
            w.write("\t".join(cols[:9] + [s]) + "\n")
            counts[k] = 0
            for chrom in order:
                lines = set()
                for name in (f"{chrom}.json", f"carrier_{chrom}.json"):
                    p = cache_dir / name
                    if p.exists():
                        lines |= set(json.loads(p.read_text())[k])
                lines = sorted(lines, key=lambda l: int(l.split("\t")[1]))
                w.writelines(l + "\n" for l in lines)
                counts[k] += len(lines)
    (SAMPLES_DIR / "README.md").write_text(
        f"Samples picked from 1000 Genomes SAS with random seed {SEED}.\n"
        f"- dev_sample.vcf.gz: {dev} ({ped.loc[dev, 'Population']}), {counts['dev']} variant sites\n"
        f"- heldout_sample.vcf.gz: not opened until the user test ({counts['heldout']} variant sites)\n"
        "Source: 1000G high-coverage phased panel (singletons are not in this callset).\n"
    )
    print("variant sites written:", counts)


if __name__ == "__main__":
    main()
