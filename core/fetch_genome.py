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
import struct
import zlib
from pathlib import Path

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
    raw = requests.get(url + ".tbi", timeout=120).content
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
        raw = requests.get(url, headers={"Range": f"bytes={c_start}-{c_end}"}, timeout=600).content
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
    raw = requests.get(url, headers={"Range": "bytes=0-4000000"}, timeout=300).content
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


def main() -> None:
    SAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    dev, held, ped = pick_samples()
    picks = {"dev": dev, "heldout": held}
    print(f"seed {SEED}: dev={dev} ({ped.loc[dev, 'Population']}), held-out sample chosen (population hidden)")

    panel = pd.read_csv(DATA / "panel_genes.csv")
    writers, header_done, counts = {}, False, {k: 0 for k in picks}
    for chrom, genes in panel.groupby("chrom", sort=False):
        url = VCF.format(chrom=chrom)
        index = _read_tbi(url)
        if not header_done:
            head = _header(url)
            cols = head[-1].split("\t")
            col_of = {k: cols.index(s) for k, s in picks.items()}
            for k, s in picks.items():
                w = gzip.open(SAMPLES_DIR / f"{k}_sample.vcf.gz", "wt")
                w.writelines(l + "\n" for l in head[:-1] if l.startswith("##"))
                w.write(f"##source_sample=1000 Genomes high coverage {s}, ACMG SF v3.2 panel regions only\n")
                w.write("\t".join(cols[:9] + [s]) + "\n")
                writers[k] = w
            header_done = True
        for g in genes.itertuples():
            for line in _fetch_region(url, index, chrom, int(g.start) - FLANK, int(g.end) + FLANK):
                f = line.split("\t")
                for k, col in col_of.items():
                    gt = f[col].split(":")[0]
                    if gt.replace("|", "/") in ("0/0", "./.", "0", "."):
                        continue
                    writers[k].write("\t".join(f[:9] + [f[col]]) + "\n")
                    counts[k] += 1
        print(f"  {chrom}: done ({len(genes)} genes)")
    for w in writers.values():
        w.close()
    (SAMPLES_DIR / "README.md").write_text(
        f"Samples picked from 1000 Genomes SAS with random seed {SEED}.\n"
        f"- dev_sample.vcf.gz: {dev} ({ped.loc[dev, 'Population']}), {counts['dev']} variant sites\n"
        f"- heldout_sample.vcf.gz: not opened until the user test ({counts['heldout']} variant sites)\n"
    )
    print("variant sites written:", counts)


if __name__ == "__main__":
    main()
