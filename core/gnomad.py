"""gnomAD v4 allele frequencies (overall vs South Asian), cached on disk so the app does not
depend on the gnomAD API being up."""
import json
import time
from pathlib import Path

import requests

CACHE = Path(__file__).resolve().parent.parent / "data" / "gnomad_cache.json"
API = "https://gnomad.broadinstitute.org/api"
QUERY = """
query($id: String!) {
  variant(variantId: $id, dataset: gnomad_r4) {
    exome { ac an populations { id ac an } }
    genome { ac an populations { id ac an } }
  }
}"""


def _load() -> dict:
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def lookup(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """Return {'af_all', 'af_sas', 'an_sas', 'found'} combining gnomAD v4 exomes and genomes."""
    vid = f"{chrom}-{pos}-{ref}-{alt}"
    cache = _load()
    if vid in cache:
        return cache[vid]
    for attempt in range(3):
        r = requests.post(API, json={"query": QUERY, "variables": {"id": vid}}, timeout=60)
        if r.status_code == 429:
            time.sleep(5 * (attempt + 1))
            continue
        break
    body = r.json()
    v = (body.get("data") or {}).get("variant")
    if not v:
        out = {"found": False, "af_all": 0.0, "af_sas": 0.0, "an_all": 0, "an_sas": 0}
    else:
        ac = an = ac_sas = an_sas = 0
        for part in (v.get("exome"), v.get("genome")):
            if not part:
                continue
            ac += part["ac"]
            an += part["an"]
            for p in part["populations"]:
                if p["id"] == "sas":
                    ac_sas += p["ac"]
                    an_sas += p["an"]
        out = {
            "found": True,
            "af_all": ac / an if an else 0.0,
            "af_sas": ac_sas / an_sas if an_sas else 0.0,
            "an_all": an,
            "an_sas": an_sas,
        }
    cache[vid] = out
    CACHE.write_text(json.dumps(cache, indent=1))
    return out


if __name__ == "__main__":
    # BRCA2 c.68-7T>A is a known example; just checks the API answers.
    print(lookup("13", 32316508, "T", "A"))
