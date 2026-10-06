"""Hugging Face Space entry point: the Gradio UI, calling the backend functions in-process.

The FastAPI app in api/main.py exposes the same functions over HTTP for a separate deploy
(uvicorn api.main:app); in the Space the UI calls them directly, since the Space runner owns the port.
Data files too large for git (ClinVar SQLite, the two sample VCFs) are downloaded once from the
GitHub release on first start.

Run locally:  python app.py   then open http://localhost:7860
"""
import gzip
import shutil
from pathlib import Path

try:  # Hugging Face ZeroGPU hosting needs one GPU-decorated function; this app never calls it
    import spaces

    @spaces.GPU
    def _unused_gpu():
        return None
except ImportError:  # local runs
    pass

import requests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RELEASE = "https://github.com/Manas-repo/100x-Capstone-Genome-Intelligence/releases/download/data-v1/"
FILES = {  # release asset -> local path
    "clinvar_app.sqlite.gz": DATA / "clinvar_app.sqlite",
    "dev_sample.vcf.gz": DATA / "samples" / "dev_sample.vcf.gz",
    "heldout_sample.vcf.gz": DATA / "samples" / "heldout_sample.vcf.gz",
}


def fetch_data() -> None:
    for asset, dest in FILES.items():
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {asset}")
        with requests.get(RELEASE + asset, stream=True, timeout=(30, 1800)) as r:
            r.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as f:
                src = gzip.GzipFile(fileobj=r.raw) if asset == "clinvar_app.sqlite.gz" else r.raw
                shutil.copyfileobj(src, f, 1 << 20)
            tmp.rename(dest)


fetch_data()

from frontend.app import CSS, demo  # noqa: E402

if __name__ == "__main__":
    demo.launch(css=CSS)
