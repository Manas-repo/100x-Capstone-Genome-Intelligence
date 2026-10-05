"""Hugging Face Space entry point: the FastAPI backend and the Gradio UI in one process.

The API stays at /health, /samples, /analyse, /feedback; the UI is served at /.
Data files too large for git (ClinVar SQLite, the two sample VCFs) are downloaded once from the
GitHub release on first start.

Run locally:  python app.py   then open http://localhost:7860
"""
import gzip
import os
import shutil
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PORT = int(os.getenv("PORT", "7860"))
RELEASE = "https://github.com/Manas-repo/100x-Capstone-Genome-Intelligence/releases/download/data-v1/"
FILES = {  # release asset -> local path
    "clinvar_app.sqlite.gz": DATA / "clinvar_app.sqlite",
    "dev_sample.vcf.gz": DATA / "samples" / "dev_sample.vcf.gz",
    "heldout_sample.vcf.gz": DATA / "samples" / "heldout_sample.vcf.gz",
    "gnomad_cache.json": DATA / "gnomad_cache.json",
    "llm_cache.json": DATA / "llm_cache.json",
}


def fetch_data() -> None:
    for asset, dest in FILES.items():
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {asset}")
        with requests.get(RELEASE + asset, stream=True, timeout=(30, 1800)) as r:
            if r.status_code == 404 and asset.endswith(".json"):
                continue  # caches are optional
            r.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as f:
                src = gzip.GzipFile(fileobj=r.raw) if asset == "clinvar_app.sqlite.gz" else r.raw
                shutil.copyfileobj(src, f, 1 << 20)
            tmp.rename(dest)


fetch_data()
os.environ.setdefault("BACKEND_URL", f"http://127.0.0.1:{PORT}")
sys.path.insert(0, str(ROOT / "frontend"))

import gradio as gr  # noqa: E402
import uvicorn  # noqa: E402

from api.main import app as api  # noqa: E402
from frontend.app import demo  # noqa: E402

app = gr.mount_gradio_app(api, demo, path="/")

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=PORT)
