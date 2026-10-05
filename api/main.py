"""FastAPI backend (T12a): wraps the pipeline and stores the user's agree/disagree answers.

Run locally:  uvicorn api.main:app --reload   (from the repo root)

Endpoints
- GET  /health            liveness plus whether the model and Supabase are configured
- GET  /samples           the bundled 1000 Genomes samples
- POST /analyse/{sample}  baseline list and system findings for a bundled sample
- POST /analyse           same, for an uploaded VCF
- POST /feedback          one agree/disagree answer on one finding
"""
import json
import math
import os
import sys
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

from dotenv import load_dotenv
from fastapi import FastAPI, File, HTTPException, UploadFile
from pydantic import BaseModel

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "core"))
load_dotenv(ROOT / ".env")

from adjudicate import analyse  # noqa: E402
from baseline import run_baseline  # noqa: E402

SAMPLES = {"dev": "dev_sample.vcf.gz", "heldout": "heldout_sample.vcf.gz"}
LOCAL_FEEDBACK = ROOT / "results" / "feedback_local.jsonl"
_cache: dict[str, dict] = {}

app = FastAPI(title="Genome Intelligence: challengeable findings")


def _clean(x):
    """Make pipeline output JSON-safe (numpy scalars, NaN)."""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if hasattr(x, "item"):  # numpy scalar
        x = x.item()
    if isinstance(x, float) and math.isnan(x):
        return None
    return x


def _run(path: Path) -> dict:
    base = run_baseline(path)
    return _clean({"baseline": base.to_dict(orient="records"), "system": analyse(path)})


@app.get("/health")
def health():
    return {
        "ok": True,
        "model": bool(os.getenv("GROQ_API_KEY")),
        "supabase": bool(os.getenv("SUPABASE_URL") and os.getenv("SUPABASE_SERVICE_ROLE_KEY")),
    }


@app.get("/samples")
def samples():
    return [{"id": k, "label": "Development sample" if k == "dev" else "Held-out sample"} for k in SAMPLES]


def analyse_named(sample: str) -> dict:
    """Plain function behind POST /analyse/{sample}; the Gradio UI calls it directly when co-hosted."""
    if sample not in _cache:
        _cache[sample] = _run(ROOT / "data" / "samples" / SAMPLES[sample])
    return {"sample": sample, **_cache[sample]}


def analyse_file(path: Path, name: str) -> dict:
    return {"sample": name, **_run(path)}


@app.post("/analyse/{sample}")
def analyse_sample(sample: str):
    if sample not in SAMPLES:
        raise HTTPException(404, f"unknown sample '{sample}'")
    return analyse_named(sample)


@app.post("/analyse")
async def analyse_upload(file: UploadFile = File(...)):
    name = file.filename or "upload.vcf"
    if not (name.endswith(".vcf") or name.endswith(".vcf.gz")):
        raise HTTPException(400, "Upload a .vcf or .vcf.gz file (GRCh38).")
    suffix = ".vcf.gz" if name.endswith(".gz") else ".vcf"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp.write(await file.read())
        path = Path(tmp.name)
    try:
        return analyse_file(path, name)
    except Exception as e:  # a malformed file should say so, not crash the server
        raise HTTPException(422, f"Could not read this VCF: {type(e).__name__}: {e}")
    finally:
        path.unlink(missing_ok=True)


class Feedback(BaseModel):
    session_id: str
    sample: str
    variation_id: int
    gene: str | None = None
    source: str  # "system" or "baseline"
    shown_verdict: str | None = None  # Defensible / Contested / Declined, or the baseline label
    answer: str  # "agree" or "disagree"
    reason: str | None = None


@app.post("/feedback")
def feedback(fb: Feedback):
    return save_feedback(fb.model_dump())


def save_feedback(answer: dict) -> dict:
    row = {**Feedback(**answer).model_dump(), "id": str(uuid.uuid4()),
           "created_at": datetime.now(timezone.utc).isoformat()}
    url, key = os.getenv("SUPABASE_URL"), os.getenv("SUPABASE_SERVICE_ROLE_KEY")
    if url and key:
        from supabase import create_client

        create_client(url, key).table("feedback").insert(row).execute()
        return {"saved": "supabase", "id": row["id"]}
    LOCAL_FEEDBACK.parent.mkdir(exist_ok=True)
    with LOCAL_FEEDBACK.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row) + "\n")
    return {"saved": "local", "id": row["id"]}
