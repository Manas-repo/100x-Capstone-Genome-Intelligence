---
title: Genome Intelligence
emoji: 🧬
colorFrom: green
colorTo: gray
sdk: gradio
app_file: app.py
pinned: false
---

# 100x Capstone: Genome Intelligence

Takes a real variant file (VCF) and returns a few challengeable findings for germline monogenic disease genes (ACMG secondary findings list). Each finding shows supporting evidence, disputing evidence, evidence age and population match. Findings whose evidence does not hold up are declined, not reported.

Sense-making, not diagnosis. Not a substitute for a clinician or genetic counselor.

Built for the 100X Engineers C7 Capstone (Genome Intelligence brief).

## Stack

- Backend: FastAPI on Render
- Frontend: Gradio on Hugging Face Spaces
- LLM: Groq
- Storage (user feedback on findings): Supabase (Postgres)

## Status

Work in progress.
