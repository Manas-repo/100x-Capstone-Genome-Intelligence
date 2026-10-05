"""One-off check that the Groq key works and returns valid JSON."""
import json
import os

from dotenv import load_dotenv
from groq import Groq

load_dotenv()

key = os.getenv("GROQ_API_KEY")
model = os.getenv("GROQ_MODEL", "openai/gpt-oss-120b")
if not key:
    raise SystemExit("GROQ_API_KEY is empty. Paste it into the .env file first.")

client = Groq(api_key=key)
resp = client.chat.completions.create(
    model=model,
    response_format={"type": "json_object"},
    messages=[
        {"role": "system", "content": "Reply in JSON only."},
        {
            "role": "user",
            "content": 'A lab classified a variant as Likely Pathogenic citing "found in trans with a '
            'pathogenic variant in an affected child". Does that reasoning apply to a healthy, '
            'unrelated adult? Reply as {"applies": "yes|no|cant_tell", "reason": "..."}',
        },
    ],
)
out = json.loads(resp.choices[0].message.content)
print(f"model: {model}")
print(json.dumps(out, indent=2))
