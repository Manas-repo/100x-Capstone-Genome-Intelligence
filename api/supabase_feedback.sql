-- T12b: run once in the Supabase SQL editor.
create table if not exists feedback (
  id uuid primary key,
  created_at timestamptz not null default now(),
  session_id text not null,
  sample text not null,
  variation_id bigint not null,
  gene text,
  source text not null check (source in ('system', 'baseline')),
  shown_verdict text,
  answer text not null check (answer in ('agree', 'disagree')),
  reason text
);

-- Only the backend (service-role key) writes; nobody reads through the public API.
alter table feedback enable row level security;
