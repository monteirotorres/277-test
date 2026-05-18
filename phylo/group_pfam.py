#!/usr/bin/env python3
"""
group_pfam.py — assign each Swiss-Prot EC 2.7.7.* entry to a primary Pfam family.

Strategy: each protein is assigned to the Pfam (from its annotations) that has the
HIGHEST total count across the dataset. This is a faithful "split-by-Pfam" mapping;
proteins without any Pfam go to a separate bucket.

Also fetches Pfam human-readable names via InterPro API.

Outputs in ~/phylo277/:
  protein_to_family.tsv
  family_summary.tsv
  families/<PFxxxxx>.fasta  (one FASTA per family with ≥ MIN_SIZE members)
"""
import json, time
from pathlib import Path
from collections import defaultdict, Counter
import requests

WORKDIR = Path.home() / "phylo277"
FAM_DIR = WORKDIR / "families"
FAM_DIR.mkdir(exist_ok=True)

MIN_SIZE = 8   # families with fewer members go to "small_other"

# ── Read TSV ────────────────────────────────────────────────
tsv_path   = WORKDIR / "ec277_swissprot.tsv"
fasta_path = WORKDIR / "ec277_swissprot.fasta"
rows = []
with tsv_path.open(encoding="utf-8") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for ln in fh:
        cells = ln.rstrip("\n").split("\t")
        if len(cells) < len(header):
            cells += [""] * (len(header) - len(cells))
        rows.append(dict(zip(header, cells)))
print(f"[*] {len(rows):,} rows loaded")

# Parse Pfam column
def parse_pfams(s):
    return [p for p in s.rstrip(";").split(";") if p]

for r in rows:
    r["pfams"] = parse_pfams(r.get("Pfam", ""))

# Global Pfam count
pfam_count = Counter()
for r in rows:
    for p in r["pfams"]:
        pfam_count[p] += 1

print(f"[*] {len(pfam_count)} distinct Pfams, top 10:")
for p, c in pfam_count.most_common(10):
    print(f"    {p}  {c:>5}")

# Assign primary Pfam = most globally-abundant Pfam in this protein's set
for r in rows:
    if not r["pfams"]:
        r["primary"] = "NO_PFAM"
    else:
        r["primary"] = max(r["pfams"], key=lambda p: pfam_count[p])

family_count = Counter(r["primary"] for r in rows)
print(f"\n[*] {len(family_count)} primary families.")
print("[*] families above MIN_SIZE:")
big = [(p, c) for p, c in family_count.most_common() if c >= MIN_SIZE and p != "NO_PFAM"]
for p, c in big[:30]:
    print(f"    {p}  {c:>5}")
print(f"\n[*] total in big families: {sum(c for _, c in big):,}/{len(rows):,}")
print(f"[*] in small_other: {sum(c for p, c in family_count.items() if p != 'NO_PFAM' and c < MIN_SIZE)}")
print(f"[*] NO_PFAM: {family_count.get('NO_PFAM', 0)}")

big_pfams = {p for p, _ in big}

# ── Pfam name lookup (InterPro API) ─────────────────────────
name_cache_path = WORKDIR / "pfam_names.json"
if name_cache_path.exists():
    pfam_names = json.loads(name_cache_path.read_text())
else:
    pfam_names = {}

unique_pfams = list(big_pfams - set(pfam_names))
print(f"\n[*] fetching {len(unique_pfams)} new Pfam names from InterPro...")
for i, pf in enumerate(unique_pfams):
    try:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}/"
        r = requests.get(url, timeout=20, headers={"Accept": "application/json"})
        if r.status_code == 200:
            d = r.json()["metadata"]
            pfam_names[pf] = {"name": d.get("name", {}).get("name") if isinstance(d.get("name"), dict) else d.get("name"), "type": d.get("type")}
        else:
            pfam_names[pf] = {"name": "?", "type": "?"}
    except Exception as e:
        pfam_names[pf] = {"name": f"err:{e.__class__.__name__}", "type": "?"}
    if (i + 1) % 10 == 0:
        print(f"    ... {i+1}/{len(unique_pfams)}", flush=True)
    time.sleep(0.1)

name_cache_path.write_text(json.dumps(pfam_names, indent=2))
print(f"[*] cached {len(pfam_names)} Pfam names → {name_cache_path}")

# ── Write protein_to_family.tsv ─────────────────────────────
p2f_path = WORKDIR / "protein_to_family.tsv"
with p2f_path.open("w", encoding="utf-8") as fh:
    fh.write("accession\tid\tec\tlength\torganism\tprimary_pfam\tprimary_pfam_name\tall_pfams\n")
    for r in rows:
        primary = r["primary"]
        name = pfam_names.get(primary, {}).get("name", "?") if primary != "NO_PFAM" else "no Pfam annotation"
        fh.write("\t".join([
            r["Entry"], r["Entry Name"], r.get("EC number", ""),
            r.get("Length", ""), r.get("Organism", ""),
            primary, str(name), ";".join(r["pfams"]),
        ]) + "\n")
print(f"[+] {p2f_path}")

# ── Family summary ──────────────────────────────────────────
fam_sum_path = WORKDIR / "family_summary.tsv"
with fam_sum_path.open("w", encoding="utf-8") as fh:
    fh.write("pfam_acc\tpfam_name\tcount\tsample_member\tsample_ec\n")
    for p, c in family_count.most_common():
        if c < MIN_SIZE or p == "NO_PFAM":
            continue
        sample = next(r for r in rows if r["primary"] == p)
        fh.write(f"{p}\t{pfam_names.get(p, {}).get('name', '?')}\t{c}\t{sample['Entry Name']}\t{sample.get('EC number', '')}\n")
print(f"[+] {fam_sum_path}")

# ── Write per-family FASTAs ─────────────────────────────────
print(f"\n[*] writing per-family FASTA files (≥ {MIN_SIZE} members) ...")

# Index FASTA
seqs = {}
acc = None
buf = []
with fasta_path.open(encoding="utf-8") as fh:
    for ln in fh:
        if ln.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc = ln[1:].split("|")[1] if "|" in ln else ln[1:].split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if acc:
        seqs[acc] = "".join(buf)
print(f"[*] {len(seqs):,} sequences indexed from FASTA")

# Group by primary
fam_members = defaultdict(list)
for r in rows:
    if r["primary"] in big_pfams:
        fam_members[r["primary"]].append(r["Entry"])

for pf, accs in fam_members.items():
    name = pfam_names.get(pf, {}).get("name", "?")
    fp = FAM_DIR / f"{pf}.fasta"
    with fp.open("w", encoding="utf-8") as fh:
        for a in accs:
            if a in seqs:
                fh.write(f">{a}\n{seqs[a]}\n")
    print(f"    {pf} ({name}): {len(accs)} → {fp.name}")

print(f"\n[done] {len(big_pfams)} families ready in {FAM_DIR}")
