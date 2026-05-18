#!/usr/bin/env python3
"""
fetch_uniprot.py — fetch all Swiss-Prot entries with EC 2.7.7.*
Saves: ec277_swissprot.tsv (metadata + Pfam) and ec277_swissprot.fasta.
"""
import sys, time, gzip, io
from pathlib import Path
import requests

WORKDIR = Path.home() / "phylo277"
WORKDIR.mkdir(exist_ok=True)

QUERY = "(ec:2.7.7.*) AND (reviewed:true)"

FIELDS = ",".join([
    "accession",
    "id",
    "protein_name",
    "ec",
    "length",
    "organism_name",
    "organism_id",
    "lineage",
    "xref_pfam",
    "sequence",
])

STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"

def stream_tsv():
    print(f"[*] Streaming TSV (query: {QUERY})", flush=True)
    params = {
        "query": QUERY,
        "format": "tsv",
        "fields": FIELDS,
        "compressed": "true",
    }
    r = requests.get(STREAM_URL, params=params, stream=True, timeout=300)
    r.raise_for_status()
    raw = r.content
    print(f"[*] downloaded {len(raw):,} bytes (gzip)", flush=True)
    return gzip.decompress(raw).decode("utf-8")


def stream_fasta():
    print(f"[*] Streaming FASTA", flush=True)
    params = {
        "query": QUERY,
        "format": "fasta",
        "compressed": "true",
    }
    r = requests.get(STREAM_URL, params=params, stream=True, timeout=300)
    r.raise_for_status()
    raw = r.content
    print(f"[*] downloaded {len(raw):,} bytes (gzip)", flush=True)
    return gzip.decompress(raw).decode("utf-8")


def main():
    tsv_path = WORKDIR / "ec277_swissprot.tsv"
    fasta_path = WORKDIR / "ec277_swissprot.fasta"

    if not tsv_path.exists():
        tsv = stream_tsv()
        tsv_path.write_text(tsv, encoding="utf-8")
    else:
        tsv = tsv_path.read_text(encoding="utf-8")
        print(f"[*] reusing existing {tsv_path}", flush=True)

    if not fasta_path.exists():
        fasta = stream_fasta()
        fasta_path.write_text(fasta, encoding="utf-8")
    else:
        print(f"[*] reusing existing {fasta_path}", flush=True)

    lines = tsv.strip().split("\n")
    print(f"\n[+] {len(lines)-1:,} Swiss-Prot entries (excluding header)")
    print(f"[+] TSV   → {tsv_path}")
    print(f"[+] FASTA → {fasta_path}")

    # Quick stats: Pfam distribution
    hdr = lines[0].split("\t")
    pfam_idx = hdr.index("Pfam")
    ec_idx = hdr.index("EC number")
    pfam_count, no_pfam = {}, 0
    ec_count = {}
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) <= pfam_idx:
            no_pfam += 1
            continue
        pfams = [p for p in cells[pfam_idx].rstrip(";").split(";") if p]
        if not pfams:
            no_pfam += 1
        for p in pfams:
            pfam_count[p] = pfam_count.get(p, 0) + 1
        for ec in cells[ec_idx].rstrip(";").split("; "):
            ec = ec.strip()
            if ec:
                ec_count[ec] = ec_count.get(ec, 0) + 1

    print(f"\n[+] entries without Pfam: {no_pfam}")
    print(f"[+] distinct Pfam accessions: {len(pfam_count)}")
    print(f"\n[+] top 20 Pfam families by count:")
    for p, c in sorted(pfam_count.items(), key=lambda x: -x[1])[:20]:
        print(f"    {p}  {c:>4}")
    print(f"\n[+] EC sub-classes (top 20):")
    for e, c in sorted(ec_count.items(), key=lambda x: -x[1])[:20]:
        print(f"    {e:<14} {c:>4}")


if __name__ == "__main__":
    main()
