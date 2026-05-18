#!/usr/bin/env python3
"""Unified summary of all built trees (IQ-TREE and FastTree)."""
import csv, re
from pathlib import Path
from statistics import median

WORK = Path.home() / "phylo277"
TREES = WORK / "trees"

# Load family summary for names + sizes
fam_meta = {}
with open(WORK / "family_summary.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        fam_meta[r["pfam_acc"]] = {
            "name": r["pfam_name"],
            "count": int(r["count"]),
            "ec": r.get("sample_ec", ""),
        }

def count_taxa_in_newick(s):
    # Count commas + 1, more or less — but better just count leaf labels
    # Leaf labels appear between ( or , and the next non-name char (: or , or ))
    return len(re.findall(r"[(,]([^():,]+):", s))

def parse_supports(s):
    """Extract support values: IQ-TREE format ')alrt/ufb:' OR FastTree ')support:'"""
    iq = [float(m) for m in re.findall(r"\)\d+\.?\d*/(\d+\.?\d*):", s)]
    if iq:
        return iq
    # FastTree: ')0.95:' style
    ft = [float(m) for m in re.findall(r"\)(\d+\.?\d*):", s)]
    return ft

rows = []
for d in sorted(TREES.iterdir()):
    if not d.is_dir():
        continue
    pf = d.name
    iq_tree = d / "tree.treefile"
    ft_tree = d / "tree.nwk"
    iq_log = d / "tree.iqtree"

    tool, tree_path = None, None
    if ft_tree.exists() and ft_tree.stat().st_size > 0:
        tool, tree_path = "fasttree", ft_tree
    elif iq_tree.exists() and iq_tree.stat().st_size > 0:
        tool, tree_path = "iqtree", iq_tree
    if not tool:
        rows.append({"pf": pf, "tool": "NONE", "status": "no_tree"})
        continue

    tree_s = tree_path.read_text()
    n_taxa = count_taxa_in_newick(tree_s)
    sups = parse_supports(tree_s)
    if tool == "iqtree":
        # tree.iqtree has stats
        info = iq_log.read_text() if iq_log.exists() else ""
        m_sites = re.search(r"([0-9]+) amino-acid sites", info)
        n_sites = int(m_sites.group(1)) if m_sites else None
        ufb_med = median(sups) if sups else None
        support_type = "UFB%"
    else:
        # FastTree: SH-aLRT support
        n_sites = None  # could parse from log
        ufb_med = round(median(sups), 3) if sups else None
        support_type = "SH-aLRT"

    name = fam_meta.get(pf, {}).get("name", "?")
    count = fam_meta.get(pf, {}).get("count", 0)
    rows.append({
        "pf": pf,
        "name": name[:50],
        "n_pfam_total": count,
        "tool": tool,
        "n_taxa": n_taxa,
        "n_sites": n_sites or "?",
        "support_type": support_type,
        "median_support": ufb_med if ufb_med is not None else "?",
        "status": "ok",
    })

# sort by Pfam acc
rows.sort(key=lambda r: r["pf"])

# Write CSV
out = WORK / "all_trees_summary.tsv"
with out.open("w") as fh:
    cols = ["pf", "name", "n_pfam_total", "tool", "n_taxa", "n_sites",
            "support_type", "median_support", "status"]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

# Print table
print(f"{'Pfam':<10} {'Tool':<9} {'N_taxa':>7} {'N_sites':>8} {'Support':>9} {'Name'}")
print("-" * 100)
for r in rows:
    sup = r.get("median_support", "?")
    sup_s = f"{sup:.2f}" if isinstance(sup, float) and sup <= 1.0 else str(sup)
    print(f"{r['pf']:<10} {r['tool']:<9} {r.get('n_taxa','?'):>7} {r.get('n_sites','?'):>8} {sup_s:>9} {r.get('name','')}")

# Identify problems
print("\n--- low-quality trees (< 30 sites or < 0.5 support) ---")
for r in rows:
    if r["status"] != "ok":
        continue
    sites = r.get("n_sites")
    sup = r.get("median_support")
    bad = False
    if isinstance(sites, int) and sites < 30:
        bad = True
    if isinstance(sup, float):
        if r["support_type"] == "SH-aLRT" and sup < 0.5:
            bad = True
        if r["support_type"] == "UFB%" and sup < 70:
            bad = True
    if bad:
        print(f"  {r['pf']:<10} taxa={r['n_taxa']:<4} sites={sites} sup={sup}  {r['name']}")

print(f"\n--- failed (no tree) ---")
for r in rows:
    if r["status"] != "ok":
        print(f"  {r['pf']}")
print(f"\n[*] full TSV → {out}")
print(f"[*] total ok: {sum(1 for r in rows if r['status']=='ok')}/{len(rows)}")
