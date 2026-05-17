#!/usr/bin/env python3
"""
monttinen_core.py — Monttinen 36-residue structural core via TM-align
Reference: 1VFG chain A (NTase beta-palm, CCA-adding enzyme, Aquifex aeolicus)
Queries:   3BDP-A (PolA), 1WAJ-A (PolB), 1RTD-A (RT/p66),
           1BPX-A (PolX/polbeta), 1HQM-C (RNAP-beta), 1RDR-A (RdRp/3Dpol)

CORRECTION vs original report: 1HQM chain C (beta subunit, 997 Ca) is used.
Chain B (alpha subunit, 229 Ca) was used in the original — that is wrong.
"""
import subprocess, sys, re
from pathlib import Path

try:
    from Bio.PDB import PDBParser, PDBIO, Select
    import pandas as pd
except ImportError:
    sys.exit("pip install biopython pandas")

WORKDIR = Path.home() / "monttinen_core"

STRUCTURES = [
    ("1VFG", "A", "NTase", "NTase beta-palm (CCA-adding, Aquifex)"),
    ("3BDP", "A", "PolA",  "DNA pol A (Bst pol I, B.stearothermophilus)"),
    ("1WAJ", "A", "PolB",  "DNA pol B (RB69 gp43)"),
    ("1RTD", "A", "RT",    "Reverse transcriptase p66 (HIV-1)"),
    ("1BPX", "A", "PolX",  "DNA pol X/beta (H.sapiens)"),
    ("1HQM", "C", "RNAP",  "RNA pol beta subunit (T.aquaticus) [chain C, CORRECTED]"),
    ("1RDR", "A", "RdRp",  "Viral RdRp 3Dpol (Poliovirus)"),
]


class ChainSelect(Select):
    def __init__(self, cid):
        self.cid = cid
    def accept_chain(self, c):
        return c.id == self.cid


def extract_chain(src, chain_id, dest):
    p = PDBParser(QUIET=True)
    io = PDBIO()
    io.set_structure(p.get_structure("s", str(src)))
    io.save(str(dest), ChainSelect(chain_id))


def get_residues(src, chain_id):
    p = PDBParser(QUIET=True)
    s = p.get_structure("s", str(src))
    return [(r.id[1], r.id[2].strip(), r.resname)
            for ch in s[0] if ch.id == chain_id
            for r in ch if r.id[0] == " "]


def run_tmalign(p1, p2):
    return subprocess.run(
        ["TMalign", str(p1), str(p2)],
        capture_output=True, text=True, timeout=120
    ).stdout


def parse_tmalign(out):
    """Parse TM-align v20220412 output.
    Labels are Chain_1/Chain_2; alignment is 3 long lines (no block-wrap).
    """
    lines = out.splitlines()
    tm1 = tm2 = rmsd = aligned = None
    for ln in lines:
        m = re.search(r"TM-score=\s*([\d.]+).*Chain_1", ln)
        if m:
            tm1 = float(m.group(1))
        m = re.search(r"TM-score=\s*([\d.]+).*Chain_2", ln)
        if m:
            tm2 = float(m.group(1))
        m = re.search(r"Aligned length=\s*(\d+).*RMSD=\s*([\d.]+)", ln)
        if m:
            aligned, rmsd = int(m.group(1)), float(m.group(2))

    start = next(
        (i + 1 for i, ln in enumerate(lines)
         if '":" denotes' in ln or "denotes residue pairs" in ln),
        None
    )
    pairs = []
    if start is not None and start + 2 < len(lines):
        s1  = lines[start].rstrip()
        sym = lines[start + 1].rstrip()
        s2  = lines[start + 2].rstrip()
        i1 = i2 = 0
        for pos in range(max(len(s1), len(s2))):
            c1 = s1[pos]  if pos < len(s1)  else " "
            sc = sym[pos] if pos < len(sym) else " "
            c2 = s2[pos]  if pos < len(s2)  else " "
            if c1 not in ("-", " ") and c2 not in ("-", " ") and sc in (":", "."):
                pairs.append((i1, i2))
            if c1 not in ("-", " "):
                i1 += 1
            if c2 not in ("-", " "):
                i2 += 1
    return dict(tm1=tm1, tm2=tm2, rmsd=rmsd, aligned=aligned, pairs=pairs)


COLORS = dict(
    NTase="red", PolA="orange", PolB="yellow",
    RT="forest", PolX="cyan", RNAP="blue", RdRp="magenta"
)


def write_pml(path, label, pdb_path, resnums, n_core):
    pct = 100 * len(resnums) / n_core if n_core else 0
    sel = "+".join(map(str, resnums)) if resnums else ""
    lines = [
        f"# {label}: Monttinen core  {len(resnums)}/{n_core} ({pct:.0f}%)",
        f"load {pdb_path}, {label}",
        "color gray70, all",
        "hide everything, all",
        "show cartoon, all",
    ]
    if sel:
        lines += [
            f"color {COLORS.get(label, 'red')}, {label} and resi {sel}",
            f"show sticks, {label} and resi {sel}",
        ]
    lines += ["bg_color white", f"zoom {label}"]
    Path(path).write_text("\n".join(lines))


def main():
    print("=" * 60)
    print("  Monttinen structural core — TM-align")
    print("=" * 60)

    # 1. Extract chains
    print("\n[1] Extracting chains ...")
    pdbs, residues = {}, {}
    for pdb_id, chain, label, desc in STRUCTURES:
        dest = WORKDIR / f"{pdb_id}_{chain}.pdb"
        if not dest.exists():
            extract_chain(WORKDIR / f"{pdb_id}.pdb", chain, dest)
        pdbs[label]     = dest
        residues[label] = get_residues(WORKDIR / f"{pdb_id}.pdb", chain)
        print(f"  {label:6s}  {pdb_id}-{chain}  {len(residues[label]):4d} res  {desc}")

    ref     = "NTase"
    ref_res = residues[ref]
    queries = [s[2] for s in STRUCTURES[1:]]

    # 2. TM-align
    print(f"\n[2] TM-align {ref} vs each query ...")
    all_pairs, tm_rows = {}, []
    for q in queries:
        raw = run_tmalign(pdbs[ref], pdbs[q])
        p   = parse_tmalign(raw)
        all_pairs[q] = p["pairs"]
        tm_rows.append(dict(
            fold=q, tm_ref=p["tm1"], tm_tgt=p["tm2"],
            rmsd=p["rmsd"], n_aligned=p["aligned"], n_pairs=len(p["pairs"])
        ))
        print(f"  {ref} vs {q:6s}  TM(ref)={p['tm1']:.3f}  TM(tgt)={p['tm2']:.3f}"
              f"  RMSD={p['rmsd']:.2f}  aligned={p['aligned']}  pairs={len(p['pairs'])}")

    # 3. Core
    print("\n[3] Core sizes at multiple thresholds:")
    n_ref    = len(ref_res)
    ref_sets = {q: {r for r, _ in all_pairs[q]} for q in queries}
    covered  = {idx: {q for q in queries if idx in ref_sets[q]} for idx in range(n_ref)}

    for thr in [6, 5, 4, 3]:
        n = sum(1 for v in covered.values() if len(v) >= thr)
        print(f"  >= {thr}/6 folds : {n:3d} residues")

    core5 = sorted(i for i, v in covered.items() if len(v) >= 5)
    print(f"\n  Working core (>=5/6): {len(core5)} res  [Monttinen: 36 res, 89 structs]")

    # 4. Presence / absence
    print("\n[4] Presence/absence per fold:")
    pa_rows = [dict(fold=ref + "(ref)", present=len(core5), total=len(core5),
                    pct=100.0, absent_resnums_1VFG="")]
    for q in queries:
        present = [i for i in core5 if i in ref_sets[q]]
        absent  = [ref_res[i][0] for i in core5 if i not in ref_sets[q] and i < len(ref_res)]
        pct     = 100 * len(present) / len(core5) if core5 else 0
        status  = "PRESENT" if pct >= 70 else ("PARTIAL" if pct >= 40 else "ABSENT")
        pa_rows.append(dict(
            fold=q, present=len(present), total=len(core5), pct=round(pct, 1),
            absent_resnums_1VFG=(";".join(map(str, absent[:15])) +
                                  ("..." if len(absent) > 15 else ""))
        ))
        print(f"  {q:6s}: {len(present):2d}/{len(core5)}  ({pct:5.1f}%)  {status}")

    # 5. Core residue table on 1VFG
    core_rows = [
        dict(idx=i, resnum=ref_res[i][0], resname=ref_res[i][2],
             n_folds=len(covered[i]), folds=";".join(sorted(covered[i])))
        for i in core5 if i < len(ref_res)
    ]

    # 6. Write files
    print("\n[5] Writing outputs ...")
    pd.DataFrame(tm_rows).to_csv(WORKDIR / "tmscore_summary.tsv",   sep="\t", index=False)
    pd.DataFrame(pa_rows).to_csv(WORKDIR / "presence_absence.tsv",  sep="\t", index=False)
    pd.DataFrame(core_rows).to_csv(WORKDIR / "core_residues.tsv",   sep="\t", index=False)
    print("  tmscore_summary.tsv   presence_absence.tsv   core_residues.tsv")

    # Combined PyMOL script
    lines = [f"# Monttinen core — all 7 folds  core={len(core5)} res (>=5/6)"]
    for lbl, path in pdbs.items():
        lines.append(f"load {path}, {lbl}")
    lines += ["color gray70, all", "hide everything, all", "show cartoon, all", ""]
    for lbl in [ref] + queries:
        tgt_res = residues[lbl]
        if lbl == ref:
            rn = [ref_res[i][0] for i in core5 if i < len(ref_res)]
        else:
            r2t = {r: t for r, t in all_pairs[lbl]}
            rn  = [tgt_res[r2t[i]][0] for i in core5 if i in r2t and r2t[i] < len(tgt_res)]
        if rn:
            sel = "+".join(map(str, rn))
            lines += [f"color {COLORS[lbl]}, {lbl} and resi {sel}",
                      f"show sticks, {lbl} and resi {sel}"]
    lines += ["bg_color white", "zoom NTase"]
    (WORKDIR / "color_core.pml").write_text("\n".join(lines))

    for lbl in [ref] + queries:
        tgt_res = residues[lbl]
        if lbl == ref:
            rn = [ref_res[i][0] for i in core5 if i < len(ref_res)]
        else:
            r2t = {r: t for r, t in all_pairs[lbl]}
            rn  = [tgt_res[r2t[i]][0] for i in core5 if i in r2t and r2t[i] < len(tgt_res)]
        write_pml(WORKDIR / f"color_{lbl}.pml", lbl, pdbs[lbl], rn, len(core5))
    print("  color_core.pml   color_<fold>.pml (one per structure)")

    # Final summary
    print("\n" + "=" * 60)
    print(f"  {'Fold':<12} {'Present':>8} {'%':>7}  Result")
    print(f"  {'-'*12} {'-'*8} {'-'*7}  {'-'*20}")
    for r in pa_rows:
        tag = ("reference" if "(ref)" in r["fold"] else
               "PRESENT" if r["pct"] >= 70 else
               "PARTIAL" if r["pct"] >= 40 else "ABSENT")
        print(f"  {r['fold']:<12} {r['present']:>4}/{r['total']:<4}  {r['pct']:>5.1f}%  {tag}")
    print()
    print(f"  Core (>=5/6):   {len(core5)} residues on 1VFG chain A")
    print(f"  Monttinen 2016: 36 residues (5 superfamilies, 89 experimental structures)")
    print(f"  Key fix:        1HQM-C (beta subunit) replaces 1HQM-B (alpha) from report")
    print(f"\n  Outputs: {WORKDIR}")


if __name__ == "__main__":
    main()
