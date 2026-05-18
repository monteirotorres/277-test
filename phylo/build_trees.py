#!/usr/bin/env python3
"""
build_trees.py — per-family ML trees with cross-family outgroups.

For each Pfam family with ≥ MIN_SIZE members in ~/phylo277/families/:
  1. CD-HIT to collapse near-duplicates (≥ CDHIT_ID) and cap at CAP_SIZE
  2. Pick one outgroup per OTHER family (median-length sequence)
  3. MAFFT alignment
  4. trimAl --gappyout
  5. IQ-TREE LG+G4 + UFBoot 1000 + aLRT 1000

Run with: python3 build_trees.py [--only PFxxx[,PFyyy...]] [--threads N]
"""
import argparse, json, os, random, subprocess, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

WORKDIR  = Path.home() / "phylo277"
FAM_DIR  = WORKDIR / "families"
TREE_DIR = WORKDIR / "trees"
TREE_DIR.mkdir(exist_ok=True)
LOG_DIR  = WORKDIR / "logs"
LOG_DIR.mkdir(exist_ok=True)

MIN_SIZE  = 30
CAP_SIZE  = 120
CDHIT_ID  = 0.80
RANDOM_SEED = 42

random.seed(RANDOM_SEED)


def read_fasta(path):
    seqs = {}
    acc = None
    buf = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc = ln[1:].split()[0]
                buf = []
            else:
                buf.append(ln)
        if acc:
            seqs[acc] = "".join(buf)
    return seqs


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")


def run(cmd, log_path, env=None, timeout=None):
    """Run a command, capture stdout+stderr to log."""
    with open(log_path, "ab") as log:
        log.write(("\n$ " + " ".join(str(c) for c in cmd) + "\n").encode())
        log.flush()
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            env=env, timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError(f"FAILED: {' '.join(str(c) for c in cmd)} (see {log_path})")


def cdhit_reduce(infa, outfa, log_path, identity=CDHIT_ID, cap=CAP_SIZE):
    seqs = read_fasta(infa)
    if len(seqs) <= cap:
        # Just copy
        write_fasta(outfa, seqs)
        return list(seqs.keys())
    # Pick a word size compatible with identity (cd-hit rules of thumb)
    n = 5 if identity >= 0.7 else 4 if identity >= 0.6 else 3 if identity >= 0.5 else 2
    cmd = ["cd-hit", "-i", str(infa), "-o", str(outfa),
           "-c", str(identity), "-n", str(n),
           "-M", "0", "-T", "1", "-d", "0"]
    run(cmd, log_path)
    reduced = read_fasta(outfa)
    if len(reduced) > cap:
        # random downsample to cap
        keys = list(reduced.keys())
        random.shuffle(keys)
        keep = set(keys[:cap])
        reduced = {k: v for k, v in reduced.items() if k in keep}
        write_fasta(outfa, reduced)
    return list(reduced.keys())


def median_length_seq(seqs):
    """Return (acc, seq) with median length."""
    items = sorted(seqs.items(), key=lambda kv: len(kv[1]))
    return items[len(items) // 2]


def process_family(pf, all_pfs, env_path):
    """Worker for one focal family."""
    t0 = time.time()
    log_path = LOG_DIR / f"{pf}.log"
    log_path.write_bytes(b"")  # reset

    src_fa  = FAM_DIR / f"{pf}.fasta"
    work    = TREE_DIR / pf
    work.mkdir(exist_ok=True)
    focal_fa = work / "focal.fasta"
    focal_red = work / "focal_reduced.fasta"
    combined = work / "combined.fasta"
    aln      = work / "combined.aln"
    trimmed  = work / "combined.trim"

    # 1. CD-HIT reduce focal
    cdhit_reduce(src_fa, focal_red, log_path)
    focal_seqs = read_fasta(focal_red)

    # 2. Outgroups: median-length entry from each other family
    out_seqs = {}
    for op in all_pfs:
        if op == pf:
            continue
        op_seqs = read_fasta(FAM_DIR / f"{op}.fasta")
        if not op_seqs:
            continue
        acc, sq = median_length_seq(op_seqs)
        out_seqs[f"OG_{op}_{acc}"] = sq

    # 3. Combined FASTA
    write_fasta(combined, {**focal_seqs, **out_seqs})

    n_focal = len(focal_seqs)
    n_out   = len(out_seqs)
    print(f"[{pf}] focal={n_focal} out={n_out} → MAFFT", flush=True)

    # 4. MAFFT
    mafft = env_path / "bin" / "mafft"
    with open(aln, "w") as out:
        sub = subprocess.run([str(mafft), "--auto", "--thread", "1", str(combined)],
                              stdout=out, stderr=subprocess.PIPE)
    if sub.returncode != 0:
        log_path.write_bytes(log_path.read_bytes() + sub.stderr)
        raise RuntimeError(f"MAFFT failed for {pf}")

    # 5. trimAl
    trimal = env_path / "bin" / "trimal"
    run([str(trimal), "-in", str(aln), "-out", str(trimmed), "-gappyout"], log_path)

    # 6. IQ-TREE
    iqtree = env_path / "bin" / "iqtree"
    pre = work / "tree"
    cmd = [str(iqtree), "-s", str(trimmed), "-m", "LG+G4",
           "--ufboot", "1000", "--alrt", "1000",
           "-T", "2", "--prefix", str(pre), "-redo",
           "-quiet"]
    run(cmd, log_path)

    dt = time.time() - t0
    print(f"[{pf}] DONE in {dt:.0f}s · focal={n_focal} out={n_out}", flush=True)
    return {"pf": pf, "n_focal": n_focal, "n_out": n_out, "elapsed_s": dt}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma-separated Pfam list")
    ap.add_argument("--threads", "-j", type=int, default=4)
    ap.add_argument("--min-size", type=int, default=MIN_SIZE)
    args = ap.parse_args()

    env_path = Path(os.environ.get("CONDA_PREFIX") or
                    "/home/torres/miniconda3/envs/phylo277")

    # Determine eligible families
    fam_summary = (WORKDIR / "family_summary.tsv").read_text().strip().split("\n")
    eligible = []
    for ln in fam_summary[1:]:
        cells = ln.split("\t")
        if int(cells[2]) >= args.min_size:
            eligible.append(cells[0])
    print(f"[*] {len(eligible)} families with ≥ {args.min_size} members")

    if args.only:
        wanted = set(args.only.split(","))
        to_process = [p for p in eligible if p in wanted]
    else:
        to_process = eligible

    print(f"[*] processing {len(to_process)} families, {args.threads} parallel\n")

    results = []
    errors  = []
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_family, p, eligible, env_path): p
                   for p in to_process}
        for fut in as_completed(futures):
            pf = futures[fut]
            try:
                results.append(fut.result())
            except Exception as e:
                print(f"[{pf}] ERROR: {e}", flush=True)
                errors.append({"pf": pf, "err": str(e)})

    # Save summary
    summary = {"results": results, "errors": errors,
               "params": {"min_size": args.min_size, "cap": CAP_SIZE,
                          "cdhit": CDHIT_ID, "threads": args.threads}}
    (WORKDIR / "tree_build_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[done] {len(results)} OK, {len(errors)} ERR")
    print(f"[done] summary → {WORKDIR}/tree_build_summary.json")


if __name__ == "__main__":
    main()
