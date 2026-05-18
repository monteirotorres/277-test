#!/usr/bin/env python3
"""
build_trees_v2.py — fixed + faster per-family ML trees.

Changes vs v1:
  * MAFFT: --anysymbol  (handles selenocysteine 'U')
  * trimAl: -automated1 (preserves divergent regions vs gappyout)
  * Optional --tool fasttree   (FastTree -lg -gamma + SH-aLRT)  ~50x faster
  * Skip families that already have a tree.treefile / tree.nwk

Usage:
    python3 build_trees_v2.py --tool fasttree --threads 4
    python3 build_trees_v2.py --tool iqtree --only PF02696,PF00680 --force
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
    acc, buf = None, []
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
    with open(log_path, "ab") as log:
        log.write(("\n$ " + " ".join(str(c) for c in cmd) + "\n").encode())
        log.flush()
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                            env=env, timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError(f"FAILED ({p.returncode}): {' '.join(str(c) for c in cmd)}")


def cdhit_reduce(infa, outfa, log_path, identity=CDHIT_ID, cap=CAP_SIZE):
    seqs = read_fasta(infa)
    if len(seqs) <= cap:
        write_fasta(outfa, seqs)
        return
    n = 5 if identity >= 0.7 else 4 if identity >= 0.6 else 3 if identity >= 0.5 else 2
    cmd = ["cd-hit", "-i", str(infa), "-o", str(outfa),
           "-c", str(identity), "-n", str(n),
           "-M", "0", "-T", "1", "-d", "0"]
    run(cmd, log_path)
    reduced = read_fasta(outfa)
    if len(reduced) > cap:
        keys = list(reduced.keys())
        random.shuffle(keys)
        keep = set(keys[:cap])
        reduced = {k: v for k, v in reduced.items() if k in keep}
        write_fasta(outfa, reduced)


def median_length_seq(seqs):
    items = sorted(seqs.items(), key=lambda kv: len(kv[1]))
    return items[len(items) // 2]


def process_family(pf, all_pfs, env_path, tool, force):
    t0 = time.time()
    log_path = LOG_DIR / f"{pf}.log"
    work = TREE_DIR / pf
    work.mkdir(exist_ok=True)

    # Skip if already done
    iq_tree = work / "tree.treefile"
    ft_tree = work / "tree.nwk"
    if not force:
        if tool == "iqtree" and iq_tree.exists() and iq_tree.stat().st_size > 0:
            return {"pf": pf, "skipped": True, "tool": tool}
        if tool == "fasttree" and ft_tree.exists() and ft_tree.stat().st_size > 0:
            return {"pf": pf, "skipped": True, "tool": tool}

    log_path.write_bytes(b"")  # reset

    src_fa    = FAM_DIR / f"{pf}.fasta"
    focal_red = work / "focal_reduced.fasta"
    combined  = work / "combined.fasta"
    aln       = work / "combined.aln"
    trimmed   = work / "combined.trim"

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
    print(f"[{pf}] focal={n_focal} out={n_out} → MAFFT (--anysymbol)", flush=True)

    # 4. MAFFT --anysymbol
    mafft = env_path / "bin" / "mafft"
    with open(aln, "w") as out:
        sub = subprocess.run([str(mafft), "--auto", "--anysymbol",
                              "--thread", "1", str(combined)],
                              stdout=out, stderr=subprocess.PIPE)
    if sub.returncode != 0:
        with open(log_path, "ab") as f:
            f.write(b"\n=== MAFFT stderr ===\n" + sub.stderr)
        raise RuntimeError(f"MAFFT failed for {pf}")

    # 5. trimAl -automated1
    trimal = env_path / "bin" / "trimal"
    run([str(trimal), "-in", str(aln), "-out", str(trimmed),
         "-automated1"], log_path)

    # Sanity: count sites
    first = next(iter(read_fasta(trimmed).values()))
    n_sites = len(first)
    print(f"[{pf}] trimmed alignment: {n_sites} sites", flush=True)

    # 6. Build tree
    if tool == "iqtree":
        iqtree = env_path / "bin" / "iqtree"
        pre = work / "tree"
        cmd = [str(iqtree), "-s", str(trimmed), "-m", "LG+G4",
               "--ufboot", "1000", "--alrt", "1000",
               "-T", "2", "--prefix", str(pre), "-redo", "-quiet"]
        run(cmd, log_path)
    else:  # fasttree
        ft = env_path / "bin" / "fasttree"
        out_tree = work / "tree.nwk"
        # FastTree with LG + gamma + SH-aLRT (default support)
        with open(out_tree, "w") as o:
            sub = subprocess.run([str(ft), "-lg", "-gamma",
                                  "-spr", "4", "-mlacc", "2", "-slownni",
                                  str(trimmed)],
                                 stdout=o, stderr=subprocess.PIPE)
        if sub.returncode != 0:
            with open(log_path, "ab") as f:
                f.write(b"\n=== FastTree stderr ===\n" + sub.stderr)
            raise RuntimeError(f"FastTree failed for {pf}")
        with open(log_path, "ab") as f:
            f.write(b"\n=== FastTree stderr ===\n" + sub.stderr)

    dt = time.time() - t0
    print(f"[{pf}] DONE [{tool}] in {dt:.0f}s · focal={n_focal} out={n_out} sites={n_sites}", flush=True)
    return {"pf": pf, "n_focal": n_focal, "n_out": n_out,
            "elapsed_s": round(dt, 1), "n_sites": n_sites, "tool": tool}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", choices=["iqtree", "fasttree"], default="fasttree")
    ap.add_argument("--only", default=None, help="comma-separated Pfam list")
    ap.add_argument("--skip", default=None, help="comma-separated Pfam list to exclude")
    ap.add_argument("--threads", "-j", type=int, default=4)
    ap.add_argument("--min-size", type=int, default=MIN_SIZE)
    ap.add_argument("--force", action="store_true", help="rebuild even if tree exists")
    args = ap.parse_args()

    env_path = Path(os.environ.get("CONDA_PREFIX") or
                    "/home/torres/miniconda3/envs/phylo277")

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
        to_process = list(eligible)

    if args.skip:
        skipset = set(args.skip.split(","))
        to_process = [p for p in to_process if p not in skipset]

    print(f"[*] tool={args.tool}, processing {len(to_process)} families, "
          f"threads={args.threads}, force={args.force}\n")

    results, errors = [], []
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_family, p, eligible, env_path,
                              args.tool, args.force): p
                   for p in to_process}
        for fut in as_completed(futures):
            pf = futures[fut]
            try:
                r = fut.result()
                results.append(r)
                if r.get("skipped"):
                    print(f"[{pf}] (skipped — already built)", flush=True)
            except Exception as e:
                print(f"[{pf}] ERROR: {e}", flush=True)
                errors.append({"pf": pf, "err": str(e)})

    summary_path = WORKDIR / f"tree_build_summary_{args.tool}.json"
    summary_path.write_text(json.dumps(
        {"results": results, "errors": errors,
         "params": {"tool": args.tool, "min_size": args.min_size,
                    "cap": CAP_SIZE, "cdhit": CDHIT_ID,
                    "threads": args.threads}},
        indent=2))
    n_done = len([r for r in results if not r.get("skipped")])
    n_skip = len([r for r in results if r.get("skipped")])
    print(f"\n[done] {n_done} built, {n_skip} skipped, {len(errors)} errors")
    print(f"[done] summary → {summary_path}")


if __name__ == "__main__":
    main()
