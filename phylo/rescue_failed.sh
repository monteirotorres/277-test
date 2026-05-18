#!/usr/bin/env bash
# Rescue PF00978 and PF02661 by running FastTree directly on the untrimmed MAFFT alignment.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phylo277
for pf in PF00978 PF02661; do
    d=/home/torres/phylo277/trees/$pf
    aln=$d/combined.aln
    out=$d/tree.nwk
    if [ ! -f "$aln" ]; then
        echo "[$pf] missing $aln — skipping"
        continue
    fi
    # Strip all-gap rows (which trimAl would remove)
    python3 -c "
from pathlib import Path
seqs={}
acc=None;buf=[]
for ln in open('$aln'):
    ln=ln.rstrip()
    if ln.startswith('>'):
        if acc: seqs[acc]=''.join(buf)
        acc=ln[1:].split()[0]; buf=[]
    else: buf.append(ln)
seqs[acc]=''.join(buf)
seqs={k:v for k,v in seqs.items() if v.replace('-','').strip()}
with open('$d/combined.untrimmed.fasta','w') as f:
    for k,v in seqs.items(): f.write(f'>{k}\n{v}\n')
print(f'kept {len(seqs)} sequences')
"
    echo "[$pf] running FastTree on untrimmed alignment..."
    fasttree -lg -gamma -spr 4 -mlacc 2 -slownni $d/combined.untrimmed.fasta > $out 2>$d/fasttree.stderr
    if [ -s "$out" ]; then
        echo "[$pf] OK → $out"
    else
        echo "[$pf] FAILED"
        cat $d/fasttree.stderr
    fi
done
