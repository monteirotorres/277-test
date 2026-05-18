#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phylo277
echo "=== which ==="
which fasttree FastTree veryfasttree 2>&1
echo "=== ls bins matching tree ==="
ls ~/miniconda3/envs/phylo277/bin/ | grep -iE 'tree|fast'
