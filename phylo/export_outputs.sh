#!/usr/bin/env bash
# Copy pipeline outputs from WSL to the Windows workspace.
src=/home/torres/phylo277
dst=/mnt/c/Users/monte/OneDrive/Documentos/claude_work/phylo/results
mkdir -p $dst/trees
echo "[*] copying top-level files..."
cp $src/ec277_swissprot.tsv         $dst/
cp $src/family_summary.tsv          $dst/
cp $src/protein_to_family.tsv       $dst/
cp $src/pfam_names.json             $dst/
cp $src/all_trees_summary.tsv       $dst/
cp $src/tree_build_summary_fasttree.json $dst/ 2>/dev/null
cp $src/tree_build_summary.json     $dst/iqtree_first_batch_summary.json 2>/dev/null
echo "[*] copying tree files..."
for d in $src/trees/*/; do
    pf=$(basename $d)
    mkdir -p $dst/trees/$pf
    for f in tree.treefile tree.nwk tree.iqtree combined.aln combined.trim; do
        [ -f "$d/$f" ] && cp "$d/$f" $dst/trees/$pf/
    done
done
echo "[*] total trees: $(ls $dst/trees | wc -l)"
echo "[*] outputs dir size: $(du -sh $dst | awk '{print $1}')"
