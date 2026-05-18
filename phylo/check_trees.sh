#!/usr/bin/env bash
echo "=== tree directory contents ==="
for d in /home/torres/phylo277/trees/*/; do
    pf=$(basename "$d")
    n=$(ls "$d" 2>/dev/null | wc -l)
    has_tree="NO"
    [ -f "$d/tree.treefile" ] && has_tree="YES"
    has_iqlog="NO"
    [ -f "$d/tree.iqtree" ] && has_iqlog="YES"
    echo "$pf [files=$n, treefile=$has_tree, iqtree=$has_iqlog]"
done

echo
echo "=== PF02696 failure log (last 25 lines) ==="
tail -25 /home/torres/phylo277/logs/PF02696.log 2>/dev/null || echo "(no log)"

echo
echo "=== example: PF00078 (RT) iqtree summary ==="
if [ -f /home/torres/phylo277/trees/PF00078/tree.iqtree ]; then
    grep -E "Number of |Log-likelihood|BIC|Best-fit model|Total tree length|^Sum of internal" /home/torres/phylo277/trees/PF00078/tree.iqtree | head -20
fi
