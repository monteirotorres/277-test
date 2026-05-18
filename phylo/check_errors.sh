#!/usr/bin/env bash
for pf in PF00978 PF02661; do
    echo "================ $pf ================"
    echo "--- log tail ---"
    tail -30 /home/torres/phylo277/logs/$pf.log 2>/dev/null
    echo "--- dir contents ---"
    ls -la /home/torres/phylo277/trees/$pf/ 2>/dev/null
    echo
done
