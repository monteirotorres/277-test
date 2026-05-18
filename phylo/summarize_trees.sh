#!/usr/bin/env bash
# Summarize per-family tree quality.
echo -e "pfam\tn_taxa\tn_sites\tn_pars_inf\tloglik\ttree_len\tmedian_ufb"
for d in /home/torres/phylo277/trees/*/; do
    pf=$(basename "$d")
    iq="$d/tree.iqtree"
    nw="$d/tree.treefile"
    [ -f "$iq" ] || continue
    n_taxa=$(grep "Input data:" "$iq" | grep -oE "[0-9]+ sequences" | head -1 | awk '{print $1}')
    n_sites=$(grep "Input data:" "$iq" | grep -oE "[0-9]+ amino-acid sites" | head -1 | awk '{print $1}')
    n_pars=$(grep -oE "Number of parsimony informative sites: [0-9]+" "$iq" | head -1 | awk '{print $NF}')
    loglik=$(grep "BEST SCORE FOUND" "$iq" | head -1 | grep -oE "\-?[0-9.]+$" | head -1)
    [ -z "$loglik" ] && loglik=$(grep "Log-likelihood of the tree:" "$iq" | head -1 | grep -oE "\-[0-9.]+" | head -1)
    tlen=$(grep "Total tree length" "$iq" | head -1 | grep -oE "[0-9]+\.[0-9]+" | head -1)
    # UFBoot median from tree string: bootstrap values in tree.treefile are like ")value/value:" — first number after )
    if [ -f "$nw" ]; then
        median_ufb=$(grep -oE "\)[0-9.]+/[0-9.]+" "$nw" | sed 's|^)||' | awk -F'/' '{print $2}' | sort -n | awk '
            { a[NR]=$1 }
            END {
                if (NR==0) {print "NA"; exit}
                if (NR%2) print a[(NR+1)/2]; else printf "%.1f\n", (a[NR/2]+a[NR/2+1])/2
            }')
    else
        median_ufb="NA"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$pf" "${n_taxa:-?}" "${n_sites:-?}" "${n_pars:-?}" "${loglik:-?}" "${tlen:-?}" "${median_ufb:-?}"
done
