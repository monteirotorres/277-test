#!/usr/bin/env bash
# Copy all Newick trees + a JSON manifest to claude_work/trees/ for the website.
src=/home/torres/phylo277
dst=/mnt/c/Users/monte/OneDrive/Documentos/claude_work/trees
rm -rf $dst
mkdir -p $dst

# Copy newick files, prefer IQ-TREE output where present
for d in $src/trees/*/; do
    pf=$(basename $d)
    if [ -f "$d/tree.treefile" ]; then
        cp "$d/tree.treefile" "$dst/$pf.nwk"
    elif [ -f "$d/tree.nwk" ]; then
        cp "$d/tree.nwk" "$dst/$pf.nwk"
    fi
done

# Build a manifest JSON from all_trees_summary.tsv + family_summary.tsv + pfam_names.json
python3 <<'PY'
import json, csv, os
work = '/home/torres/phylo277'
dst = '/mnt/c/Users/monte/OneDrive/Documentos/claude_work/trees'

# load family meta
fam = {}
with open(f'{work}/family_summary.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        fam[r['pfam_acc']] = {
            'name': r['pfam_name'],
            'count': int(r['count']),
            'sample_member': r['sample_member'],
            'sample_ec': r['sample_ec'],
        }

# load tree summary
trees = []
with open(f'{work}/all_trees_summary.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        pf = r['pf']
        sup = r.get('median_support', '?')
        sup_n = None
        if sup not in ('?', '', None):
            try:
                sup_n = float(sup)
            except: pass
        meta = fam.get(pf, {})
        trees.append({
            'pf': pf,
            'name': meta.get('name', r.get('name', '?')),
            'count_total': meta.get('count', 0),
            'tool': r.get('tool'),
            'n_taxa': int(r['n_taxa']) if r['n_taxa'].isdigit() else None,
            'n_sites': int(r['n_sites']) if r['n_sites'].isdigit() else None,
            'support_type': r.get('support_type'),
            'support': sup_n,
            'sample_member': meta.get('sample_member', ''),
            'sample_ec': meta.get('sample_ec', ''),
        })

trees.sort(key=lambda t: -t['count_total'])
manifest = {
    'generated': '2026-05-17',
    'pipeline': {
        'source': 'UniProt Swiss-Prot (reviewed:true) AND ec:2.7.7.*',
        'total_entries': 13886,
        'distinct_pfams': 647,
        'families_total': 88,
        'families_ge_30': 56,
        'tools': {
            'mafft': '--auto --anysymbol',
            'trimal': '-automated1',
            'iqtree': 'LG+G4 + UFBoot 1000 + aLRT 1000 (9 canonical families)',
            'fasttree': '-lg -gamma -spr 4 -mlacc 2 -slownni (47 others)',
        },
    },
    'trees': trees,
}
with open(f'{dst}/manifest.json', 'w') as f:
    json.dump(manifest, f, indent=2)
print(f'manifest: {len(trees)} trees')
PY

echo "[*] copied $(ls $dst/*.nwk | wc -l) trees"
echo "[*] dst size: $(du -sh $dst | awk '{print $1}')"
ls $dst | head -5
