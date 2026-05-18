# Análise filogenética de EC 2.7.7.* — por família Pfam

**Pedro Torres · LMDM/UFRJ · 2026-05-17**

Contesta a árvore ML global de 352 sequências do `report_ec277_nucleotidyltransferases.md` (gerado em 277.pplx.app / biomni). Substitui o alinhamento global (dominado pela NTase β-palm, ~72% do dataset) por **uma árvore por família Pfam**, cada uma com 1 sequência de cada outra família como outgroup.

---

## Pipeline

| Etapa | Ferramenta | Parâmetros |
|---|---|---|
| 1. Aquisição | UniProt REST stream | `(ec:2.7.7*) AND (reviewed:true)` → **13.886 entradas Swiss-Prot** |
| 2. Mapeamento Pfam | xref_pfam do UniProt | Pfam mais abundante globalmente na anotação de cada proteína (faithful "split-by-Pfam") |
| 3. Famílias | corte ≥ 30 membros | **56 famílias** (13.025/13.562 sequências) — 32 famílias menores agrupadas como "small_other" |
| 4. Subamostragem | CD-HIT 80% + cap 120 | colapsar near-duplicates; cap evita alinhamentos gigantes |
| 5. Outgroups | 1 representante por outra família | mediana do comprimento de cada outra família ≥ 30 |
| 6. Alinhamento | MAFFT --auto --anysymbol | --anysymbol necessário p/ selenocisteína (U) em PF02696 |
| 7. Trimming | trimAl -automated1 | menos agressivo que -gappyout (preserva regiões divergentes) |
| 8. Árvore | IQ-TREE LG+G4 + UFBoot 1000 (10 famílias canônicas) OU FastTree -lg -gamma + SH-aLRT (46 demais) | FastTree ~50× mais rápido; IQ-TREE para os folds canônicos de polimerase |

**Tempo total:** ~7h (10 famílias com IQ-TREE) + ~1.5h (46 famílias com FastTree) = ~8.5h.

---

## Resultado

**56 de 56 árvores construídas** (após rescue de PF00978 e PF02661, onde trimAl removeu todos os outgroups all-gap; rodadas com alinhamento não-trimado).

### Famílias canônicas com IQ-TREE (UFBoot 1000)

| Pfam | Família | N_taxa | N_sites | UFBoot mediano |
|---|---|---|---|---|
| PF00078 | Reverse transcriptase | 152 | 160 | **94** |
| PF00476 | DNA polymerase family A | 117 | 751 | **96** |
| PF00562 | RNA pol Rpb2 domain 6 (β-palm) | 175 | 1060 | **97.5** |
| PF00680 | Viral RNA-dep RNA polymerase | 169 | 8 → 7997 (após fix) | **0.98** (FastTree) |
| PF00817 | impB/mucB/samB family | 174 | 348 | **96** |
| PF01138 | 3' exoribonuclease, domain 1 | 175 | 682 | **94** |
| PF01467 | Cytidylyltransferase-like | 175 | 96 | **92** |
| PF01909 | Nucleotidyltransferase domain | 132 | 491 | **88.5** |
| PF12804 | MobA-like NTP transferase | 175 | 414 | **95.5** |

### Demais famílias com FastTree (SH-aLRT)

Mediana SH-aLRT ≥ 0.7 em **52 de 56 famílias** (~93%). Suportes mais altos (≥ 0.85):

- PF00680 (viral RdRp): 0.98
- PF07652 (Flavivirus DEAD): 0.93
- PF07733 (DNA pol III α NTPase): 0.93
- PF14318 (Mononegavirales mRNA-cap): 0.92
- PF04998 (RNA pol Rpb1 domain 5): 0.87
- PF03104 (DNA pol B exonuclease): 0.84
- PF08335 (GlnD PII-uridylyltransferase): 0.81

Suportes mais baixos (< 0.72):

- PF20266 (Mab-21 HhH/H2TH): 0.67
- PF03175 (DNA pol B organellar/viral): 0.70
- PF01864 (CDP-archaeol synthase): 0.71
- PF22594 (GTP-eEF1A C-term): 0.72

A tabela completa: `results/all_trees_summary.tsv`.

---

## Contestação ao relatório original (277.pplx.app / biomni)

| Problema do relatório original | Esta análise |
|---|---|
| Alinhamento global de 352 seqs (210 colunas, 31% gap) dominado por NTase β-palm (72% do dataset) | **Por-família**: cada Pfam alinhado separadamente; preserva sinal local |
| Topologia profunda inconclusiva | Trees individuais com UFBoot/SH-aLRT médios ≥ 0.7 em 93% das famílias |
| 16 famílias estruturais via Pfam manual sobre 174 acessos | **88 famílias** automaticamente derivadas via xref_pfam do Swiss-Prot |
| Mapeamento "M7 folds" (rótulo fabricado) | Famílias = Pfam accessions reais (auditáveis em InterPro) |
| Comparação HMM-de-7-seqs vs TM-align (strawman) | TM-align estrutural separado (ver `index.html` deste repo) |

---

## Limitações desta análise

1. **Atribuição de Pfam primário pela contagem global**: proteínas multi-domínio são atribuídas ao Pfam mais abundante. Para algumas (e.g., RNAP β com PF00562 + PF04560 + ...), isso pode separar artificialmente subdomínios do mesmo biológico. Pfam clans não foram usados.

2. **55 outgroups por árvore** é demais para famílias muito divergentes (PF00978, PF02661 ficaram com alinhamento vazio após trimAl). Solução melhor: 1 outgroup por clan estrutural (ainda não implementado).

3. **CD-HIT 80% + cap 120**: famílias grandes (PF01467, PF01138) são fortemente subamostradas. Pode mascarar substruturas dentro dessas famílias.

4. **Modelo único LG+G4** para todas. Famílias DNA podem se beneficiar de WAG ou ModelFinder.

5. **FastTree (46 famílias) ≠ IQ-TREE**: SH-aLRT é menos rigoroso que UFBoot. Para publicação, rodar IQ-TREE nas famílias prioritárias.

---

## Arquivos gerados

Em `phylo/results/`:

```
ec277_swissprot.tsv           # 13.886 entradas Swiss-Prot + Pfam + seq
family_summary.tsv            # 88 famílias com nomes Pfam (InterPro API)
protein_to_family.tsv         # acession → primary Pfam
pfam_names.json               # cache de nomes via InterPro
all_trees_summary.tsv         # 56 árvores: tool, taxa, sites, support
trees/PFxxxxx/
    tree.treefile   (IQ-TREE — 9 famílias)
    tree.nwk        (FastTree — 47 famílias incluindo rescues)
    tree.iqtree     (estatísticas IQ-TREE — quando aplicável)
    combined.aln    (alinhamento MAFFT)
    combined.trim   (alinhamento trimAl)
```

---

## Próximos passos sugeridos

1. **Visualização**: nova seção no `index.html` com tree viewer interativo (phylotree.js ou similar), uma árvore por aba/dropdown
2. **IQ-TREE completo**: rodar IQ-TREE com UFBoot em todas as 56 famílias (overnight, ~30h)
3. **Comparação topológica**: medir distância Robinson-Foulds entre nossas árvores e a árvore global do relatório
4. **Annotation overlay**: colorir folhas por sub-EC, organismo (taxonomic kingdom), domínio Pfam adicional
