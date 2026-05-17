# 277-test — Mönttinen Structural Core via TM-align

Análise local do core estrutural de Mönttinen et al. 2016 nas polimerases de ácido nucleico (EC 2.7.7.\*), com correção crítica do erro de cadeia do relatório original.

**Autor:** Pedro Torres — LMDM/UFRJ
**Data:** 2026-05-17

---

## Resumo

Refazer e contestar análises do `report_ec277_nucleotidyltransferases.md` (gerado em 277.pplx.app / biomni), com foco em:

1. Detecção do motivo conservado de Mönttinen (36 resíduos) via TM-align local.
2. Presença/ausência do motivo na estrutura 3D dos 6 folds query.
3. Correção crítica: o relatório original usou cadeia errada para RNAP (1HQM-B = α, sem palm), substituída por 1HQM-C (β, com domínio palm PF00562).

---

## Resultados principais

| Métrica | Valor |
|---|---|
| Estruturas comparadas | 7 (1 referência + 6 queries) |
| Core ≥ 5/6 folds | **48 resíduos** (Mönttinen 2016: 36, 89 estruturas) |
| Core universal (6/6) | **4 resíduos** (Phe34, Val35, Lys197, Ala199) |
| Folds PRESENTE (≥70%) | 5/6 (PolX é PARCIAL, 66.7%) |
| ASP catalíticos no core | **Asp22 (5/6), Asp33 (5/6)** — evidência do mecanismo two-metal-ion |

### Cobertura por fold

| Fold | PDB | Cα | Cobertura | Status |
|---|---|---|---|---|
| NTase (ref) | 1VFG-A | 342 | 48/48 (100%) | Referência |
| PolA | 3BDP-A | 580 | 46/48 (95.8%) | PRESENTE |
| RT   | 1RTD-A | 554 | 45/48 (93.8%) | PRESENTE |
| PolB | 1WAJ-A | 903 | 44/48 (91.7%) | PRESENTE |
| RNAP | **1HQM-C** | 997 | 40/48 (83.3%) | PRESENTE (corrigido) |
| RdRp | 1RDR-A | 316 | 37/48 (77.1%) | PRESENTE |
| PolX | 1BPX-A | 331 | 32/48 (66.7%) | PARCIAL |

---

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `monttinen_core.py` | Script principal: download PDB, TM-align, parsing, detecção de core, geração de outputs |
| `analise_tmalign_monttinen_corrigida.md` | Relatório completo com contestação quantitativa ao relatório original |
| `index.html` | Website interativo com KPIs, heatmap 7×48, strip de sequência, gráfico TM-score (servido via GitHub Pages) |
| `tmscore_summary.tsv` | TM-scores e RMSD dos 6 pares NTase↔fold |
| `core_residues.tsv` | 48 resíduos do core ≥5/6, anotados por fold |
| `presence_absence.tsv` | Presença/ausência por fold (% cobertura) |
| `color_core.pml` | Script PyMOL: carrega 7 estruturas e colore o core |
| `color_<fold>.pml` | Script PyMOL individual por fold |

---

## Dependências

- TM-align v20220412
- Python 3.12
- BioPython 1.87
- pandas

---

## Como reproduzir

```bash
# Download das 7 estruturas do RCSB PDB
for id in 1VFG 3BDP 1WAJ 1RTD 1BPX 1HQM 1RDR; do
  wget "https://files.rcsb.org/download/${id}.pdb"
done

# Roda a análise
python3 monttinen_core.py
```

Outputs gerados em `~/monttinen_core/`.

---

## Visualização

- **Website:** https://monteirotorres.github.io/277-test/ (GitHub Pages) — ou abrir `index.html` localmente em qualquer browser moderno.
- **PyMOL:** `File → Run Script → color_core.pml`.

---

## Erro crítico corrigido vs relatório original

| Métrica | Relatório original (pplx/biomni) | Esta análise |
|---|---|---|
| Cadeia RNAP (1HQM) | B (α, 229 Cα, sem palm) | **C (β, 997 Cα, com palm)** |
| TM NTase ↔ RNAP | 0.242 | **0.315 (+30%)** |
| Core universal (6/6) | "1 resíduo — artefato" | **4 resíduos** |
| Core ≥ 5/6 | 32 resíduos | **48 resíduos** |
| ASP catalíticos | Não reportados | **Asp22, Asp33** |

---

## Limitações

1. Apenas 1 representante por fold (Mönttinen 2016 usou 89 estruturas em 5 superfamílias).
2. Sem controle negativo (folds não-palm: Rossmann, TIM-barrel).
3. Single-rep produz core mais frouxo (48 vs 36) — convergiria com múltiplos representantes.

---

## Referências

- Mönttinen HAM et al. *MBE* 33:1697 (2016) — 36-residue structural core, HSF tool.
- Zhang Y & Skolnick J. *NAR* 33:2302 (2004) — TM-align/TM-score.
