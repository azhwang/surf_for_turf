# Prediction for 2:2:1 goldfish rep 3 (recorded 2026-10-15, before its data were analysed)

Exploratory finding in goldfish reps 1-2 vs human reps 1, 2, 4 (2:2:1, 11-50 ns): conserved FGF R90
(FGF2 155-aa numbering; core position 63; complex residues 466 = FGF C, 593 = FGF D) engages the
shared heparin far more in goldfish (heavy-atom contacts within 4.5 A, window mean:
FGF C 11.1 / 36.1, FGF D 9.6 / 16.1) than in human (FGF C 2.7-5.5, FGF D 0.0-0.1), and its
per-residue MM/GBSA contribution is about -11.5 vs -0.5 to -2 kcal/mol.

Prediction: in goldfish rep 3, the window-mean R90-heparin contacts exceed the human maximum
(FGF C > 5.5 and FGF D > 0.1), and R90's per-residue contribution is below -5 kcal/mol in at least
one FGF copy. If not, the R90 pattern is treated as not replicated.

## Outcome (scored 2026-10-16, after rep 3 completed all 50 segments)
- R90-heparin contacts, window mean: FGF C 18.9 (> 5.5), FGF D 15.8 (> 0.1). Met.
- R90 per-residue contribution: FGF C -10.6, FGF D -9.4 kcal/mol (< -5). Met.
- **Prediction confirmed.** The R90-heparin engagement in 2:2:1 replicates in all 3 goldfish replicates
  and both FGF copies; in human it is absent or weak in all 3. In 2:2:2 R90 contributes the same in
  both species (FGF C -2.0 vs -2.7, FGF D -0.1 vs -0.1), so the difference is specific to the shared-heparin assembly.
- Caveat: the net FGF-heparin energy barely differs (2:2:1 goldfish - human -6.1, p = 0.77) because
  canonical-site residues (K128, R129, K134, K35) lose some of their contribution in goldfish. This is a
  difference in heparin-binding mode, not a demonstrated difference in net affinity.
