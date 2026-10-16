# Analysis plan: FGF2 86–97 segment-swap chimeras (fixed before any chimera MD was run)

Written 2026-10-16, after the chimera systems were built (`make_chimeras.py`) and before any of
them were simulated. It is committed to git before the runs start. Results are reported whichever
way they come out. Runs will not be lengthened, and replicates will not be added or dropped,
because of the results.

## Question
In the 2:2:1 assembly (1E0O, one shared heparin), the conserved Arg90 (FGF2 155-aa numbering)
binds the shared heparin in goldfish FGF2 but not in human FGF2. This finding was exploratory,
and it was replicated in goldfish rep 3 against a recorded prediction
(`ext_results/PREDICTION_goldfish_rep3.md`). Eight of the 25 human→goldfish substitutions cluster
around Arg90, in FGF2 86–97: K86N E87S L92I A93G S94T K95R C96R V97T.
**Does that segment control whether Arg90 engages heparin?**

## Systems
- `chim_human_gfseg`: human core with goldfish residues at FGF2 86–97.
- `chim_goldfish_hseg`: goldfish core with human residues at FGF2 86–97.
- Each in 2:2:1 (primary) and 2:2:2 (control); 3 replicates each.
- Protocol identical to the parent extensions: equilibration, 1 ns production, then 50 × 1 ns
  segments. Seeds are distinct from the parents'.
- A replicate enters the analysis only when all 50 segments are complete.

## Measurement (identical to the parent analysis)
- Window: extension segments 11–50.
- **Primary.** Arg90–heparin heavy-atom contacts (4.5 Å) for each FGF2 copy, the window mean every
  0.5 ns. Computed with `analyze_r90.sh`, which reproduces the parent values exactly.
- **Confirmatory.** Arg90 per-residue contribution from the whole-assembly decomposition
  (idecomp = 1, every 1 ns; `analyze_exploratory.sh`).

Parent reference (`ext_results/r90_parent_reference.txt`, same script), contacts for copy C / copy D:

| | 2:2:1 | 2:2:2 |
|---|---|---|
| goldfish | 11.1/9.6, 36.1/16.1, 18.9/15.8 | 3.1/0, 4.1/0, 2.6/0 |
| human | 2.7/0, 3.4/0.1, 5.5/0 | 4.6/0, 2.5/0, 3.3/0 |

All goldfish 2:2:1 copies have ≥ 9.6 contacts; every other copy has ≤ 5.5.

## Definitions (fixed now)
- **Engaged copy:** window-mean Arg90–heparin contacts > 7.5 (the midpoint of the parent gap
  5.5–9.6) **and** an Arg90 per-residue contribution below −5 kcal/mol.
- **Engaged replicate:** at least one engaged copy.

## Hypotheses and decision rules (2:2:1)
- **H1, the segment is sufficient:** `chim_human_gfseg` has ≥ 2 of 3 replicates engaged.
- **H2, the segment is necessary:** `chim_goldfish_hseg` has ≤ 1 of 3 replicates engaged.

| H1 | H2 | Conclusion |
|---|---|---|
| yes | yes | The 86–97 segment controls the Arg90 binding mode. Recommend the chimera and its reverse as experimental constructs. |
| yes | no | Sufficient but not necessary: other goldfish positions can also produce the mode. |
| no | yes | Necessary but not sufficient: the segment acts together with other goldfish positions. |
| no | no | The segment does not control the mode; look elsewhere (for example the L135A heparin-contact substitution). |

## Control (2:2:2)
- Prediction: neither chimera has more than 1 of 3 replicates engaged, as for both parents.
- An engaged 2:2:2 replicate would mean that the mode is no longer specific to the shared heparin.

## Secondary (reported descriptively, no decision attached)
- Whole-assembly ΔG, FGF2–heparin and FGF2–FGFR (heparin removed) energies for each chimera
  against both parents (mean ± SD, Welch p).
- Interpretation caveat: the human→goldfish segment adds one net positive charge per FGF2 (C96R,
  E87S, K86N). A higher Arg90 engagement in `chim_human_gfseg` could partly reflect a general
  electrostatic attraction rather than repositioning. The FGF2–heparin energy and the
  canonical-site salt bridges are reported to help separate the two.
- With n = 3 per system the power is low. Decisions use replicate counts against a threshold set
  from the parents, not p-values.

## Pipeline after the runs finish
1. `JOBS=4 bash run_all_mmpbsa_v2.sh rebuild_v2/1e0o/chim_human_gfseg rebuild_v2/1e0o/chim_goldfish_hseg rebuild_v2/1fq9/chim_human_gfseg rebuild_v2/1fq9/chim_goldfish_hseg`
   (1 ns MM/GBSA; provides complex_chamber.parm7 and the splits)
2. `bash analyze_extension.sh <rep>`, then `bash analyze_exploratory.sh <rep>`, then
   `bash analyze_r90.sh <rep>` for each replicate.
3. Tabulate engaged copies and replicates by the definitions above.
