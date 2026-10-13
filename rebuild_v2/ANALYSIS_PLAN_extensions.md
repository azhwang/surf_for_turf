# Analysis plan: 50 ns extensions (fixed before any extension data was analysed)

Written 2026-10-13, after the 2:2:2 trajectories were copied but before any of them
were analysed, and committed to git before the analysis ran. Results are reported
whichever way they come out. Runs will not be lengthened, and replicates will not be
added or dropped, because of the results.

## Data
- 2:2:2 (1FQ9): human reps 1, 2, 3 and goldfish reps 2, 3, 5.
- 2:2:1 (1E0O): human reps 1, 2, 4 and goldfish reps 1, 2, 3.
- Each run is 1 ns of production (`step5`) followed by 50 × 1 ns extension segments
  (`ext01`–`ext50`), with frames every 100 ps.
- A replicate enters the primary analysis only when all 50 segments are complete.
  Incomplete replicates (2:2:2 human rep 1, 2:2:1 goldfish rep 3) are added
  unchanged when they finish. Until then, the n per group is stated.

## Window
- Primary analysis window: extension segments 11–50 (40 ns). Segments 1–10 are
  treated as further equilibration. Full time series are still shown.

## Metrics
1. **Heparin-site salt bridges (primary structural metric).** For each FGF copy, the
   number of heparin-site Lys/Arg side-chain N atoms (NZ, NE, NH1, NH2) within 4.0 Å of
   heparin sulfate or carboxylate O atoms (OS22–OS24, OS62–OS64, O61, O62), every frame.
   - Site: the five canonical positions FGF2 K128, R129, K134, K138, K144 (core
     positions 104, 105, 110, 114, 120 in 1FQ9; 101, 102, 107, 111, 117 in 1E0O).
     The secondary variant counts all FGF Lys/Arg.
   - Heparin: each FGF copy's own interface heparin in 2:2:2 (FGF C with CARA, FGF D
     with CARB); the single shared heparin in 2:2:1.
2. **MM/GBSA binding energy (primary energetic metric).** Same settings as the 1 ns
   analysis (GB igb = 5, mbondi2, 0.15 M, FGF = ligand, FGFR + heparin = receptor,
   interface 1, interface 2 and whole assembly), with no decomposition.
   - Frames every 0.5 ns over the window: 80 frames per replicate.
   - Within-replicate uncertainty from 4 blocks of 10 ns.
3. **FGFR–FGFR distance (secondary; tests the v1 "smaller receptor distance"
   claim).** The distance between the Cα centres of mass of the two FGFR chains,
   every frame.

## 2:2:1 heparin state (defined in advance)
- Interface 1, FGF C: a frame is "tight" if it has ≥ 14 canonical-site salt bridges.
  14 is the midpoint between the 1 ns group means (about 16.5 tight, about 11 loose).
- Reported per replicate: the fraction of window frames that are tight, and how often
  replicates switch state.

## Statistics
- The unit is the replicate mean over the window (n = 3 per group when complete).
- Goldfish vs human, separately per stoichiometry and per metric: the difference in
  means with a 95% CI and a Welch t-test.
- **Test of the v1 manuscript claim, specified exactly:** "goldfish FGF2 has a lower
  (more favourable) binding energy than human FGF2 in the 2:2:1 assembly". This counts
  as supported only if the whole-assembly 2:2:1 ΔG(goldfish) − ΔG(human) is below 0
  with Welch p < 0.05. Interface-level and structural metrics are reported as
  secondary and do not substitute for it.
- With n = 3 the power is low. A non-significant difference is reported as "not
  resolved", not as "no difference".
