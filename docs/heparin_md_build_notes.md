# Heparin MD Parameterization Notes
# For use with NAMD + CHARMM36 (par_all36_carb.prm)
# Context: FGF–FGFR–heparin quaternary assembly / MM-GBSA pipeline

---

## 1. Force Field Stack

Load all three parameter files in your NAMD .conf:

```tcl
parameters    par_all36_carb.prm
parameters    par_all36_cgenff.prm
parameters    par_water_ions.prm
```

Also load any `.str` stream files from CHARMM-GUI's `toppar/` folder — these contain
sulfated residue parameters that supplement the main .prm files.

---

## 2. Heparin Composition

Major repeating disaccharide:
  →4)-α-L-IdoA(2S)-(1→4)-α-D-GlcNS(6S)-(1→

Per disaccharide unit: ~3-4 negative charges from sulfate/carboxylate groups.
An 8-disaccharide (16-residue) chain carries roughly -24e total.

Residue names in CHARMM36:
  - α-L-IdoA(2S)  →  AIDOA  +  2-O-sulfate modification
  - α-D-GlcNS(6S) →  BGLCN  +  N-sulfate (Sulfoamino at position 2)
                              +  6-O-sulfate (Sulfate at position 6)

Alternative names for IdoA2S (useful for literature/PDB searches):
  - 2-O-sulfo-α-L-iduronic acid
  - L-iduronic acid 2-sulfate
  - IduA2S
  - UA2S (Shively & Conrad notation)
  - I2S  (CHARMM topology code)
  - IDS  (PDB CCD code)

---

## 3. Correct GRS (Glycan Residue Sequence) for CHARMM-GUI

For a 3-disaccharide fragment:

```
1              AIDOA_2SUF
2 - 14B:       BGLCN_2NSF_6SUF
3 - - 14A:     AIDOA_2SUF
4 - - - 14B:   BGLCN_2NSF_6SUF
5 - - - - 14A: AIDOA_2SUF
6 - - - - - 14B: BGLCN_2NSF_6SUF
```

Linkage pattern alternates strictly: 14B (IdoA→GlcN) / 14A (GlcN→IdoA)
The B/A prefix on BGLCN/AGLCN reflects the anomeric config of the glycosidic bond
to the NEXT residue — CHARMM-GUI assigns this automatically when you set linkage type.

Modification selections in CHARMM-GUI Glycan Modeler:
  - AIDOA: Sulfate at position 2
  - BGLCN: Sulfoamino at position 2  (N-sulfate — NOT "Sulfate", which is O-linked)
            Sulfate at position 6     (6-O-sulfate)

---

## 4. Known CHARMM-GUI Bug: 2BNSF Patch Failure

The Glycan Modeler currently fails when applying N-sulfate (Sulfoamino) to GlcN.
Error in step1_glycan_modeler.out:

```
patch 2BNSF CARB 2 setup warn
  Atom HO2 in PATCH 2BNSF was to be deleted in residue AGLCN/BGLCN 2 but is not present.
  Atom O2  in PATCH 2BNSF was to be deleted in residue AGLCN/BGLCN 2 but is not present.
Warning from PSFSUM: The sum of charges (-2.100000) is not an integer
ABNORMAL TERMINATION
```

Root cause: patch 2BNSF is written to delete O2/HO2 (a hydroxyl at C2), but GlcN has
an amine (NH2) at C2, not a hydroxyl. The GUI is selecting the wrong patch variant.
This occurs regardless of anomeric config (AGLCN or BGLCN) or linkage type.

Status: Confirmed bug, not a user error. All GRS combinations attempted:
  - AIDOA_2SUF / BGLCN_2NSF_6SUF with 14B linkage  → fails
  - AIDOA_2SUF / AGLCN_2NSF_6SUF with 14A linkage  → fails
  - Single disaccharide                              → fails

---

## 5. Workarounds (in order of preference)

### Option A: Pre-built heparin PDB from RCSB
Search RCSB for heparin oligomers, then use CHARMM-GUI PDB Reader & Manipulator
to convert to CHARMM36 PSF/PDB with NAMD output.
Useful PDB entries: 1HPN, 3IRJ (heparin-containing crystal structures).

### Option B: Manual psfgen
Requires: top_all36_carb.rtf + top_all36_cgenff.rtf

```tcl
package require psfgen
topology top_all36_carb.rtf
topology top_all36_cgenff.rtf

segment HEPN {
    residue 1 IDOA
    residue 2 AGLCNA
    residue 3 IDOA
    residue 4 AGLCNA
    # ... repeat as needed
}

# Glycosidic linkage patches (1->4)
patch 14bb HEPN:1 HEPN:2
patch 14aa HEPN:2 HEPN:3
# ... alternate

# Sulfation patches
patch SI2  HEPN:1   ;# 2-O-sulfate on IdoA
patch NSO3 HEPN:2   ;# N-sulfate on GlcN
patch SO6  HEPN:2   ;# 6-O-sulfate on GlcN
# ... repeat for each residue

coordpdb heparin.pdb HEPN
guesscoord
writepsf heparin.psf
writepdb heparin.pdb
```

NOTE: Patch names (SI2, NSO3, SO6) must be verified against your actual
top_all36_carb.rtf — names vary between CHARMM36 versions.

### Option C: Post on CHARMM-GUI forum
Submit the step1_glycan_modeler.out error + GRS string.
Developers are responsive to patching bugs.
URL: charmm-gui.org/forum

---

## 6. Solvation and Ionization (VMD)

```tcl
package require solvate
package require autoionize

# 12-15 Å padding — heparin is extended
solvate heparin.psf heparin.pdb -t 15 -o solvated

# 150 mM NaCl, neutralize
autoionize -psf solvated.psf -pdb solvated.pdb \
           -sc 0.15 -cation SOD -anion CLA \
           -o ionized
```

Use Joung-Cheatham ion parameters for Na+/Cl- — validated for polyanionic systems.
A 3-disaccharide heparin carries ~-9e; expect significant Na+ addition.

---

## 7. NAMD Configuration Notes

```tcl
# Electrostatics
PME                  yes
PMEGridSpacing       1.0

# Cutoffs
cutoff               12.0
switchdist           10.0
pairlistdist         14.0

# Timestep
timestep             1.0    ;# fs; conservative for initial runs

# Temperature (physiological)
langevin             on
langevinDamping      1.0
langevinTemp         310

# Pressure (NPT)
LangevinPiston       on
LangevinPistonTarget 1.01325
LangevinPistonPeriod 200
LangevinPistonDecay  100
LangevinPistonTemp   310
```

PME grid size must be factorable by 2, 3, 5 — set PMEGridSizeX/Y/Z manually if NAMD warns.

---

## 8. Staged Minimization Protocol

Heparin simulations often crash early due to bad sulfate geometry from guessed coords.

```tcl
# Stage 1: restrained minimization
constraints          on
consref              heparin.pdb
conskfile            heparin.pdb    ;# B-factor column = force constant
conskcol             B
constraintScaling    10.0
minimize             2000

# Stage 2: reduced restraints
constraintScaling    1.0
minimize             2000

# Stage 3: free minimization
constraints          off
minimize             5000
```

---

## 9. Post-Equilibration Validation Checklist

| Target                        | Method                                      |
|-------------------------------|---------------------------------------------|
| Glycosidic angles (phi/psi)   | Compare to NMR/crystal data                 |
| Sulfate geometry              | S-O bond ~1.46 Å, S-O-C angles             |
| IdoA ring pucker              | Cremer-Pople parameters (4C1 vs 2SO)        |
| Ion binding                   | RDF of Na+ around sulfate O (peak ~2.3 Å)  |
| RMSD stability                | Should stabilize within 10-50 ns            |

IdoA ring flexibility note: IdoA populates both 4C1 and 2SO (skew-boat) conformations.
Consider H-REMD or metadynamics to adequately sample ring pucker in production.

---

## 10. Common Pitfalls

- Missing parameter errors: check NAMD log for "ERROR: MISSING" before production
- PME grid size: adjust PMEGridSizeX/Y/Z if NAMD warns
- Patch ordering in psfgen: linkage patches before sulfation patches
- Box too small: use >= 15 Å buffer for extended heparin chains
- Under-neutralization: causes artifacts with PME electrostatics
- Don't use generic GAFF parameters for sulfates — poorly validated
- IdoA ring flexibility: do not ignore in production runs
