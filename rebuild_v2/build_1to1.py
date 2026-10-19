#!/usr/bin/env python3
"""
build_1to1.py -- 1:1 FGF2-FGFR systems for all seven orthologs, both receptors.

Why these are derived from the existing 2:2:2 and 2:2:1 builds rather than from the 1EVT and
1DJS crystals: the FGFR chains of 1EVT (chain C) and 1DJS (chain A) have 10-12 residue gaps
(1EVT 294-303 and 312-323; 1DJS 297-306, plus shorter ones) that would need loop modelling,
and a separately modelled receptor would no longer be comparable with Table 1. The rebuilt
2:2:2 and 2:2:1 systems already carry a complete, CHARMM-GUI-built FGFR with those loops
modelled, and an FGF2 core trimmed and placed identically across species. Deleting the second
FGF2, the second FGFR and all heparin from those builds leaves exactly the 1:1 complex in its
crystallographic interface geometry, built by the same protocol.

  1fq9 -> PROA (FGFR1) + PROC (FGF2)      "1to1_fgfr1"
  1e0o -> PROA (FGFR2) + PROC (FGF2)      "1to1_fgfr2"

Per system: solute, re-solvated in a box sized to this smaller complex (12 A padding),
neutralised and brought to 0.15 M KCl, then N_REPS replicate folders with NAMD inputs and
slurm scripts. No heparin, so the carbohydrate dihedral restraints are dropped.

Usage:  python3 build_1to1.py [tmpl/species ...]      (default: all 14)
"""
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUT = BASE / "onetoone"
VMD_DIR = "/Applications/VMD 2.0.0a7-pre2.app/Contents/vmd2/lib"
NAMD_BIN = "/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2"
SRC = {"fgfr1": "1fq9", "fgfr2": "1e0o"}        # receptor label -> source template
LABEL = {"fgfr1": "11r1", "fgfr2": "11r2"}
SPECIES = ["human", "human_stab", "bovine_stab", "chinook_salmon", "goldfish",
           "japanese_ricefish", "zebrafish"]
N_REPS = 5
PAD = 12.0                      # A of water between solute and box edge
BACKBONE = {"N", "CA", "C", "O"}

BUILD_TCL = r"""
package require psfgen
package require solvate
package require autoionize
resetpsf
readpsf {solute_psf}
coordpdb {solute_pdb}
# keep one FGFR (PROA) and its paired FGF2 (PROC); drop the partner pair and all heparin
delatom PROB
delatom PROD
{del_carb}
writepsf pair.psf
writepdb pair.pdb
set m [mol new pair.psf]
mol addfile pair.pdb waitfor all
set all [atomselect $m all]
$all moveby [vecscale -1.0 [measure center $all]]
$all writepdb pair_centered.pdb
puts "PAIR_MINMAX [measure minmax $all]"
mol delete $m
set h [expr {box}/2.0]
solvate pair.psf pair_centered.pdb -minmax [list [list -$h -$h -$h] [list $h $h $h]] -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -cation POT -anion CLA -seg IONS -o ionized
quit
"""


def run_vmd(tcl, cwd):
    env = dict(os.environ, VMDDIR=VMD_DIR)
    p = subprocess.run([f"{VMD_DIR}/vmd_MACOSXARM64", "-dispdev", "text", "-e", tcl],
                       cwd=cwd, env=env, capture_output=True, text=True)
    (cwd / "vmd_build.log").write_text(p.stdout + p.stderr)
    return p.stdout


def psf_atoms(psf):
    out = []
    with open(psf) as f:
        for line in f:
            if "!NATOM" in line:
                n = int(line.split()[0])
                break
        for _ in range(n):
            p = next(f).split()
            out.append((p[1], p[2], p[3], p[4]))
    return out


def write_posres(pdb_in, psf, out):
    """CHARMM-GUI convention: protein backbone 1.0, other protein heavy atoms 0.5, else 0."""
    atoms = psf_atoms(psf)
    lines = [l for l in open(pdb_in) if l.startswith(("ATOM", "HETATM"))]
    assert len(lines) == len(atoms), f"PDB/PSF mismatch {len(lines)} vs {len(atoms)}"
    with open(out, "w") as fh:
        for l, (seg, _, _, name) in zip(lines, atoms):
            k = 0.0
            if seg.startswith("PRO") and not name.startswith("H"):
                k = 1.0 if name in BACKBONE else 0.5
            fh.write(f"{l[:54]}  1.00{k:6.2f}{l[66:]}")
        fh.write("END\n")


def strip_carb_restraints(txt):
    """Remove the heparin chair-restraint block; these systems have no carbohydrate."""
    txt = re.sub(r"(?m)^# carbohydrate chair restraint\n", "", txt)
    txt = re.sub(r"(?m)^exec sed .*carbohydrate_restraint\.dat.*\n", "", txt)
    txt = re.sub(r"(?m)^exec sed .*carbohydrate\\_restraint\.dat.*\n", "", txt)
    txt = re.sub(r"(?m)^extraBondsFile .*\n", "", txt)
    txt = re.sub(r"(?m)^extraBonds\s+on.*\n", "", txt)
    return txt


def slurm(kind, job, hours):
    step = "step4_equilibration" if kind == "eq" else "step5_production"
    log, tag = ("eq.log", "e") if kind == "eq" else ("prod.log", "p")
    return f"""#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu
#SBATCH --job-name={job}
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time={hours:02d}:00:00
#SBATCH --output=slurm{tag}-%j.out

{NAMD_BIN} +p${{SLURM_CPUS_PER_TASK:-4}} {step}.inp > {log}
"""


def donor_inputs(tk, sp):
    """A replicate folder holding step4/step5 inputs to copy settings from.

    Some uploads came back without namdrep1 or without the .inp files, so fall back to any
    replicate of this species, then to any species of the same template. The box size and the
    seed are rewritten per system, so any donor of the right template is equivalent.
    """
    need = ("step4_equilibration.inp", "step5_production.inp")
    for root in [BASE / tk / sp] + sorted((BASE / tk).glob("*")):
        for rep in sorted((root / "sim").glob("namdrep[0-9]")) if (root / "sim").is_dir() else []:
            if all((rep / f).exists() for f in need):
                return rep
    raise SystemExit(f"no NAMD input donor found for {tk}")


def build(rec, sp):
    tk = SRC[rec]
    src = BASE / tk / sp / "build"
    if not (src / "solute.psf").exists():
        raise SystemExit(f"missing {src}/solute.psf -- run the main rebuild for {tk}/{sp} first")
    sysdir = OUT / rec / sp
    bdir, sdir = sysdir / "build", sysdir / "sim"
    for d in (bdir, sdir):
        shutil.rmtree(d, ignore_errors=True)
    bdir.mkdir(parents=True)

    carb = sorted({a[0] for a in psf_atoms(src / "solute.psf") if a[0].startswith("CAR")})
    del_carb = "\n".join(f"delatom {c}" for c in carb)

    # first pass with a placeholder box, to measure the pair; then rebuild at the right size
    probe = BUILD_TCL.format(solute_psf=src / "solute.psf", solute_pdb=src / "solute.pdb",
                             del_carb=del_carb, box=200.0)
    (bdir / "probe.tcl").write_text(probe)
    out = run_vmd("probe.tcl", bdir)
    m = re.search(r"PAIR_MINMAX \{(.*?)\} \{(.*?)\}", out)
    if not m:
        raise SystemExit(f"could not measure the pair in {sysdir}; see vmd_build.log")
    lo, hi = [list(map(float, g.split())) for g in m.groups()]
    box = max(h - l for l, h in zip(lo, hi)) + 2 * PAD
    box = float(int(box) + 1)

    for f in ("solvated.psf", "solvated.pdb", "ionized.psf", "ionized.pdb"):
        (bdir / f).unlink(missing_ok=True)
    tcl = BUILD_TCL.format(solute_psf=src / "solute.psf", solute_pdb=src / "solute.pdb",
                           del_carb=del_carb, box=box)
    (bdir / "build.tcl").write_text(tcl)
    run_vmd("build.tcl", bdir)
    if not (bdir / "ionized.psf").exists():
        raise SystemExit(f"VMD build failed in {bdir}; see vmd_build.log")

    sdir.mkdir()
    scaf_rep = donor_inputs(tk, sp)
    def pick(name):
        own = BASE / tk / sp / "sim" / name
        if own.exists():
            return own
        for other in sorted((BASE / tk).glob("*")):
            if (other / "sim" / name).exists():
                return other / "sim" / name
        raise SystemExit(f"no {name} found for {tk}")

    shutil.copytree(pick("toppar"), sdir / "toppar")
    tc = BASE / tk / sp / "sim" / "toppar_charmm"
    shutil.copytree(tc if tc.exists() else pick("toppar_charmm"), sdir / "toppar_charmm")
    pbc = pick("step3_pbcsetup.str").read_text()
    for k in ("A", "B", "C"):
        pbc = re.sub(rf"SET {k} = \S+", f"SET {k} = {box:g}", pbc)
    (sdir / "step3_pbcsetup.str").write_text(pbc)

    for r in range(1, N_REPS + 1):
        rd = sdir / f"namdrep{r}"
        (rd / "restraints").mkdir(parents=True)
        shutil.copy(bdir / "ionized.psf", rd / "step3_input.psf")
        shutil.copy(bdir / "ionized.pdb", rd / "step3_input.pdb")
        write_posres(bdir / "ionized.pdb", bdir / "ionized.psf", rd / "restraints" / "prot_posres.ref")
        # distinct per receptor, species and replicate, so no two runs share a velocity seed
        seed = 30000 + 10000 * (0 if rec == "fgfr1" else 1) + 100 * SPECIES.index(sp) + 10 * r
        for inp in ("step4_equilibration.inp", "step5_production.inp"):
            txt = (scaf_rep / inp).read_text()
            txt = strip_carb_restraints(txt)
            txt = re.sub(r"(?m)^parameters( +)(\.\./)?toppar/", r"parameters\1../toppar/", txt)
            # the source input already carries a seed from the main rebuild; replace it
            # rather than adding a second line
            txt, n = re.subn(r"(?m)^seed\s+\d+;", f"seed                    {seed + (inp[4] == '5')};", txt)
            if n == 0:
                txt = txt.replace("set temp                310;",
                                  f"set temp                310;\nseed                    {seed + (inp[4] == '5')};", 1)
                n = 1
            assert n == 1, f"expected one seed line, found {n}"
            (rd / inp).write_text(txt)
        job = f"{LABEL[rec]}_{sp}_r{r}"
        (rd / "eq_slurm.sh").write_text(slurm("eq", job + "_eq", 4))
        (rd / "prod_slurm.sh").write_text(slurm("prod", job + "_pr", 8))
        (rd / "submit.sh").write_text(
            '#!/bin/bash\n# equilibration, then production only if equilibration finished cleanly\n'
            'cd "$(dirname "$0")"\njid=$(sbatch --parsable eq_slurm.sh)\n'
            'sbatch --dependency=afterok:$jid prod_slurm.sh\n')
        os.chmod(rd / "submit.sh", 0o755)

    atoms = psf_atoms(bdir / "ionized.psf")
    ions = {k: sum(1 for a in atoms if a[2] == k) for k in ("POT", "CLA")}
    segs = sorted({a[0] for a in atoms if a[0].startswith("PRO")})
    summary = (f"{rec}/{sp}: {len(atoms)} atoms, box {box:g} A, segments {','.join(segs)}, "
               f"POT {ions['POT']} CLA {ions['CLA']}")
    (sysdir / "build_summary.txt").write_text(summary + "\n")
    return summary


def main():
    targets = sys.argv[1:] or [f"{r}/{s}" for r in SRC for s in SPECIES]
    for t in targets:
        rec, sp = t.split("/")
        print(build(rec, sp), flush=True)
    sub = OUT / "submit_1to1.sh"
    sub.write_text(
        '#!/bin/bash\n'
        '# Run ON THE CLUSTER from onetoone/. Submits equilibration for every replicate, with\n'
        '# production chained on --dependency=afterok.\n'
        '#   bash submit_1to1.sh            # both receptors, 14 systems x 5 replicates\n'
        '#   bash submit_1to1.sh fgfr1      # one receptor only\n'
        'cd "$(dirname "$0")"\n'
        'targets=("$@"); [ ${#targets[@]} -eq 0 ] && targets=(fgfr1 fgfr2)\n'
        'for t in "${targets[@]}"; do\n'
        '  for rep in "$t"/*/sim/namdrep*; do\n'
        '    [ -f "$rep/submit.sh" ] || continue\n'
        '    echo "== $rep"; bash "$rep/submit.sh"\n'
        '  done\n'
        'done\n')
    sub.chmod(0o755)
    print("wrote", sub)


if __name__ == "__main__":
    main()
