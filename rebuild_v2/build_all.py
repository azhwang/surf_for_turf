#!/usr/bin/env python3
"""
build_all.py -- step 3 of the consistent FGF rebuild: build, solvate, ionize,
and write ready-to-submit NAMD replicate folders for every system.

Per system (<tmpl>/<species>/):
  build/   psfgen: FGFR + heparin read from the scaffold PSF (keeps CHARMM-GUI's
           heparin sulfation/linkage patches), FGF PROC/PROD rebuilt from the
           trimmed, placed cores (NTER/CTER, missing atoms by guesscoord);
           centred at the origin, solvated in the SAME cubic box as the original
           builds (1FQ9 124 A, 1E0O 129 A), neutralised + 0.15 M KCl (POT/CLA).
  sim/     step3_pbcsetup.str, toppar/ (from the scaffold build), and
           namdrep1..5/ with step3_input.psf/pdb, restraints (regenerated with
           CHARMM-GUI's conventions), step4/step5 .inp (the scaffold's own files,
           plus an explicit per-replicate `seed`), eq/prod slurm scripts and
           submit.sh, which chains prod after eq with --dependency=afterok.

Usage:  python3 build_all.py [tmpl/species ...]     (default: all 14)
"""
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
MD = BASE.parent
VMD_DIR = "/Applications/VMD 2.0.0a7-pre2.app/Contents/vmd2/lib"
SCAFFOLD = {"1fq9": MD / "all_metrics_analysis/1fq9_goldfish", "1e0o": MD / "all_metrics_analysis_221/1e0o_chinook_salmon"}
BOX = {"1fq9": 124.0, "1e0o": 129.0}
LABEL = {"1fq9": "222", "1e0o": "221"}
N_REPS = 5
BACKBONE = {"N", "CA", "C", "O"}
SUGAR_RING = {"C1", "C2", "C3", "C4", "C5", "O5"}
NAMD_BIN = "/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2"

BUILD_TCL = r"""
package require psfgen
package require solvate
package require autoionize
resetpsf
topology {toppar}/top_all36_prot.rtf
topology {toppar}/top_all36_carb.rtf
readpsf {scaf_psf}
coordpdb {scaf_pdb}
delatom PROC
delatom PROD
pdbalias residue HIS HSD
foreach seg {{PROC PROD}} {{
    segment $seg {{ first NTER; last CTER; pdb fgf_$seg.pdb }}
    coordpdb fgf_$seg.pdb $seg
}}
guesscoord
writepsf solute.psf
writepdb solute.pdb
set m [mol new solute.psf]
mol addfile solute.pdb waitfor all
set all [atomselect $m all]
$all moveby [vecscale -1.0 [measure center $all]]
set mm [measure minmax $all]
puts "SOLUTE_MINMAX $mm"
$all writepdb solute_centered.pdb
mol delete $m
set h [expr {box}/2.0]
solvate solute.psf solute_centered.pdb -minmax [list [list -$h -$h -$h] [list $h $h $h]] -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -cation POT -anion CLA -seg IONS -o ionized
quit
"""


def run_vmd(tcl, cwd):
    env = dict(os.environ, VMDDIR=VMD_DIR)
    p = subprocess.run([f"{VMD_DIR}/vmd_MACOSXARM64", "-dispdev", "text", "-e", tcl],
                       cwd=cwd, env=env, capture_output=True, text=True)
    (cwd / "vmd_build.log").write_text(p.stdout + p.stderr)
    if "ionized.psf" not in os.listdir(cwd):
        raise RuntimeError(f"VMD build failed in {cwd}; see vmd_build.log")
    return p.stdout


def psf_atoms(psf):
    """[(seg, resid, resname, atomname)] in PSF order."""
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
    """CHARMM-GUI convention (column B): protein backbone N/CA/C/O 1.0, other protein
    heavy atoms 0.5; carbohydrate ring C1-C5/O5 1.0, other carbohydrate heavy 0.5; else 0."""
    atoms = psf_atoms(psf)
    lines = [l for l in open(pdb_in) if l.startswith(("ATOM", "HETATM"))]
    assert len(lines) == len(atoms), "PDB/PSF atom count mismatch"
    with open(out, "w") as fh:
        for l, (seg, _, _, name) in zip(lines, atoms):
            heavy = not name.startswith("H")
            if seg.startswith("PRO") and heavy:
                k = 1.0 if name in BACKBONE else 0.5
            elif seg.startswith("CAR") and heavy:
                k = 1.0 if name in SUGAR_RING else 0.5
            else:
                k = 0.0
            fh.write(f"{l[:54]}  1.00{k:6.2f}{l[66:]}")
        fh.write("END\n")


def remap_carb_restraints(old_dat, old_psf, new_psf, out):
    """Carbohydrate chair restraints: translate atom indices via (seg, resid, atom) identity.
    NAMD extraBondsFile indices are 0-based."""
    old = psf_atoms(old_psf)
    new_index = {(s, r, a): i for i, (s, r, _, a) in enumerate(psf_atoms(new_psf))}
    n = 0
    with open(out, "w") as fh:
        for line in open(old_dat):
            p = line.split()
            if not p or p[0] != "DIHEDRAL":
                continue
            idx = [new_index[(old[int(i)][0], old[int(i)][1], old[int(i)][3])] for i in p[1:5]]
            fh.write(f"DIHEDRAL {idx[0]:6d} {idx[1]:6d} {idx[2]:6d} {idx[3]:6d} {' '.join(p[5:])}\n")
            n += 1
    return n


def slurm(kind, job, hours):
    step = "step4_equilibration" if kind == "eq" else "step5_production"
    log = "eq.log" if kind == "eq" else "prod.log"
    tag = "e" if kind == "eq" else "p"
    return f"""#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu     # partition name
#SBATCH --job-name={job}
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time={hours:02d}:00:00
#SBATCH --output=slurm{tag}-%j.out

{NAMD_BIN} {step}.inp > {log}
"""


def build(tk, sp):
    sysdir = BASE / tk / sp
    bdir, sdir = sysdir / "build", sysdir / "sim"
    for d in (bdir, sdir):
        shutil.rmtree(d, ignore_errors=True)
    bdir.mkdir(parents=True)
    scaf = Path(glob.glob(str(SCAFFOLD[tk] / "charmm-gui-*"))[0])
    for s in ("PROC", "PROD"):
        shutil.copy(sysdir / f"fgf_{s}.pdb", bdir)
    tcl = BUILD_TCL.format(toppar=scaf / "toppar", scaf_psf=scaf / "step1_pdbreader.psf",
                           scaf_pdb=scaf / "step1_pdbreader.pdb", box=BOX[tk])
    (bdir / "build.tcl").write_text(tcl)
    out = run_vmd("build.tcl", bdir)
    mm = re.search(r"SOLUTE_MINMAX \{(.*?)\} \{(.*?)\}", out)
    lo, hi = [list(map(float, g.split())) for g in mm.groups()]
    pad = min(min(BOX[tk] / 2 + l, BOX[tk] / 2 - h) for l, h in zip(lo, hi))

    sdir.mkdir()
    # NAMD-format parameters: the toppar the original replicates ran with. CHARMM-GUI's
    # top-level toppar/ is the CHARMM version (contains `set nat ?NATC` etc.), which
    # NAMD rejects with "UNKNOWN PARAMETER IN CHARMM PARAMETER FILE".
    shutil.copytree(scaf / "namdrep1" / "toppar", sdir / "toppar")
    pbc = (scaf / "step3_pbcsetup.str").read_text()
    for k in ("A", "B", "C"):
        pbc = re.sub(rf"SET {k} = \S+", f"SET {k} = {BOX[tk]:g}", pbc)
    (sdir / "step3_pbcsetup.str").write_text(pbc)
    scaf_rep = scaf / "namdrep1"
    chim = sp.startswith("chim_")              # segment-swap chimeras (make_chimeras.py): 3 replicates
    for r in range(1, (3 if chim else N_REPS) + 1):
        rd = sdir / f"namdrep{r}"
        (rd / "restraints").mkdir(parents=True)
        shutil.copy(bdir / "ionized.psf", rd / "step3_input.psf")
        shutil.copy(bdir / "ionized.pdb", rd / "step3_input.pdb")
        write_posres(bdir / "ionized.pdb", bdir / "ionized.psf", rd / "restraints" / "prot_posres.ref")
        n_dihe = remap_carb_restraints(scaf_rep / "restraints" / "carbohydrate_restraint.dat",
                                       scaf_rep / "step3_input.psf", bdir / "ionized.psf",
                                       rd / "restraints" / "carbohydrate_restraint.dat")
        seed = 1000 * (1 if tk == "1fq9" else 2) + 10 * r + (500 if chim else 0)   # chimeras: seeds distinct from the parents'
        for inp in ("step4_equilibration.inp", "step5_production.inp"):
            # parameters read from the shared sim/toppar/ (no per-replicate symlink:
            # links get dropped or broken by scp/SFTP/Globus transfers)
            txt = re.sub(r"(?m)^parameters( +)toppar/", r"parameters\1../toppar/", (scaf_rep / inp).read_text())
            txt = txt.replace("set temp                310;",
                              f"set temp                310;\nseed                    {seed + (inp[4] == '5')};", 1)
            assert f"seed" in txt
            (rd / inp).write_text(txt)
        job = f"{LABEL[tk]}_{sp}_r{r}"
        (rd / "eq_slurm.sh").write_text(slurm("eq", job + "_eq", 6))
        (rd / "prod_slurm.sh").write_text(slurm("prod", job + "_pr", 12))
        (rd / "submit.sh").write_text(
            "#!/bin/bash\n# equilibration, then production only if equilibration finished cleanly\n"
            "cd \"$(dirname \"$0\")\"\njid=$(sbatch --parsable eq_slurm.sh)\n"
            "sbatch --dependency=afterok:$jid prod_slurm.sh\n")
        os.chmod(rd / "submit.sh", 0o755)

    natom = len(psf_atoms(bdir / "ionized.psf"))
    ions = {k: sum(1 for a in psf_atoms(bdir / "ionized.psf") if a[2] == k) for k in ("POT", "CLA")}
    summary = (f"{tk}/{sp}: {natom} atoms, box {BOX[tk]:g} A, min solute-to-box padding {pad:.1f} A, "
               f"POT {ions['POT']} CLA {ions['CLA']}, carb dihedral restraints {n_dihe}")
    (sysdir / "build_summary.txt").write_text(summary + "\n")
    return summary


def main():
    targets = sys.argv[1:] or [f"{tk}/{sp.name}" for tk in ("1fq9", "1e0o") for sp in sorted((BASE / tk).iterdir())]
    for t in targets:
        tk, sp = t.split("/")
        print(build(tk, sp), flush=True)


if __name__ == "__main__":
    main()
