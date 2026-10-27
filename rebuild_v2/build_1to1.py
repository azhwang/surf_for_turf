#!/usr/bin/env python3
"""
build_1to1.py -- 1:1 FGF2-FGFR systems for all seven orthologs, both receptors.

Built on the dedicated 1:1 crystal structures, 1EVT (FGF1-FGFR1) and 1DJS (FGF1-FGFR2), rather
than on the 2:2:2 / 2:2:1 assemblies. A 2:2 receptor conformation is not the right starting
point for a 1:1 simulation: 1DJS and 1E0O differ by a large D2-D3 hinge rotation (whole-chain
CA RMSD 12.2 A, D2 1.0 A, D3 1.7 A) and engage FGF very differently (575 vs 244 heavy-atom
contacts within 4.5 A). 1EVT and the 1FQ9 receptor are equivalent (0.87 A), so for FGFR1 the
two routes agree, but using both crystals keeps one protocol for both receptors.

The receptor chains of both crystals have a 10-12 residue gap in the same D3 loop; those are
closed with the structurally equivalent loop of the complete 1FQ9 chain C by
prepare_1to1_scaffold.py, and prepare_1to1_fgf.py places each ortholog's FGF2 core on the
crystal FGF chain. Run both before this script (run_1to1.sh does all three in order).

Per system this script:
  1. builds the PSF with psfgen (segments PROA = FGFR, PROC = FGF2, both D3 disulfides patched,
     missing side-chain atoms and all hydrogens from guesscoord);
  2. relaxes the grafted loops with NAMD in implicit solvent, with every atom outside the
     loops and their 5 A surroundings held fixed. This closes the one remaining mid-loop break
     per receptor (1.96 A in FGFR1, 2.89 A in FGFR2) and the few overlaps the rigid graft leaves
     against FGF, while leaving all crystal coordinates elsewhere untouched;
  3. re-solvates at 12 A padding, neutralises to 0.15 M KCl, and writes N_REPS replicate folders
     with NAMD inputs and slurm scripts. No heparin, so the carbohydrate dihedral restraints are
     dropped.

  fgfr1 <- 1EVT chain C + chain A      fgfr2 <- 1DJS chain A + chain B

Usage:  python3 build_1to1.py [rec/species ...]      (default: all 14)
"""
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
OUT = BASE / "onetoone"
VMD_DIR = "/Applications/VMD 2.0.0a7-pre2.app/Contents/vmd2/lib"
NAMD_BIN = "/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2"
# the grafted D3 loops, relaxed in step 2 (prepare_1to1_scaffold.py closes these gaps)
LOOPS = {"fgfr1": [(294, 303), (312, 323)], "fgfr2": [(297, 306)]}
FREE_PAD = 5.0                  # A: residues this close to a grafted loop also move
CLASH = 2.8                     # A: FGFR-FGF heavy-atom contacts this close are relaxed too
CLASH_PAD = 4.0                 # A: and residues this close to a clashing one
TOPPAR_SRC = "1fq9/goldfish/sim/toppar_charmm"
# only for copying the NAMD input settings and the toppar tree, so the 1:1 runs use exactly the
# same protocol as the corresponding assembly set
DONOR_TMPL = {"fgfr1": "1fq9", "fgfr2": "1e0o"}
NAMD_LOCAL = Path.home() / "Desktop/LAMM/Surf_for_Turf/NAMD_3.0.2_MacOS-universal-multicore/namd3"
LABEL = {"fgfr1": "11r1", "fgfr2": "11r2"}
SPECIES = ["human", "human_stab", "bovine_stab", "chinook_salmon", "goldfish",
           "japanese_ricefish", "zebrafish"]
N_REPS = 5
PAD = 12.0                      # A of water between solute and box edge
BACKBONE = {"N", "CA", "C", "O"}

PSFGEN_TCL = r"""
package require psfgen
resetpsf
topology {top}
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
segment PROA {{ first NTER; last CTER; pdb {rec_pdb} }}
coordpdb {rec_pdb} PROA
segment PROC {{ first NTER; last CTER; pdb {fgf_pdb} }}
coordpdb {fgf_pdb} PROC
{disu}
regenerate angles dihedrals
guesscoord
writepsf solute_raw.psf
writepdb solute_raw.pdb
quit
"""

MIN_CONF = """\
# local relaxation of the grafted D3 loops; everything else is held at its crystal position
structure               solute_raw.psf
coordinates             solute_raw.pdb
paraTypeCharmm          on
parameters              {par}
outputName              min
exclude                 scaled1-4
1-4scaling              1.0
cutoff                  16.0
switching               on
switchdist              15.0
pairlistdist            18.0
GBIS                    on
alphaCutoff             14.0
ionConcentration        0.15
solventDielectric       78.5
fixedAtoms              on
fixedAtomsFile          free.pdb
fixedAtomsCol           B
temperature             0
minimize                {steps}
"""

BUILD_TCL = r"""
package require psfgen
package require solvate
package require autoionize
resetpsf
readpsf {solute_psf}
coordpdb {solute_pdb}
writepsf pair.psf
writepdb pair.pdb
set h [expr {box}/2.0]
solvate pair.psf pair.pdb -minmax [list [list -$h -$h -$h] [list $h $h $h]] -o solvated
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


def orient(src, dst):
    """Rewrite a solute PDB rotated onto its principal axes and centred on the origin.

    The solute tumbles freely in water, so its orientation in the crystal frame is arbitrary,
    but the box is cubic and sized from the axis-aligned extent. In the 1EVT and 1DJS frames the
    long axis of the receptor runs diagonally, so measuring without reorienting would either
    inflate the box or, worse, leave the longest dimension of the complex (about 101 A) with
    less clearance than PAD once the complex rotates. Aligning first makes the box depend only
    on the shape of the complex: side = longest principal extent + 2 * PAD.
    """
    rows = [l for l in open(src) if l.startswith(("ATOM", "HETATM"))]
    X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rows])
    heavy = np.array([not l[12:16].strip().startswith("H") for l in rows])
    C = X - X[heavy].mean(0)
    _, V = np.linalg.eigh(C[heavy].T @ C[heavy])
    # eigh returns an orthonormal basis that may be left-handed; using it as given mirrors the
    # structure into D-amino acids, which NAMD runs happily and silently
    if np.linalg.det(V) < 0:
        V[:, 0] = -V[:, 0]
    P = C @ V                                   # columns ordered by ascending extent
    P -= (P[heavy].max(0) + P[heavy].min(0)) / 2
    with open(dst, "w") as fh:
        for l, (x, y, z) in zip(rows, P):
            fh.write(f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}")
        fh.write("END\n")
    assert chirality_violations(dst) == 0, f"{dst} is mirrored"
    return float((P[heavy].max(0) - P[heavy].min(0)).max())


def chirality_violations(pdb):
    """Residues whose CB sits on the D side of the N-CA-C plane. Must be zero."""
    res = {}
    for l in open(pdb):
        if l.startswith(("ATOM", "HETATM")):
            res.setdefault((l[72:76].strip(), int(l[22:26])), {})[l[12:16].strip()] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return sum(1 for a in res.values() if {"N", "CA", "C", "CB"} <= set(a)
               and np.dot(np.cross(a["N"] - a["CA"], a["C"] - a["CA"]), a["CB"] - a["CA"]) < 0)


def minimise_loops(bdir, rec, par, steps=2000):
    """Relax the grafted loops and their surroundings; everything else stays fixed.

    free.pdb carries 0.00 in the B column for atoms that may move and 1.00 for the rest
    (fixedAtomsCol B treats non-zero as fixed).
    """
    rows = [l for l in open(bdir / "solute_raw.pdb") if l.startswith(("ATOM", "HETATM"))]
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rows])
    key = [(l[72:76].strip(), int(l[22:26])) for l in rows]
    heavy = np.array([not l[12:16].strip().startswith("H") for l in rows])
    loop = np.array([k[0] == "PROA" and any(a <= k[1] <= b for a, b in LOOPS[rec]) for k in key])
    near = np.linalg.norm(xyz[:, None, :] - xyz[None, loop, :], axis=2).min(1) <= FREE_PAD
    # the rigid placement of the FGF2 core also leaves a few tight side-chain contacts across
    # the interface (1.1-2.4 A here), the same thing fix_rotamers.py handles in the assembly
    # rebuild; let those residues and their immediate surroundings relax as well
    a = heavy & np.array([k[0] == "PROA" for k in key])
    c = heavy & np.array([k[0] == "PROC" for k in key])
    d = np.linalg.norm(xyz[a][:, None, :] - xyz[c][None, :, :], axis=2)
    clash = np.zeros(len(rows), bool)
    clash[np.flatnonzero(a)[np.where(d < CLASH)[0]]] = True
    clash[np.flatnonzero(c)[np.where(d < CLASH)[1]]] = True
    if clash.any():
        clash |= np.linalg.norm(xyz[:, None, :] - xyz[None, clash, :], axis=2).min(1) <= CLASH_PAD
    free_res = {k for k, f in zip(key, loop | near | clash) if f}
    free = np.array([k in free_res for k in key])          # whole residues, not single atoms
    with open(bdir / "free.pdb", "w") as fh:
        for l, f in zip(rows, free):
            fh.write(f"{l[:54]}  1.00{0.0 if f else 1.0:6.2f}{l[66:]}")
        fh.write("END\n")
    (bdir / "min.conf").write_text(MIN_CONF.format(par=par, steps=steps))
    p = subprocess.run([str(NAMD_LOCAL), "+p4", "min.conf"], cwd=bdir, capture_output=True, text=True)
    (bdir / "min.log").write_text(p.stdout + p.stderr)
    if not (bdir / "min.coor").exists():
        raise SystemExit(f"loop minimisation failed in {bdir}; see min.log")
    return int(free.sum()), len(free_res)


def geometry_report(psf, pdb, rec):
    """Peptide-bond lengths at the graft junctions, worst close contact, ring piercings."""
    atoms = psf_atoms(psf)
    rows = [l for l in open(pdb) if l.startswith(("ATOM", "HETATM"))]
    crd = {}
    for (seg, resid, _, name), l in zip(atoms, rows):
        if not name.startswith("H"):
            crd[(seg, int(resid), name)] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    res = sorted({(s, r) for s, r, _ in crd})
    gaps = []
    for (s1, r1), (s2, r2) in zip(res, res[1:]):
        if s1 != s2 or r2 != r1 + 1:
            continue
        d = np.linalg.norm(crd[(s1, r1, "C")] - crd[(s2, r2, "N")])
        if d > 1.6:
            gaps.append(f"{s1}:{r1}-{r2} {d:.2f}A")
    loopset = {r for a, b in LOOPS[rec] for r in range(a, b + 1)}
    keys = list(crd)
    X = np.array([crd[k] for k in keys])
    li = [i for i, k in enumerate(keys) if k[0] == "PROA" and k[1] in loopset]
    d = np.linalg.norm(X[li][:, None, :] - X[None, :, :], axis=2)
    for a, i in enumerate(li):
        for b, k in enumerate(keys):
            if k[0] == keys[i][0] and abs(k[1] - keys[i][1]) <= 1:
                d[a, b] = 9.9
    ia = [i for i, k in enumerate(keys) if k[0] == "PROA"]
    ic = [i for i, k in enumerate(keys) if k[0] == "PROC"]
    iface = float(np.linalg.norm(X[ia][:, None, :] - X[ic][None, :, :], axis=2).min())
    out = subprocess.run([sys.executable, str(BASE.parent / "check_ring_piercing.py"), str(pdb)],
                         capture_output=True, text=True).stdout.splitlines()
    return gaps, float(d.min()), iface, max(0, len(out) - 1)


def build(rec, sp):
    tk = DONOR_TMPL[rec]
    fgf = OUT / rec / sp / "fgf_PROC.pdb"
    recpdb = OUT / "scaffold" / f"receptor_{rec}.pdb"
    for f in (fgf, recpdb):
        if not f.exists():
            raise SystemExit(f"missing {f} -- run prepare_1to1_scaffold.py and "
                             f"prepare_1to1_fgf.py first")
    sysdir = OUT / rec / sp
    bdir, sdir = sysdir / "build", sysdir / "sim"
    for d in (bdir, sdir):
        shutil.rmtree(d, ignore_errors=True)
    bdir.mkdir(parents=True)

    def pick(name):
        own = BASE / tk / sp / "sim" / name
        if own.exists():
            return own
        for other in sorted((BASE / tk).glob("*")):
            if (other / "sim" / name).exists():
                return other / "sim" / name
        raise SystemExit(f"no {name} found for {tk}")

    toppar_charmm = BASE / TOPPAR_SRC
    disu = "\n".join(f"patch DISU PROA:{a} PROA:{b}" for a, b in
                      (l.split() for l in
                       (OUT / "scaffold" / f"disu_{rec}.txt").read_text().split("\n") if l.strip()))
    (bdir / "psfgen.tcl").write_text(PSFGEN_TCL.format(
        top=toppar_charmm / "top_all36_prot.rtf", rec_pdb=recpdb, fgf_pdb=fgf, disu=disu))
    log = run_vmd("psfgen.tcl", bdir)
    if not (bdir / "solute_raw.psf").exists():
        raise SystemExit(f"psfgen failed in {bdir}; see vmd_build.log")

    nfree, nres = minimise_loops(bdir, rec, toppar_charmm / "par_all36m_prot.prm")
    shutil.copy(bdir / "solute_raw.psf", bdir / "solute.psf")
    # min.coor is a NAMD binary restart; convert it back to PDB on the raw PSF
    (bdir / "tocoor.tcl").write_text(
        f'mol new solute_raw.psf\nmol addfile min.coor waitfor all\n'
        f'[atomselect top all] writepdb solute.pdb\nquit\n')
    run_vmd("tocoor.tcl", bdir)
    if not (bdir / "solute.pdb").exists():
        raise SystemExit(f"could not convert min.coor in {bdir}")
    gaps, worst, iface, pierce = geometry_report(bdir / "solute.psf", bdir / "solute.pdb", rec)

    extent = orient(bdir / "solute.pdb", bdir / "solute_oriented.pdb")
    box = float(int(extent + 2 * PAD) + 1)
    (bdir / "build.tcl").write_text(BUILD_TCL.format(
        solute_psf=bdir / "solute.psf", solute_pdb=bdir / "solute_oriented.pdb", box=box))
    run_vmd("build.tcl", bdir)
    if not (bdir / "ionized.psf").exists():
        raise SystemExit(f"VMD build failed in {bdir}; see vmd_build.log")

    sdir.mkdir()
    scaf_rep = donor_inputs(tk, sp)
    shutil.copytree(pick("toppar"), sdir / "toppar")
    shutil.copytree(toppar_charmm, sdir / "toppar_charmm")
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
            txt = strip_carb_restraints((scaf_rep / inp).read_text())
            txt = re.sub(r"(?m)^parameters( +)(\.\./)?toppar/", r"parameters\1../toppar/", txt)
            # the donor input already carries a seed from the main rebuild; replace it
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

    assert chirality_violations(bdir / "ionized.pdb") == 0, f"{bdir}/ionized.pdb is mirrored"
    atoms = psf_atoms(bdir / "ionized.psf")
    ions = {k: sum(1 for a in atoms if a[2] == k) for k in ("POT", "CLA")}
    segs = sorted({a[0] for a in atoms if a[0].startswith("PRO")})
    summary = (f"{rec}/{sp}: {len(atoms)} atoms, box {box:g} A, segments {','.join(segs)}, "
               f"POT {ions['POT']} CLA {ions['CLA']}; loop relaxation freed {nres} residues "
               f"({nfree} atoms), peptide gaps {gaps or 'none'}, closest loop contact "
               f"{worst:.2f} A, closest FGFR-FGF contact {iface:.2f} A, ring piercings {pierce}")
    (sysdir / "build_summary.txt").write_text(summary + "\n")
    return summary


def main():
    targets = sys.argv[1:] or [f"{r}/{s}" for r in LOOPS for s in SPECIES]
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
