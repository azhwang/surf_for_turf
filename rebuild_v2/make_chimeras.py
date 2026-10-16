#!/usr/bin/env python3
"""
make_chimeras.py -- FGF2 segment-swap chimeras to test the Arg90 lead (2026-10-16).

Exploratory finding (ext_results/explo_summary.txt, PREDICTION_goldfish_rep3.md): in the 2:2:1
assembly the conserved Arg90 binds the shared heparin in goldfish but not in human FGF2. Arg90
itself is identical in both, but 8 of the 25 goldfish substitutions cluster around it, in
FGF2 86-97 (155-aa numbering): K86N E87S L92I A93G S94T K95R C96R V97T. The chimeras swap that
segment:
    chim_human_gfseg     human core, goldfish residues at 86-97
    chim_goldfish_hseg   goldfish core, human residues at 86-97

Construction: the host's placed core (<tmpl>/<host>/fgf_PROC.pdb, fgf_PROD.pdb, the same files
the parent systems were built from) keeps its backbone (N, CA, C, O) everywhere. At each swapped
position the donor's side chain (from the donor's placed core) is superposed on the host residue's
N/CA/C. Both cores have identical length and register (127 / 129 residues, no indels), so
positions map one to one. fix_rotamers.fix_system then removes any piercing or severe clash, and
precheck.py must pass before build_all.py builds the system exactly like the parents (3 replicates,
distinct seeds). Production extensions (50 x 1 ns, same as the parents) are prepared too and
chained after production by submit_chimeras.sh.

Usage: python3 make_chimeras.py [1e0o|1fq9 ...]      (default: both templates)
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from prepare_fgf_core import THREE, kabsch            # noqa: E402
from fix_rotamers import fix_system                    # noqa: E402
import make_extensions as mx                           # noqa: E402

SEGMENT = range(86, 98)                                 # FGF2 86-97, 155-aa numbering
OFF = {"1e0o": 27, "1fq9": 24}                          # FGF2 number = core position + OFF
CHIMERAS = {"chim_human_gfseg": ("human", "goldfish"), "chim_goldfish_hseg": ("goldfish", "human")}
N_REPS = 3
BB = {"N", "CA", "C", "O", "OT1", "OT2", "OXT"}


def residues(pdb):
    out, order = {}, []
    for l in open(pdb):
        if not l.startswith("ATOM"):
            continue
        rid = int(l[22:26])
        if rid not in out:
            out[rid] = [l[17:21].strip(), []]
            order.append(rid)
        out[rid][1].append((l[12:16].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return [(r, *out[r]) for r in order]


def swap(host, donor, tk):
    """host/donor: [(resid, resname, atoms)] of the same length. Returns new residues and the mutation list."""
    assert len(host) == len(donor), (len(host), len(donor))
    new, muts = [], []
    for k, (h, d) in enumerate(zip(host, donor), 1):
        fgf2 = k + OFF[tk]
        if fgf2 not in SEGMENT or THREE[h[1]] == THREE[d[1]]:
            new.append(h)
            continue
        hm, dm = dict(h[2]), dict(d[2])
        R, t = kabsch(np.array([dm[a] for a in ("N", "CA", "C")]), np.array([hm[a] for a in ("N", "CA", "C")]))
        atoms = [(n, hm[n]) for n, _ in h[2] if n in BB] + [(n, x @ R.T + t) for n, x in d[2] if n not in BB]
        new.append((h[0], d[1], atoms))
        muts.append(f"{THREE[h[1]]}{fgf2}{THREE[d[1]]}")
    return new, muts


def write_pdb(res, seg, path):
    with open(path, "w") as fh:
        n = 0
        for rid, rn, atoms in res:
            for name, (x, y, z) in atoms:
                n += 1
                fh.write(f"ATOM  {n:5d} {name:<4s} {rn:<4s}X{rid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      {seg:<4s}\n")
        fh.write("END\n")


def main():
    tks = sys.argv[1:] or ["1e0o", "1fq9"]
    for tk in tks:
        for name, (host, donor) in CHIMERAS.items():
            out = BASE / tk / name
            if out.exists():
                shutil.rmtree(out)
            out.mkdir()
            report = [f"chimera {name}: host {host} core (backbone everywhere), donor {donor} side chains at "
                      f"FGF2 {SEGMENT.start}-{SEGMENT.stop - 1}"]
            for seg in ("PROC", "PROD"):
                h = residues(BASE / tk / host / f"fgf_{seg}.pdb")
                d = residues(BASE / tk / donor / f"fgf_{seg}.pdb")
                new, muts = swap(h, d, tk)
                write_pdb(new, seg, out / f"fgf_{seg}.pdb")
                seq = "".join(THREE[r[1]] for r in new)
                report.append(f"{seg}: {len(muts)} substitutions {' '.join(muts)}; sequence {seq}")
            (out / "core_report.txt").write_text("\n".join(report) + "\n")
            print("\n".join(report))
    # same order as run_rebuild.sh: precheck (writes precheck.pdb) -> rotamer fix -> precheck must pass
    subprocess.run([sys.executable, str(BASE / "precheck.py")], capture_output=True, text=True)
    for tk in tks:
        for name in CHIMERAS:
            for n in fix_system(BASE / tk / name):
                print(f"{tk}/{name}: {n}")
    # precheck covers every system in the template folders (read-only apart from precheck.pdb)
    p = subprocess.run([sys.executable, str(BASE / "precheck.py")], capture_output=True, text=True)
    chk = [l for l in p.stdout.splitlines() if "chim_" in l or l.startswith("    ")]
    print("\n".join(chk))
    if any("chim_" in l and ("piercings 0" not in l or "severe<1.8A 0" not in l) for l in chk):
        raise SystemExit("chimera precheck failed -- not building")
    for tk in tks:
        for name in CHIMERAS:
            print(subprocess.run([sys.executable, str(BASE / "build_all.py"), f"{tk}/{name}"],
                                 capture_output=True, text=True, check=True).stdout.strip())
            sim = BASE / tk / name / "sim"
            for r in range(1, N_REPS + 1):
                d = sim / f"namdrep{r}"
                prod = (d / "step5_production.inp").read_text()
                seed = int(re.search(r"(?m)^seed\s+(\d+);", prod).group(1)) * 100
                (d / "ext_template.inp").write_text(mx.make_template(prod))
                job = f"ext_{mx.LABEL[tk]}_{name}_r{r}"
                (d / "ext_slurm.sh").write_text(mx.SLURM.format(job=job, total=mx.TOTAL_NS, seg=mx.SEG_STEPS,
                                                                seed=seed, namd=mx.NAMD))
                (d / "ext_slurm.sh").chmod(0o755)
                (d / "submit_chain.sh").write_text(
                    "#!/bin/bash\n# equilibration -> 1 ns production -> 50 x 1 ns extension chain, each only if the "
                    "previous step finished cleanly\ncd \"$(dirname \"$0\")\"\n"
                    "e=$(sbatch --parsable eq_slurm.sh)\np=$(sbatch --parsable --dependency=afterok:$e prod_slurm.sh)\n"
                    "sbatch --dependency=afterok:$p ext_slurm.sh\n")
                (d / "submit_chain.sh").chmod(0o755)


if __name__ == "__main__":
    main()
