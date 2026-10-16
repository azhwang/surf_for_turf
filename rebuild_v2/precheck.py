#!/usr/bin/env python3
"""precheck.py -- combine scaffold + placed FGF heavy atoms and report ring
piercing, SEVERE clashes (FGF side-chain heavy atom < 1.8 A from a scaffold or
non-adjacent FGF heavy atom), milder FGF-scaffold contacts < 2.2 A, and FGF
peptide-bond gaps (C(i)-N(i+1) > 1.6 A). Exit status 1 if any piercing or
severe clash remains. Run after prepare_fgf_core.py, before building."""
import collections, glob, math, subprocess, sys
from pathlib import Path
BASE = Path(__file__).resolve().parent
SCAF = {"1fq9": "../all_metrics_analysis/1fq9_goldfish", "1e0o": "../all_metrics_analysis_221/1e0o_chinook_salmon"}
bad = 0
for tk in ("1fq9", "1e0o"):
    scaf = glob.glob(str(BASE / SCAF[tk] / "charmm-gui-*" / "step1_pdbreader.pdb"))[0]
    lines = [l for l in open(scaf) if l.startswith("ATOM") and l[72:76].strip().startswith(("PROA", "PROB", "CAR"))
             and not l[12:16].strip().startswith("H")]
    for sd in sorted(d for d in (BASE / tk).iterdir() if (d / "fgf_PROC.pdb").exists()):   # skip .DS_Store etc.
        fgf = [l for s in ("PROC", "PROD") for l in open(sd / f"fgf_{s}.pdb") if l.startswith("ATOM")]
        pre = sd / "precheck.pdb"
        pre.write_text("".join(lines + fgf))
        pierce = subprocess.run([sys.executable, str(BASE.parent / "check_ring_piercing.py"), str(pre)],
                                capture_output=True, text=True).stdout.splitlines()
        atoms = [(l[72:76].strip(), (float(l[30:38]), float(l[38:46]), float(l[46:54]))) for l in lines + fgf]
        grid = collections.defaultdict(list)
        for s, p in atoms:
            grid[tuple(int(c // 2.5) for c in p)].append((s, p))
        clash = 0
        for s, p in atoms:
            if s in ("PROC", "PROD"):
                g = tuple(int(c // 2.5) for c in p)
                clash += sum(1 for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
                             for s2, p2 in grid.get((g[0] + dx, g[1] + dy, g[2] + dz), [])
                             if s2 not in ("PROC", "PROD") and math.dist(p, p2) < 2.2)
        fgf_rows = [(l[72:76].strip(), int(l[22:26]), l[12:16].strip(),
                     (float(l[30:38]), float(l[38:46]), float(l[46:54]))) for l in fgf]
        scaf_rows = [(l[72:76].strip(), int(l[22:26]), l[12:16].strip(),
                      (float(l[30:38]), float(l[38:46]), float(l[46:54]))) for l in lines]
        g2 = collections.defaultdict(list)
        for r in scaf_rows + fgf_rows:
            g2[tuple(int(c // 2) for c in r[3])].append(r)
        severe = set()
        for r in fgf_rows:
            if r[2] in ("N", "CA", "C", "O"):
                continue
            g = tuple(int(c // 2) for c in r[3])
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for o in g2.get((g[0] + dx, g[1] + dy, g[2] + dz), []):
                            if o[0] == r[0] and abs(o[1] - r[1]) <= 1:
                                continue
                            if math.dist(r[3], o[3]) < 1.8:
                                severe.add(f"{r[0]}:{r[1]}:{r[2]}~{o[0]}:{o[1]}:{o[2]}")
        gaps = []
        for s in ("PROC", "PROD"):
            res = collections.OrderedDict()
            for l in fgf:
                if l[72:76].strip() == s:
                    res.setdefault(int(l[22:26]), {})[l[12:16].strip()] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            ids = list(res)
            gaps += [f"{s}:{a}-{b} {math.dist(res[a]['C'], res[b]['N']):.2f}A" for a, b in zip(ids, ids[1:])
                     if math.dist(res[a]["C"], res[b]["N"]) > 1.6]
        n_p = len(pierce) - 1
        bad += n_p > 0 or bool(severe)
        print(f"{tk}/{sd.name:18s} piercings {n_p}  severe<1.8A {len(severe)}  clashes<2.2A {clash:3d}  "
              f"peptide gaps {gaps or 'none'}")
        for s_ in sorted(severe):
            print("    severe:", s_)
        for p in pierce[1:]:
            print("   ", p.strip())
sys.exit(1 if bad else 0)
