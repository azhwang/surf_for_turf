#!/usr/bin/env python3
"""
derive_mmpbsa_split.py

Derives the MM/GBSA receptor/ligand split for an FGF-FGFR-heparin assembly
directly from THIS replicate's own CHARMM-GUI PSF/PDB, instead of using
hardcoded residue ranges. Used by run_mmpbsa_pipeline_v2.sh.

Why this exists: the v1 pipeline hardcoded 1FQ9-goldfish residue ranges and
reused them for every species and for the 1E0O (2:2:1) systems, where the
chain lengths differ. It also wrote interface-2's ligand mask in
full-system numbering, but ante-MMPBSA.py applies -n AFTER the -s strip has
renumbered residues, so interface 2's "ligand" was part of an FGF chain plus
heparin in every system.

What it does:
  1. Reads the PSF's residue order. Chamber numbers residues sequentially in
     PSF order, so this is the global numbering. Checks it against the
     chamber prmtop's residue labels.
  2. Classifies protein segments by sequence:
       FGF  = contains the beta-trefoil motif "CKNGG"
       FGFR = contains the D2 motif "PYW"
     Heparin = every CAR* segment. Solvent/ions are recognised by resname.
  3. Pairs each FGF with the FGFR it contacts most (8 A CA-CA contacts in
     the starting PDB). Aborts if the pairing is not unambiguous.
  4. Heparin:
       - one CAR segment (2:2:1)  -> shared, kept in BOTH interfaces
       - several CAR segments     -> each fragment goes to the interface
                                     (FGF+FGFR pair) it contacts most;
                                     aborts if ambiguous
     Heparin is always grouped with the FGFR (receptor) side.
     FGF is always the ligand.
  5. Computes masks in the numbering each AmberTools step actually sees:
       - parmed/cpptraj `strip` masks: full-system numbering
       - ante-MMPBSA.py -n ligand mask: numbering AFTER stripping the
         other interface AND solvent
     and checks the ligand is contiguous in that numbering.

Output: shell variable assignments on stdout (for `eval`), and a readable
split_summary.txt. With --check, instead verifies that the component
prmtops ante-MMPBSA.py produced contain the expected residues.

Overrides (environment variables):
  PAIRING="PROC:PROA,PROD:PROB"       force FGF:FGFR pairing
  HEPARIN_ASSIGN="CARA:1,CARB:2"      force heparin fragment -> interface
  HEPARIN_SHARED=1                    keep all heparin in both interfaces
"""

import argparse
import math
import os
import sys
from collections import OrderedDict

SOLVENT_RESNAMES = {"TIP3", "WAT", "HOH", "SPC"}
ION_RESNAMES = {"POT", "CLA", "SOD", "MG", "CAL", "ZN2", "CES", "LIT", "RUB"}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "HSD": "H", "HSE": "H", "HSP": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
FGF_MOTIF = "CKNGG"
FGFR_MOTIF = "PYW"
CONTACT_CUTOFF = 8.0     # Angstrom, between CA / heparin ring+sulfur atoms
PAIR_MARGIN = 1.5        # best partner must have >= 1.5x the contacts of the runner-up
MIN_PAIR_CONTACTS = 20
HEPARIN_CONTACT_ATOMS = {"C1", "C2", "C3", "C4", "C5", "S", "S2", "S6", "SN", "S1", "S3"}


class SplitError(Exception):
    pass


class AmbiguousInterface(SplitError):
    """Pairing/heparin assignment is ambiguous: per-interface split is not
    meaningful, but the whole-assembly split (which needs no pairing) is."""


def read_psf_residues(psf_path):
    """Ordered residues as (segid, resid, resname), in PSF (= chamber) order."""
    residues = []
    last_key = None
    with open(psf_path) as f:
        for line in f:
            if "!NATOM" in line:
                natom = int(line.split()[0])
                break
        else:
            raise SplitError(f"{psf_path}: no !NATOM section")
        for _ in range(natom):
            parts = next(f).split()
            segid, resid, resname = parts[1], parts[2], parts[3]
            key = (segid, resid)
            if key != last_key:
                residues.append((segid, resid, resname))
                last_key = key
    return residues


def read_prmtop_labels(prmtop_path):
    with open(prmtop_path) as f:
        text = f.read()
    block = text.split("%FLAG RESIDUE_LABEL")[1].split("%FLAG")[0]
    body = block.split("\n", 2)[2].replace("\n", "")
    return [body[i:i + 4].strip() for i in range(0, len(body), 4) if body[i:i + 4].strip()]


def read_pdb_contact_atoms(pdb_path):
    """{segid: [(x,y,z), ...]} using CA for protein and ring/sulfur atoms for sugars."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            segid = line[72:76].strip()
            atom = line[12:16].strip()
            if (segid.startswith("PRO") and atom == "CA") or \
               (segid.startswith("CAR") and atom in HEPARIN_CONTACT_ATOMS):
                coords.setdefault(segid, []).append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


def count_contacts(a, b):
    return sum(1 for p in a for q in b if math.dist(p, q) < CONTACT_CUTOFF)


def ranges_to_mask(indices):
    """Sorted 1-based residue indices -> ':a-b,c-d' mask string."""
    if not indices:
        return ""
    indices = sorted(indices)
    spans = []
    start = prev = indices[0]
    for i in indices[1:]:
        if i == prev + 1:
            prev = i
            continue
        spans.append((start, prev))
        start = prev = i
    spans.append((start, prev))
    return ":" + ",".join(f"{s}-{e}" if s != e else f"{s}" for s, e in spans)


def renumber(all_indices, removed):
    """Map full-system index -> index after removing `removed` (1-based)."""
    mapping = {}
    n = 0
    for i in all_indices:
        if i in removed:
            continue
        n += 1
        mapping[i] = n
    return mapping


def derive(psf, pdb, prmtop):
    residues = read_psf_residues(psf)
    n_res = len(residues)

    # --- check PSF order == chamber numbering ---
    if prmtop:
        labels = read_prmtop_labels(prmtop)
        if len(labels) != n_res:
            raise SplitError(f"residue count mismatch: PSF {n_res} vs {prmtop} {len(labels)}")
        for i, ((seg, rid, rn), lab) in enumerate(zip(residues, labels), start=1):
            if rn in SOLVENT_RESNAMES or rn in ION_RESNAMES:
                continue
            if lab != rn[:4]:
                raise SplitError(f"residue {i} ({seg}:{rid}) is {rn} in PSF but {lab} in {prmtop}")

    # --- segments -> global residue indices ---
    seg_res = OrderedDict()
    seg_seq = {}
    solvent_idx = set()
    solvent_names = set()
    for i, (seg, rid, rn) in enumerate(residues, start=1):
        if rn in SOLVENT_RESNAMES or rn in ION_RESNAMES:
            solvent_idx.add(i)
            solvent_names.add("WAT" if rn in SOLVENT_RESNAMES else rn)
            continue
        seg_res.setdefault(seg, []).append(i)
        seg_seq[seg] = seg_seq.get(seg, "") + THREE_TO_ONE.get(rn, "x")

    fgf, fgfr, heparin = [], [], []
    for seg in seg_res:
        if seg.startswith("CAR"):
            heparin.append(seg)
        elif seg.startswith("PRO"):
            s = seg_seq[seg]
            is_fgf, is_fgfr = FGF_MOTIF in s, FGFR_MOTIF in s
            if is_fgf == is_fgfr:
                raise SplitError(f"cannot classify {seg} (FGF motif: {is_fgf}, FGFR motif: {is_fgfr})")
            (fgf if is_fgf else fgfr).append(seg)
        else:
            raise SplitError(f"unexpected non-solvent segment {seg}")
    if len(fgf) != 2 or len(fgfr) != 2 or not heparin:
        raise SplitError(f"expected 2 FGF + 2 FGFR + heparin, found FGF={fgf} FGFR={fgfr} heparin={heparin}")

    coords = read_pdb_contact_atoms(pdb)
    contacts = {(a, b): count_contacts(coords.get(a, []), coords.get(b, []))
                for a in seg_res for b in seg_res if a != b}
    notes = []

    fgf_sorted = sorted(fgf)
    try:
        ifaces, hep_assign = pair_interfaces(fgf, fgfr, heparin, contacts, notes)
    except AmbiguousInterface as e:
        ifaces, hep_assign = None, None
        notes.append(f"PER-INTERFACE SPLIT SKIPPED -- {e}")
        notes.append("only the whole-assembly split (both FGF vs both FGFR + all heparin) is computed")
    return build_masks(residues, seg_res, solvent_idx, solvent_names, fgf_sorted, fgfr,
                       heparin, ifaces, hep_assign, contacts, notes, psf, pdb, prmtop)


def pair_interfaces(fgf, fgfr, heparin, contacts, notes):
    # --- FGF <-> FGFR pairing ---
    if os.environ.get("PAIRING"):
        pairs = dict(p.split(":") for p in os.environ["PAIRING"].split(","))
        notes.append(f"pairing forced by PAIRING={os.environ['PAIRING']}")
    else:
        pairs = {}
        for g in fgf:
            ranked = sorted(fgfr, key=lambda r: contacts[(g, r)], reverse=True)
            best, second = contacts[(g, ranked[0])], contacts[(g, ranked[1])]
            if best < MIN_PAIR_CONTACTS or best < PAIR_MARGIN * second:
                raise AmbiguousInterface(
                    f"ambiguous pairing for {g}: " +
                    ", ".join(f"{r}={contacts[(g, r)]}" for r in ranked) +
                    f" contacts (need >= {MIN_PAIR_CONTACTS} and >= {PAIR_MARGIN}x runner-up). "
                    "Inspect the structure; override with PAIRING=\"PROC:PROA,PROD:PROB\" only if justified.")
            pairs[g] = ranked[0]
    if sorted(pairs) != sorted(fgf) or sorted(pairs.values()) != sorted(fgfr):
        raise SplitError(f"pairing is not one-to-one: {pairs}")
    fgf_sorted = sorted(fgf)
    ifaces = [(g, pairs[g]) for g in fgf_sorted]   # iface1 = lower FGF segid

    # --- heparin assignment ---
    hep_assign = {}   # seg -> set of interface numbers (1, 2)
    if os.environ.get("HEPARIN_ASSIGN"):
        for p in os.environ["HEPARIN_ASSIGN"].split(","):
            seg, k = p.split(":")
            hep_assign[seg] = {int(k)}
        notes.append(f"heparin assignment forced by HEPARIN_ASSIGN={os.environ['HEPARIN_ASSIGN']}")
    elif len(heparin) == 1 or os.environ.get("HEPARIN_SHARED") == "1":
        for h in heparin:
            hep_assign[h] = {1, 2}
        notes.append("heparin treated as SHARED: kept in both interfaces (2:2:1 convention)")
    else:
        for h in heparin:
            score = [contacts[(h, g)] + contacts[(h, r)] for g, r in ifaces]
            hi, lo = max(score), min(score)
            if hi == 0 or hi < PAIR_MARGIN * lo:
                raise AmbiguousInterface(
                    f"ambiguous heparin fragment {h}: interface contacts {score}. "
                    "Override with HEPARIN_ASSIGN=\"CARA:1,...\" or HEPARIN_SHARED=1 if justified.")
            hep_assign[h] = {score.index(hi) + 1}
    if set(hep_assign) != set(heparin):
        raise SplitError(f"heparin assignment {hep_assign} does not cover {heparin}")
    return ifaces, hep_assign


def build_masks(residues, seg_res, solvent_idx, solvent_names, fgf_sorted, fgfr,
                heparin, ifaces, hep_assign, contacts, notes, psf, pdb, prmtop):
    all_idx = list(range(1, len(residues) + 1))
    out = OrderedDict()
    summary = [f"PSF: {psf}", f"PDB: {pdb}", f"prmtop checked: {prmtop or 'no'}",
               f"FGF (ligand): {fgf_sorted}", f"FGFR: {sorted(fgfr)}", f"heparin: {heparin}"]
    summary += [f"note: {n}" for n in notes]
    summary.append("contacts (8 A): " + ", ".join(
        f"{g}-{r}={contacts[(g, r)]}" for g in fgf_sorted for r in sorted(fgfr)))

    strip_mask = ":" + ",".join(sorted(solvent_names))
    out["STRIP_MASK"] = strip_mask

    def component(name, lig_segs, rec_segs, other_segs):
        other_idx = {i for s in other_segs for i in seg_res[s]}
        lig_idx = [i for s in lig_segs for i in seg_res[s]]
        rec_idx = [i for s in rec_segs for i in seg_res[s]]
        # numbering as ante-MMPBSA.py sees it for -n: after removing the
        # other interface (already gone from -p) AND solvent (-s)
        dry = renumber(all_idx, other_idx | solvent_idx)
        lig_dry = sorted(dry[i] for i in lig_idx)
        if lig_dry != list(range(lig_dry[0], lig_dry[-1] + 1)):
            raise SplitError(f"{name}: ligand not contiguous after stripping ({ranges_to_mask(lig_dry)})")
        out[f"{name.upper()}_OTHER_MASK"] = ranges_to_mask(sorted(other_idx))
        out[f"{name.upper()}_LIGAND_MASK"] = ranges_to_mask(lig_dry)
        out[f"{name.upper()}_N_LIGAND"] = str(len(lig_idx))
        out[f"{name.upper()}_N_RECEPTOR"] = str(len(rec_idx))
        summary.append(
            f"{name}: ligand {lig_segs} ({len(lig_idx)} res) | receptor {rec_segs} ({len(rec_idx)} res) | "
            f"stripped {sorted(other_segs) or '-'} | ante -n {ranges_to_mask(lig_dry)}"
            + (f" | strip (full numbering) {ranges_to_mask(sorted(other_idx))}" if other_idx else ""))

    out["IFACE_SPLIT"] = "1" if ifaces else "0"
    for k, (g, r) in enumerate(ifaces or [], start=1):
        other_k = 2 if k == 1 else 1
        og, orr = ifaces[other_k - 1]
        own_hep = [h for h in heparin if k in hep_assign[h]]
        other_hep = [h for h in heparin if k not in hep_assign[h]]
        component(f"iface{k}", [g], [r] + own_hep, [og, orr] + other_hep)

    component("whole", fgf_sorted, sorted(fgfr) + heparin, [])
    return out, summary


def check(prmtop_dir, expected):
    """Verify ante-MMPBSA.py produced components with the expected residue counts."""
    bad = []
    names = ("iface1", "iface2", "whole") if expected["IFACE_SPLIT"] == "1" else ("whole",)
    for name in names:
        for part, key in (("ligand", "N_LIGAND"), ("receptor", "N_RECEPTOR")):
            fn = os.path.join(prmtop_dir, f"{name}_{part}.parm7")
            n = len(read_prmtop_labels(fn))
            want = int(expected[f"{name.upper()}_{key}"])
            status = "ok" if n == want else "MISMATCH"
            print(f"  {name}_{part}.parm7: {n} residues (expected {want}) {status}")
            if n != want:
                bad.append(fn)
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psf", default="step3_input.psf")
    ap.add_argument("--pdb", default="step3_input.pdb")
    ap.add_argument("--prmtop", default="", help="chamber prmtop to cross-check numbering against")
    ap.add_argument("--summary", default="split_summary.txt")
    ap.add_argument("--check", metavar="DIR", help="verify component prmtops in DIR instead of printing vars")
    args = ap.parse_args()

    try:
        out, summary = derive(args.psf, args.pdb, args.prmtop)
    except SplitError as e:
        print(f"SPLIT ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    if args.check:
        bad = check(args.check, out)
        sys.exit(1 if bad else 0)

    with open(args.summary, "w") as f:
        f.write("\n".join(summary) + "\n")
    for line in summary:
        print("# " + line, file=sys.stderr)
    for k, v in out.items():
        print(f'{k}="{v}"')


if __name__ == "__main__":
    main()
