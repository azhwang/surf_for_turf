#!/usr/bin/env python3
"""
render_structures.py -- cartoon renders of the two assemblies for Fig. 4a, drawn in PyMOL.

Run with the PyMOL environment's interpreter, not the main one:
    ~/anaconda3/envs/pymolenv/bin/python figures/render_structures.py

Reads each template's rebuilt goldfish starting structure, colours FGFR / FGF2 / heparin in the
same lipari slots as the rest of the figure, and ray-traces a transparent PNG per assembly into
figures/out/. Heparin is drawn as sticks so the chain reads against the protein cartoons.
"""
from pathlib import Path

from PIL import Image
from pymol import cmd

ROOT = Path(__file__).resolve().parent.parent
RB = ROOT / "rebuild_v2"
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

COL = {"fgfrA": "#4c5a7a", "fgfrB": "#8c606d", "fgfC": "#cb685f", "fgfD": "#e5ad82",
       "hep": "#031326"}
W, H = 1500, 1250


def crop(path, pad=8):
    """Trim the transparent margin so the render fills its panel."""
    im = Image.open(path).convert("RGBA")
    bbox = im.split()[3].getbbox()
    if bbox:
        x0, y0, x1, y1 = bbox
        im = im.crop((max(0, x0 - pad), max(0, y0 - pad),
                      min(im.width, x1 + pad), min(im.height, y1 + pad)))
        im.save(path)


def normalise(src, dst):
    """Right-justify one-letter atom names into column 14 and write an explicit element column.

    The rebuild writes FGF2 atom names left-justified from column 13, so a PDB reader takes
    "CA  " as calcium rather than a backbone carbon, the chain is never recognised as a polymer
    and no cartoon is built. The scaffold lines are already standard and pass through unchanged.
    """
    out = []
    for l in open(src):
        if l.startswith(("ATOM", "HETATM")):
            nm = l[12:16].strip()
            el = nm[0] if nm[:2] not in ("CL", "NA", "FE", "ZN", "MG") else nm[:2]
            field = f" {nm:<3s}" if len(nm) < 4 and len(el) == 1 else f"{nm:<4s}"
            l = f"{l[:12]}{field}{l[16:72]}{l[72:76]:<4s}  {el:>2s}\n"
        out.append(l)
    Path(dst).write_text("".join(out))


def render(tk, name):
    cmd.reinitialize()
    src = RB / tk / "goldfish" / "precheck.pdb"
    tmp = OUT / f"_norm_{name}.pdb"
    normalise(src, tmp)
    cmd.load(str(tmp), "asm")
    cmd.remove("hydrogens")
    # every FGF2 copy is written as chain X, so PyMOL cannot build separate cartoons;
    # give each segment its own chain, then let PyMOL assign secondary structure
    for seg, ch in (("PROA", "A"), ("PROB", "B"), ("PROC", "C"), ("PROD", "D")):
        cmd.alter(f"segi {seg}", f"chain='{ch}'")
    cmd.alter("segi CAR*", "chain='H'")
    cmd.sort()
    cmd.rebuild()
    cmd.dss()
    cmd.hide("everything")

    for key, hexcol in COL.items():
        cmd.set_color(f"c_{key}", [int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5)])

    sel = {"fgfrA": "segi PROA", "fgfrB": "segi PROB", "fgfC": "segi PROC", "fgfD": "segi PROD",
           "hep": "segi CAR*"}
    for key, s in sel.items():
        cmd.select(key, f"asm and {s}")

    cmd.show("cartoon", "fgfrA or fgfrB or fgfC or fgfD")
    cmd.show("sticks", "hep")
    cmd.set("stick_radius", 0.30, "hep")
    for key in sel:
        cmd.color(f"c_{key}", key)

    cmd.set("cartoon_transparency", 0.12, "fgfrA or fgfrB")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("ambient", 0.28)
    cmd.set("specular", 0.12)
    cmd.set("direct", 0.55)
    cmd.set("reflect", 0.25)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("orthoscopic", 1)

    cmd.orient("asm")
    cmd.zoom("asm", 1.0)
    cmd.ray(W, H)
    path = OUT / f"struct_{name}.png"
    cmd.png(str(path), dpi=400)
    crop(path)
    print("wrote", path)


if __name__ == "__main__":
    render("1fq9", "222")
    render("1e0o", "221")
