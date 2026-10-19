#!/usr/bin/env python3
"""
build_docx.py -- assemble the rendered figures and their legends into one Word document.

Run the figure scripts first (fig4.py, fig5.py, fig6.py, si_a.py, si_b.py), then this.
Writes manuscript/Surf_for_Turf_Figures.docx. Figures 1-3 and Supp. Figs 1-5 are the
experimental panels and are listed as placeholders, so the numbering matches the manuscript.
"""
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

OUT = Path(__file__).resolve().parent / "out"
DOC = Path(__file__).resolve().parent.parent / "manuscript" / "Surf_for_Turf_Figures.docx"
INK = RGBColor(0x03, 0x13, 0x26)
GREY = RGBColor(0x65, 0x5e, 0x77)
WIDE, NARROW = 6.6, 5.2          # inches on a Letter page with 1" margins

MAIN = [
    ("Figure 1", None, "Short-term iBSC growth with FGF2 orthologs and heparin co-factors. "
     "a, Nuclei counts after 5 days in Beefy-9 containing 40 ng/mL of each FGF2 variant, versus a "
     "no-FGF2 control. b, Effect of 2 µg/mL heparin, 100 ng/mL heparinase I, or both, on the "
     "human and goldfish FGF2 conditions."),
    ("Figure 2", None, "Multi-passage growth curves. Cumulative population doublings for iBSCs in "
     "Beefy-9 with human or goldfish FGF2, with or without heparinase, across three trials starting "
     "at passages 32, 45 and 79."),
    ("Figure 3", None, "Myogenicity of iBSCs in the different Beefy-9 media. PAX7 staining of "
     "proliferating cells and myosin heavy chain staining of differentiated cells."),
    ("Figure 4", "fig4", "Molecular dynamics of FGF2–FGFR–heparin assemblies across orthologs. "
     "a, The two simulated assemblies, shown as Cα traces of the rebuilt goldfish systems: the "
     "2:2:2 complex with two heparin chains (PDB 1FQ9) and the 2:2:1 complex with one shared heparin "
     "chain (PDB 1E0O). b,c, MM/GBSA whole-assembly binding energies for seven FGF2 orthologs in the "
     "2:2:2 and 2:2:1 assemblies. Points are independent replicates; bars give the mean ± SD; "
     "one-way ANOVA p = 0.92 and 0.21. More negative means stronger binding. d, Positions differing "
     "from human FGF2 in each ortholog, with FGFR-contacting and heparin-contacting positions marked "
     "beneath (contacts within 4.5 Å in the 1FQ9 model). 28 of 40 FGFR-contacting and 11 of 12 "
     "heparin-contacting positions are identical in all seven orthologs, against 65% identity across "
     "the whole core. The 86–97 segment swapped in Fig. 6 and the conserved Arg90 are marked."),
    ("Figure 5", "fig5", "Extended simulations of human and goldfish FGF2. a, Whole-assembly binding "
     "energy averaged over the 11–51 ns window for both assemblies; points are replicates, bars "
     "give the mean ± SD, and Δ is the goldfish − human difference with its 95% "
     "confidence interval (Welch t-test on replicate means, n = 3 per group). This is the "
     "pre-specified test of the v1 claim. b, Heparin-site salt bridges over time in the 2:2:1 "
     "assembly (1 ns running mean); shading marks the first 10 ns, excluded as equilibration, and "
     "the dash-dotted line is the tight-state threshold fixed in advance. c,d, Arg90–heparin "
     "heavy-atom contacts and Arg90 per-residue energy contribution, for each FGF2 copy and "
     "replicate; the shaded band spans the gap between the parent orthologs and the dashed line is "
     "the engagement threshold. e, Per-residue goldfish − human difference along the sequence in "
     "the 2:2:1 assembly, summed over both FGF2 copies, with substituted positions coloured. Panels "
     "c–e are exploratory and were not part of the pre-specified plan."),
    ("Figure 6", "fig6", "Segment-swap chimeras testing the Arg90 binding mode. a, Design. FGF2 "
     "86–97 carries eight of the 25 human-to-goldfish substitutions (K86N, E87S, L92I, A93G, "
     "S94T, K95R, C96R, V97T) and flanks the conserved Arg90. Each chimera keeps its host backbone "
     "throughout and takes the donor side chains across this segment. b, Arg90–heparin contacts "
     "in the 2:2:1 assembly for the two parent orthologs, against the engagement threshold fixed "
     "before the chimera runs. The pre-registered rules are that the goldfish segment is sufficient "
     "if the human + goldfish 86–97 chimera engages in at least two of three replicates, and "
     "necessary if the goldfish + human 86–97 chimera engages in at most one of three. Chimera "
     "simulations are in progress."),
]

SI = [
    ("Supplementary Figure 1", None, "Heparin and heparinase across all FGF2 orthologs. Short-term "
     "growth for every ortholog with no addition, heparin, heparinase, or both."),
    ("Supplementary Figure 2", None, "Cumulative doublings across the three long-term trials."),
    ("Supplementary Figure 3", None, "Additional immunostaining images."),
    ("Supplementary Figure 4", None, "Quantification of PAX7+ and MHC+ nuclei."),
    ("Supplementary Figure 5", None, "Western blots and quantification. a, Semi-quantitation for "
     "ERK, p-ERK, FGFR and p-FGFR. b, Passage 47 blots. c, Passage 34 blots."),
    ("Supplementary Figure 6", "figS6", "FGF2 core definition. Residues retained from each AlphaFold "
     "model after trimming to the template-resolved core, per ortholog and template. The trimmed "
     "N-terminal segment is unresolved in every FGF1/FGF2 crystal structure and, when retained, "
     "threaded through the receptor in several builds."),
    ("Supplementary Figure 7", "figS7", "Build quality control. Largest single-residue van der Waals "
     "contribution per replicate, for the original CHARMM-GUI builds (open symbols) and the "
     "consistent rebuild (filled). Values of hundreds of kcal/mol in the original builds indicate "
     "FGF2 N-terminal tails threaded through FGFR rings; every rebuilt system is at or below "
     "2.4 kcal/mol. Note the logarithmic scale."),
    ("Supplementary Figure 8", "figS8", "Per-interface binding energies for all orthologs. Interface "
     "1 and interface 2 in both assemblies, with one-way ANOVA p-values. No pairwise comparison is "
     "significant under Tukey's HSD."),
    ("Supplementary Figure 9", "figS9", "Convergence of the binding energies. a,b, Whole-assembly "
     "energy at 1 ns versus over the 11–51 ns window, for each extended replicate. c, Within-"
     "replicate uncertainty (mean SEM of four 10 ns blocks) against between-replicate SD. The "
     "between-replicate term dominates, so additional replicates help more than longer trajectories."),
    ("Supplementary Figure 10", "figS10", "Full extension time series. Heparin-site salt bridges for "
     "both FGF2 copies and the FGFR–FGFR distance over 51 ns, for both assemblies (1 ns running "
     "mean; shading marks the excluded first 10 ns). Line style distinguishes replicates."),
    ("Supplementary Figure 11", "figS11", "Component binding energies. FGF2–heparin (FGFR removed) "
     "and FGF2–FGFR (heparin removed) energies for human and goldfish in both assemblies. "
     "Exploratory; p-values are uncorrected."),
    ("Supplementary Figure 12", "figS12", "Per-residue decomposition in the 2:2:2 assembly, the "
     "companion to Fig. 5e. Arg90 does not stand out when each FGF2 has its own heparin."),
    ("Supplementary Figure 13", "figS13", "Static electrostatics of the FGF2 cores. a, Net charge of "
     "each ortholog core. b, Screened Coulomb (Debye–Hückel, ε = 80, 0.15 M) interaction "
     "with heparin in the starting structures. Goldfish and human have the same net charge and "
     "near-identical electrostatic attraction, so the Arg90 difference is not a simple charge effect."),
    ("Supplementary Figure 14", "figS14", "Statistical power of the MM/GBSA comparison. Smallest "
     "goldfish–human difference detectable at 80% power against the number of replicates per "
     "species, given the observed between-replicate spread. At n = 3 the comparison resolves about "
     "76 kcal/mol, far above the roughly 1 kcal/mol implied by the growth data."),
]


def add_heading(doc, text, size=13):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(size)
    r.font.color.rgb = INK
    p.space_after = Pt(6)
    return p


def add_entry(doc, label, stem, legend, width):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    png = OUT / f"{stem}.png" if stem else None
    if png and png.exists():
        p.add_run().add_picture(str(png), width=Inches(width))
    else:
        r = p.add_run("[experimental figure — to be supplied by the Tufts group]")
        r.italic = True
        r.font.size = Pt(9)
        r.font.color.rgb = GREY
    cap = doc.add_paragraph()
    cr = cap.add_run(f"{label} | ")
    cr.bold = True
    cr.font.size = Pt(9)
    cr.font.color.rgb = INK
    tr = cap.add_run(legend)
    tr.font.size = Pt(9)
    tr.font.color.rgb = INK
    cap.space_after = Pt(18)


def main():
    doc = Document()
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(0.9)
        s.top_margin = s.bottom_margin = Inches(0.9)
    st = doc.styles["Normal"]
    st.font.name = "Helvetica Neue"
    st.font.size = Pt(10)

    add_heading(doc, "Surf for Turf — figures", 16)
    p = doc.add_paragraph()
    r = p.add_run("Simulation figures rendered from the repository data; colours are cmcrameri "
                  "lipari throughout. Figures 1–3 and Supplementary Figures 1–5 are the "
                  "experimental panels and are left as placeholders so the numbering matches the "
                  "manuscript. Figure 6b fills in when the chimera simulations finish.")
    r.font.size = Pt(9)
    r.font.color.rgb = GREY
    p.space_after = Pt(16)

    for label, stem, legend in MAIN:
        add_entry(doc, label, stem, legend, WIDE if stem else NARROW)

    doc.add_page_break()
    add_heading(doc, "Supplementary Information", 14)
    for label, stem, legend in SI:
        add_entry(doc, label, stem, legend, WIDE if stem else NARROW)

    DOC.parent.mkdir(exist_ok=True)
    doc.save(DOC)
    print("wrote", DOC)


if __name__ == "__main__":
    main()
