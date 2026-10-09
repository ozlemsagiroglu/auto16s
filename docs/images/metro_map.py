#!/usr/bin/env python3
"""Pipeline overview (metro map).  python3 docs/images/metro_map.py
Writes auto16s_metro_map.svg (and .png if cairosvg is installed) next to this script."""
import os

W, H = 2360, 800
FONT = "Helvetica, Arial, 'DejaVu Sans', sans-serif"
C = {"qc": "#0072B2", "asv": "#009E73", "tax": "#E69F00", "div": "#CC79A7", "da": "#D55E00", "rep": "#8C8C8C"}
LW = 12
Y = 340          # main line
YR = 640         # report line
out = []
add = out.append


def text(x, y, s, size=15, weight="normal", anchor="middle", fill="#222", extra=""):
    add(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
        f'text-anchor="{anchor}" fill="{fill}" {extra}>{s}</text>')


def line(pts, color, width=LW, dash=None):
    d = "M " + " L ".join(f"{x},{y}" for x, y in pts)
    da = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" '
        f'stroke-linejoin="round"{da}/>')


def station(x, y, r=14):
    add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="white" stroke="#1a1a1a" stroke-width="4"/>')


def label(x, y, title, sub, where="below"):
    """title + grey sub-lines, fully clear of the line"""
    if where == "below":
        text(x, y + 46, title, 17, "bold")
        for i, s in enumerate(sub):
            text(x, y + 68 + 19 * i, s, 14, "normal", "middle", "#555")
    else:
        n = len(sub)
        text(x, y - 34 - 19 * n, title, 17, "bold")
        for i, s in enumerate(sub):
            text(x, y - 30 - 19 * (n - 1 - i), s, 14, "normal", "middle", "#555")


# ---------------- background, title ----------------
add(f'<rect width="{W}" height="{H}" fill="white"/>')
text(60, 72, "auto16s", 36, "bold", "start", "#111")
text(60, 104, "Paired-end 16S rRNA amplicon reads  →  ASVs  →  taxonomy  →  diversity and differential abundance  →  one report",
     17, "normal", "start", "#555")

# ---------------- coordinates ----------------
X = {"in": 110, "fastqc": 280, "pdet": 460, "ptrim": 640, "filt": 840, "err": 1020, "merge": 1200,
     "tax": 1390, "ps": 1570}
SPLIT = X["ps"] + 110            # branches leave the main line here
BX = SPLIT + 120                 # ... and are horizontal from here
YB = {"comp": 178, "da": 268, "rare": 470}
B = {"comp": 1880, "da": 1880, "rare": 1850, "alpha": 2010, "beta": 2170}
XR = 2290                        # report station

# DADA2 box (behind everything)
add(f'<rect x="{X["filt"]-95}" y="{Y-150}" width="{X["merge"]-X["filt"]+190}" height="290" rx="20" '
    f'fill="{C["asv"]}" fill-opacity="0.06" stroke="{C["asv"]}" stroke-width="2" stroke-dasharray="7 6"/>')
text(X["filt"] - 75, Y - 120, "DADA2", 18, "bold", "start", C["asv"])

# ---------------- report connectors (drawn first, under the lines) ----------------
line([(X["fastqc"], Y), (X["fastqc"], YR), (XR, YR)], C["rep"], 8)
line([(X["err"], Y), (X["err"], YR)], C["rep"], 4, "1 10")
YF = 530                         # FastQC after filtering
line([(X["filt"], YF), (X["filt"], YR)], C["rep"], 4, "1 10")
line([(B["comp"] + 120, YB["comp"]), (XR, YB["comp"]), (XR, YR)], C["rep"], 4, "1 10")
line([(B["da"] + 120, YB["da"]), (XR, YB["da"])], C["rep"], 4, "1 10")
line([(B["beta"], YB["rare"]), (XR, YB["rare"])], C["rep"], 4, "1 10")

# ---------------- metro lines ----------------
line([(X["in"], Y), (X["ptrim"] + 90, Y)], C["qc"])
line([(X["ptrim"] + 90, Y), (X["merge"] + 95, Y)], C["asv"])
line([(X["merge"] + 95, Y), (SPLIT, Y)], C["tax"])
line([(X["filt"], Y), (X["filt"], YF)], C["qc"], 8)
line([(SPLIT, Y), (BX, YB["comp"]), (B["comp"] + 120, YB["comp"])], C["tax"])
line([(SPLIT, Y), (BX, YB["da"]), (B["da"] + 120, YB["da"])], C["da"])
line([(SPLIT, Y), (BX, YB["rare"]), (B["beta"], YB["rare"])], C["div"])

# ---------------- input ----------------
x = X["in"]
add(f'<path d="M {x-19},{Y-26} h 27 l 11,11 v 41 h -38 z" fill="white" stroke="{C["qc"]}" stroke-width="3"/>')
add(f'<path d="M {x+8},{Y-26} v 11 h 11" fill="none" stroke="{C["qc"]}" stroke-width="3"/>')
text(x, Y + 7, "fastq", 11, "bold", "middle", C["qc"])
label(x, Y, "Raw reads", ["R1 / R2", "+ samplesheet"], "below")

# ---------------- main stations ----------------
main = [
    ("fastqc", "FastQC", ["raw reads"], "above"),
    ("pdet", "Primer detection", ["library of common", "16S primers"], "below"),
    ("ptrim", "Primer removal", ["per read, then checked;", "spacers, reverse pairs"], "above"),
    ("filt", "Filter + truncate", ["truncLen (automatic)", "maxEE = 2"], "above"),
    ("err", "Error model + denoise", ["learnErrors, dada", "R1 and R2 separately"], "above"),
    ("merge", "Merge + chimeras", ["mergePairs", "removeBimeraDenovo"], "below"),
    ("tax", "Taxonomy", ["SILVA 138.1", "assignTaxonomy"], "above"),
    ("ps", "phyloseq", ["drop Eukaryota, chloroplasts,", "mitochondria"], "below"),
]
for k, t, sub, where in main:
    station(X[k], Y)
    label(X[k], Y, t, sub, where)

# ---------------- branch stations ----------------
station(X["filt"], YF)
text(X["filt"] + 26, YF - 3, "FastQC", 17, "bold", "start")
text(X["filt"] + 26, YF + 17, "filtered reads", 14, "normal", "start", "#555")
station(B["comp"], YB["comp"]); label(B["comp"], YB["comp"], "Composition", ["bar plots (phylum, family, genus),", "genus heatmap (all reads)"], "above")
station(B["da"], YB["da"]); label(B["da"], YB["da"], "Differential abundance", ["MaAsLin 3: abundance + prevalence,", "genus level (all reads)"], "below")
station(B["rare"], YB["rare"]); label(B["rare"], YB["rare"], "Rarefaction", ["automatic depth"], "below")
station(B["alpha"], YB["rare"]); label(B["alpha"], YB["rare"], "Alpha diversity", ["Observed, Shannon,", "Simpson + tests"], "below")
station(B["beta"], YB["rare"]); label(B["beta"], YB["rare"], "Beta diversity", ["Bray-Curtis PCoA,", "PERMANOVA, betadisper"], "below")


# ---------------- report ----------------
add(f'<rect x="{XR-38}" y="{YR-26}" width="76" height="52" rx="11" fill="white" stroke="#1a1a1a" stroke-width="4"/>')
text(XR, YR + 6, "HTML", 14, "bold")
text(XR - 50, YR + 58, "MultiQC report", 17, "bold", "end")
text(XR - 50, YR + 79, "QC, warnings, all figures and statistics", 14, "normal", "end", "#555")

# ---------------- legend ----------------
lx, ly = 60, 755
items = [("qc", "Read QC, primers"), ("asv", "ASV inference"), ("tax", "Taxonomy, phyloseq, composition"),
         ("da", "Differential abundance"), ("div", "Diversity (rarefied)"), ("rep", "Report")]
x = lx
for k, s in items:
    line([(x, ly), (x + 42, ly)], C[k], 10)
    text(x + 56, ly + 5, s, 15, "normal", "start", "#333")
    x += 56 + 9.2 * len(s) + 46

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
       + "\n".join(out) + "\n</svg>\n")
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "auto16s_metro_map.svg"), "w") as f:
    f.write(svg)
try:
    import cairosvg
    cairosvg.svg2png(bytestring=svg.encode(), write_to=os.path.join(here, "auto16s_metro_map.png"),
                     output_width=int(W * 1.25))
except ImportError:
    pass
