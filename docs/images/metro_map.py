#!/usr/bin/env python3
"""Pipeline overview (metro map).  python3 docs/images/metro_map.py
Writes auto16s_metro_map.svg (and .png if cairosvg is installed) next to this script."""
import os

W, H = 2400, 980
FONT = "'Helvetica Neue', Helvetica, Arial, 'DejaVu Sans', sans-serif"
C = {"qc": "#0072B2", "asv": "#009E73", "tax": "#E69F00", "div": "#CC79A7", "da": "#D55E00", "rep": "#6B7280"}
INK, MUTED, BAND = "#111827", "#6B7280", "#F3F4F6"
LW = 14
Y = 430          # main line
YR = 800         # report line
out = []
add = out.append


def text(x, y, s, size=15, weight="normal", anchor="middle", fill=INK, extra=""):
    add(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
        f'text-anchor="{anchor}" fill="{fill}" {extra}>{s}</text>')


def line(pts, color, width=LW, dash=None, opacity=1, r=36):
    """polyline with rounded corners (quadratic curve of radius ~r at every bend)"""
    import math
    d = f"M {pts[0][0]},{pts[0][1]}"
    for i in range(1, len(pts) - 1):
        (x0, y0), (x1, y1), (x2, y2) = pts[i - 1], pts[i], pts[i + 1]
        l1, l2 = math.hypot(x1 - x0, y1 - y0), math.hypot(x2 - x1, y2 - y1)
        rr = min(r, l1 / 2, l2 / 2)
        ax, ay = x1 - (x1 - x0) / l1 * rr, y1 - (y1 - y0) / l1 * rr
        bx, by = x1 + (x2 - x1) / l2 * rr, y1 + (y2 - y1) / l2 * rr
        d += f" L {ax:.1f},{ay:.1f} Q {x1},{y1} {bx:.1f},{by:.1f}"
    d += f" L {pts[-1][0]},{pts[-1][1]}"
    da = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" '
        f'stroke-linejoin="round" stroke-opacity="{opacity}"{da}/>')


def station(x, y, color, r=13):
    add(f'<circle cx="{x}" cy="{y}" r="{r + 5}" fill="white"/>')
    add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="white" stroke="{color}" stroke-width="6"/>')


def label(x, y, title, sub, where="below", anchor="middle"):
    if where == "below":
        text(x, y + 50, title, 18, "bold", anchor)
        for i, s in enumerate(sub):
            text(x, y + 73 + 20 * i, s, 14.5, "normal", anchor, MUTED)
    else:
        n = len(sub)
        text(x, y - 38 - 20 * n, title, 18, "bold", anchor)
        for i, s in enumerate(sub):
            text(x, y - 34 - 20 * (n - 1 - i), s, 14.5, "normal", anchor, MUTED)


def band(x0, x1, num, title, color):
    add(f'<rect x="{x0}" y="150" width="{x1 - x0}" height="{YR - 150 - 60}" rx="22" fill="{BAND}"/>')
    add(f'<rect x="{x0 + 22}" y="150" width="{x1 - x0 - 44}" height="5" rx="2.5" fill="{color}"/>')
    text(x0 + 22, 192, num, 15, "bold", "start", color)
    text(x0 + 42, 192, title, 15, "bold", "start", INK)


def doc_icon(x, y, color, tag):
    add(f'<path d="M {x-24},{y-32} h 34 l 14,14 v 50 h -48 z" fill="white" stroke="{color}" stroke-width="3.5" stroke-linejoin="round"/>')
    add(f'<path d="M {x+10},{y-32} v 14 h 14" fill="none" stroke="{color}" stroke-width="3.5" stroke-linejoin="round"/>')
    text(x, y + 9, tag, 11, "bold", "middle", color)


# ---------------- background, title ----------------
add(f'<rect width="{W}" height="{H}" fill="white"/>')
text(60, 80, "auto16s", 42, "bold", "start")
text(60, 118, "Raw reads to diversity statistics, differential abundance and a written report  ·  "
     "primers, truncation and rarefaction depth chosen from the data", 16, "normal", "start", MUTED)

# ---------------- coordinates ----------------
X = {"in": 120, "fastqc": 290, "pdet": 470, "ptrim": 650,
     "filt": 860, "err": 1040, "merge": 1220,
     "tax": 1430, "ps": 1600}
SPLIT = X["ps"] + 100
B = {"comp": 1930, "da": 1930, "rare": 1915, "alpha": 2050, "beta": 2210}
YB = {"comp": Y - 140, "da": Y, "rare": Y + 140}
XR1, XR2 = 1530, 2230        # report stations: QC summary, final report
XC = 2330                    # collector of the result lines

# stage bands (behind everything)
band(40, 740, "1", "Quality control, primers and adapters", C["qc"])
band(760, 1320, "2", "ASV inference (DADA2)", C["asv"])
band(1340, 1700, "3", "Taxonomy", C["tax"])
band(1720, 2370, "4", "Statistics and figures", C["div"])

# ---------------- report connectors (under the lines) ----------------
line([(X["fastqc"], Y), (X["fastqc"], YR), (XR2, YR)], C["rep"], 8, opacity=0.9)
line([(X["merge"], Y), (X["merge"], YR)], C["rep"], 3.5, "1 10")
YF = 600                     # FastQC after filtering (spur below the filter station)
line([(X["filt"], YF), (X["filt"], YR)], C["rep"], 3.5, "1 10")
line([(B["comp"] + 130, YB["comp"]), (XC, YB["comp"]), (XC, YR), (XR2, YR)], C["rep"], 3.5, "1 10")
line([(B["da"] + 130, YB["da"]), (XC, YB["da"])], C["rep"], 3.5, "1 10")
line([(B["beta"] + 40, YB["rare"]), (XC, YB["rare"])], C["rep"], 3.5, "1 10")

# ---------------- metro lines ----------------
line([(X["in"], Y), (X["ptrim"] + 100, Y)], C["qc"])
line([(X["ptrim"] + 100, Y), (X["merge"] + 110, Y)], C["asv"])
# 45-degree branches with rounded bends, drawn under the end of the main line
line([(SPLIT - 30, Y), (SPLIT + 20, Y), (SPLIT + 160, YB["rare"]), (B["beta"] + 40, YB["rare"])], C["div"], r=60)
line([(SPLIT - 30, Y), (SPLIT + 20, Y), (SPLIT + 160, YB["comp"]), (B["comp"] + 130, YB["comp"])], C["tax"], r=60)
line([(SPLIT - 30, Y), (B["da"] + 130, YB["da"])], C["da"])
line([(X["merge"] + 110, Y), (SPLIT, Y)], C["tax"])
line([(X["filt"], Y), (X["filt"], YF)], C["qc"], 9)

# ---------------- input ----------------
doc_icon(X["in"], Y, C["qc"], "fastq")
text(X["in"], Y + 66, "Raw reads", 18, "bold")
text(X["in"], Y + 88, "R1 / R2 + samplesheet", 14.5, "normal", "middle", MUTED)

# ---------------- main stations ----------------
main = [
    ("fastqc", C["qc"], "FastQC", ["raw reads"], "above"),
    ("pdet", C["qc"], "Primer detection", ["library of common", "16S primers"], "below"),
    ("ptrim", C["qc"], "Primer + adapter", ["removal, then checked"], "above"),
    ("filt", C["asv"], "Filter + truncate", ["truncLen keeping most reads,", "maxEE = 2"], "above"),
    ("err", C["asv"], "Denoise", ["error models R1 / R2"], "below"),
    ("merge", C["asv"], "Merge + chimeras", ["mergePairs,", "removeBimeraDenovo"], "above"),
    ("tax", C["tax"], "Taxonomy", ["SILVA 138.1"], "below"),
    ("ps", C["tax"], "phyloseq", ["drop chloroplasts,", "mitochondria"], "above"),
]
for k, col, t, sub, where in main:
    station(X[k], Y, col)
    label(X[k], Y, t, sub, where)
station(X["filt"], YF, C["qc"], 11)
text(X["filt"] + 28, YF - 2, "FastQC", 17, "bold", "start")
text(X["filt"] + 28, YF + 18, "filtered reads", 14.5, "normal", "start", MUTED)

# ---------------- branch stations ----------------
station(B["comp"], YB["comp"], C["tax"])
text(B["comp"], YB["comp"] - 52, "Composition", 18, "bold")
text(B["comp"], YB["comp"] - 30, "bar plots, heatmap", 14.5, "normal", "middle", MUTED)
station(B["da"], YB["da"], C["da"])
text(B["da"], YB["da"] - 52, "Differential abundance", 18, "bold")
text(B["da"], YB["da"] - 30, "MaAsLin 3: abundance + prevalence", 14.5, "normal", "middle", MUTED)
station(B["rare"], YB["rare"], C["div"])
label(B["rare"], YB["rare"], "Rarefaction", ["automatic depth"], "below")
station(B["alpha"], YB["rare"], C["div"])
label(B["alpha"], YB["rare"], "Alpha", ["Shannon, Simpson"], "below")
station(B["beta"], YB["rare"], C["div"])
label(B["beta"], YB["rare"], "Beta", ["PCoA, PERMANOVA"], "below")

# ---------------- reports ----------------
add(f'<rect x="{XR1 - 34}" y="{YR - 24}" width="68" height="48" rx="12" fill="white" stroke="{C["rep"]}" stroke-width="5"/>')
text(XR1, YR + 6, "QC", 15, "bold", "middle", C["rep"])
text(XR1, YR + 54, "QC summary (MultiQC)", 17, "bold")
text(XR1, YR + 75, "FastQC, primers, DADA2 steps", 14.5, "normal", "middle", MUTED)
add(f'<circle cx="{XR2}" cy="{YR}" r="44" fill="white"/>')
doc_icon(XR2, YR, INK, "HTML")
text(XR2 - 56, YR + 54, "Final report", 19, "bold", "end")
text(XR2 - 56, YR + 76, "samples · methods with the values used · results", 14.5, "normal", "end", MUTED)

# ---------------- legend ----------------
lx, ly = 60, 940
items = [("qc", "Read QC, primers"), ("asv", "ASV inference"), ("tax", "Taxonomy, composition"),
         ("da", "Differential abundance"), ("div", "Diversity (rarefied)"), ("rep", "Reports")]
x = lx
for k, s in items:
    line([(x, ly), (x + 40, ly)], C[k], 10)
    text(x + 54, ly + 5, s, 15, "normal", "start", "#374151")
    x += 54 + 8.6 * len(s) + 50

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
