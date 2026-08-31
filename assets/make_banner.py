"""Generate the repository banner as a self-contained SVG.

No <style> blocks, no scripts, no prefers-color-scheme: GitHub's sanitiser
strips the first two, and the third follows the browser rather than the
GitHub theme. Every colour is a presentation attribute.
"""
import math
from pathlib import Path

W, H = 1280, 320
SLATE = "#12233b"
SLATE_EDGE = "#24405f"
BLUE = "#4B8BBE"
BLUE_DEEP = "#3776AB"
YELLOW = "#FFD43B"
WHITE = "#ffffff"
MUTED = "#b9c6d6"
DIM = "#8d9fb5"

# ---------------------------------------------------------------- helix
CX, AMP = 168.0, 62.0
TOP, BOT = 58.0, 268.0
PERIODS = 2.0
STEPS = 160


def strand(phase):
    pts = []
    for i in range(STEPS + 1):
        f = i / STEPS
        y = TOP + (BOT - TOP) * f
        x = CX + AMP * math.sin(2 * math.pi * PERIODS * f + phase)
        pts.append((x, y))
    return pts


def path_d(pts):
    head = f"M {pts[0][0]:.1f} {pts[0][1]:.1f}"
    return head + "".join(f" L {x:.1f} {y:.1f}" for x, y in pts[1:])


front = strand(0.0)             # blue strand, carries the snake head
back = strand(math.pi)          # yellow strand

parts = []
add = parts.append

add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
    f'width="{W}" height="{H}" role="img" '
    f'aria-label="Python for Bioinformatics - learn by running the code">')

# panel
add(f'<rect x="0" y="0" width="{W}" height="{H}" rx="26" fill="{SLATE}"/>')
add(f'<rect x="1.5" y="1.5" width="{W-3}" height="{H-3}" rx="25" '
    f'fill="none" stroke="{SLATE_EDGE}" stroke-width="3"/>')

# faint base-letter texture behind the helix, so it reads as sequence
letters = "ACGT" * 5
for i, ch in enumerate(letters[:17]):
    add(f'<text x="{330 + i * 56}" y="46" font-family="Verdana, DejaVu Sans, sans-serif" '
        f'font-size="15" fill="{SLATE_EDGE}" opacity="0.5">{ch}</text>')

# base-pair rungs, drawn under the strands
for i in range(10, STEPS - 8, 21):
    x1, y1 = front[i]
    x2, y2 = back[i]
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="#7d90a8" stroke-width="2" stroke-linecap="round" opacity="0.45"/>')

# strands
add(f'<path d="{path_d(back)}" fill="none" stroke="{YELLOW}" '
    f'stroke-width="9" stroke-linecap="round"/>')
add(f'<path d="{path_d(front)}" fill="none" stroke="{BLUE}" '
    f'stroke-width="9" stroke-linecap="round"/>')

# snake head on top of the blue strand, angled along its tangent
hx, hy = front[0]
x2, y2 = front[3]
angle = math.degrees(math.atan2(y2 - hy, x2 - hx)) + 180
add(f'<g transform="translate({hx:.1f} {hy:.1f}) rotate({angle:.1f})">')
add(f'<ellipse cx="10" cy="0" rx="19" ry="11.5" fill="{BLUE}"/>')
add(f'<circle cx="17" cy="-4" r="3.2" fill="{SLATE}"/>')
add(f'<path d="M 29 0 L 38 0 M 38 0 L 44 -3.5 M 38 0 L 44 3.5" '
    f'stroke="{YELLOW}" stroke-width="2" fill="none" stroke-linecap="round"/>')
add('</g>')

# ---------------------------------------------------------------- text
TX = 330
FONT = "Verdana, DejaVu Sans, Geneva, sans-serif"

add(f'<text x="{TX}" y="132" font-family="{FONT}" font-size="50" '
    f'font-weight="bold" fill="{WHITE}" letter-spacing="1.5">PYTHON FOR</text>')
add(f'<text x="{TX}" y="186" font-family="{FONT}" font-size="50" '
    f'font-weight="bold" fill="{YELLOW}" letter-spacing="1.5">BIOINFORMATICS</text>')
add(f'<rect x="{TX}" y="206" width="96" height="5" rx="2.5" fill="{BLUE}"/>')
add(f'<text x="{TX}" y="243" font-family="{FONT}" font-size="21" '
    f'fill="{MUTED}">Learn by running the code, not by reading it</text>')

# ---------------------------------------------------------------- badges
BY = 286          # badge baseline


def badge(x, label, icon):
    icon(x)
    add(f'<text x="{x + 38}" y="{BY + 6}" font-family="{FONT}" font-size="18" '
        f'fill="{DIM}">{label}</text>')


def cell_icon(x):
    add(f'<circle cx="{x + 11}" cy="{BY}" r="11" fill="none" '
        f'stroke="{BLUE}" stroke-width="2.6"/>')
    add(f'<circle cx="{x + 11}" cy="{BY}" r="4.2" fill="{YELLOW}"/>')


def hex_icon(x):
    cx, cy, r = x + 11, BY, 11.5
    pts = " ".join(
        f"{cx + r * math.cos(math.radians(a)):.1f},{cy + r * math.sin(math.radians(a)):.1f}"
        for a in range(-90, 270, 60))
    add(f'<polygon points="{pts}" fill="none" stroke="{BLUE}" stroke-width="2.6"/>')


def cap_icon(x):
    cx, cy = x + 12, BY
    add(f'<polygon points="{cx},{cy-9} {cx+15},{cy-2} {cx},{cy+5} {cx-15},{cy-2}" '
        f'fill="{YELLOW}"/>')
    add(f'<path d="M {cx+8} {cy+1} L {cx+8} {cy+9}" stroke="{BLUE}" '
        f'stroke-width="2.4" stroke-linecap="round"/>')


badge(TX, "28 runnable lessons", cell_icon)
badge(TX + 275, "Biopython + scikit-bio", hex_icon)
badge(TX + 590, "MIT licensed", cap_icon)

add('</svg>')

out = Path("assets/banner.svg")
out.parent.mkdir(exist_ok=True)
out.write_text("\n".join(parts) + "\n", encoding="utf-8")
print("wrote", out, out.stat().st_size, "bytes")
