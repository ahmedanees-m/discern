"""Shared figure style: the journal's constraints in one place.

Human Genomics specifies the geometry, so it is expressed here as constants rather than repeated
as magic numbers in each figure:

  - 85 mm for a half-page figure, 170 mm for a full-page figure
  - at most 225 mm for figure plus legend
  - about 300 dpi at final size
  - no line thinner than 0.25 pt
  - all fonts embedded

Colour follows the Okabe-Ito palette, which stays distinguishable under the common forms of colour
vision deficiency and in greyscale. Nothing in these figures is encoded by colour alone: every
series is also separated by position, and where two things must be told apart in a dense panel they
differ in shape or hatch as well.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MM = 1 / 25.4
FULL_W = 170 * MM          # full page width, inches
HALF_W = 85 * MM           # half page width, inches
MAX_H = 225 * MM           # figure plus legend

# Okabe-Ito, plus a neutral grey for de-emphasis.
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERMILION = "#D55E00"
SKY = "#56B4E9"
PURPLE = "#CC79A7"
YELLOW = "#F0E442"
GREY = "#6E6E6E"
LIGHT = "#D9D9D9"

# Semantic roles, so a change of palette does not mean hunting through every figure.
DISCERN = BLUE
REVEL = ORANGE
ALPHAMISSENSE = GREEN
EXHIBIT = VERMILION        # GeneBe, shown as a circularity exhibit
EXPERT = "#4D4D4D"         # the expert panels
ROUTED = GREEN             # routed away by the partition
NO_INPUT = ORANGE          # variant-intrinsic, no input available
NOT_IMPL = GREY            # not implemented

# Type sizes. Small, because a full-width figure is reproduced at 170 mm. MICRO is reserved for
# in-panel annotation blocks, which sit beside the data and read better a little under the label
# size; axis labels and tick labels never go below SMALL.
BASE, SMALL, TINY = 8.0, 7.0, 6.2
MICRO = 5.8
TITLE = 9.0

# Annotation text on a light background. GREY is right for rules and de-emphasised marks but too
# pale for words a reader has to read.
INK = "#333333"

RC = {
    "figure.dpi": 150,
    "savefig.dpi": 400,                 # comfortably above the 300 dpi target
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,         # crop tightly, as the guidelines ask
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": BASE,
    "axes.titlesize": TITLE,
    "axes.labelsize": BASE,
    "legend.fontsize": SMALL,
    "xtick.labelsize": SMALL,
    "ytick.labelsize": SMALL,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,              # well above the 0.25 pt floor
    "grid.linewidth": 0.4,
    "lines.linewidth": 1.3,
    "patch.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,                 # embed as TrueType, not subsetted Type 3
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def apply():
    plt.rcParams.update(RC)


def panel_label(ax, letter, dx=-0.06, dy=1.06, size=TITLE):
    """Put a bold panel letter outside the axes, where it cannot collide with the data."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=size, fontweight="bold",
            va="top", ha="left")


def wrap(text, width):
    """Wrap to a character count, for labels that would otherwise run past the axes."""
    import textwrap
    return "\n".join(textwrap.wrap(text, width=width))


def fit(fig):
    """Shrink a figure until its saved extent is inside the page limits.

    Figures are saved on a tight bounding box, so the file is as wide as the drawn content, not as
    wide as ``figsize``: a tick label hanging past the left spine makes the saved PDF wider than
    the figure was declared to be. Figure 1 came out at 173 mm from a 170 mm canvas that way. The
    check therefore has to run on the tight extent, and it lives here so no figure can bypass it.
    """
    pad = 2 * RC["savefig.pad_inches"]
    total = 1.0
    # Text is sized in points and does not shrink with the canvas, so one pass undershoots:
    # scaling a 173 mm extent by 170/173 lands at about 170.2 mm. Iterate until it is inside.
    for _ in range(6):
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        w, h = bb.width + pad, bb.height + pad
        scale = min(FULL_W / w if w > FULL_W else 1.0, MAX_H / h if h > MAX_H else 1.0)
        if scale >= 1.0:
            break
        fw, fh = fig.get_size_inches()
        fig.set_size_inches(fw * scale, fh * scale)
        total *= scale
    return total


def flatten_tiff(path, dpi):
    """Rewrite a TIFF as opaque RGB.

    matplotlib writes RGBA. The alpha channel is fully opaque, but production workflows treat an
    alpha channel as a transparency to resolve, and resolve it against whatever the page is; RGB
    composited over white removes the question and is smaller.
    """
    from PIL import Image
    with Image.open(path) as im:
        if im.mode == "RGB":
            return
        rgba = im.convert("RGBA")
        flat = Image.new("RGB", rgba.size, (255, 255, 255))
        flat.paste(rgba, mask=rgba.split()[3])
    flat.save(path, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))


def save(fig, outdir, name):
    """Write the vector PDF the journal prefers, plus a high-resolution PNG and TIFF.

    TIFF is LZW-compressed: production systems ask for it, it is lossless, and an uncompressed
    400 dpi full-page figure would be tens of megabytes.
    """
    import os
    os.makedirs(outdir, exist_ok=True)
    fit(fig)
    dpi = RC["savefig.dpi"]
    paths = []
    for ext in ("pdf", "png", "tif"):
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p, **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}))
        if ext == "tif":
            flatten_tiff(p, dpi)
        paths.append(p)
    plt.close(fig)
    return paths
