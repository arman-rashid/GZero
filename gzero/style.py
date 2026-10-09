"""Plot styles.

Screen: normal font sizes. Export: Nature figure rules, one plot per file,
89 mm wide (120 or 183 mm if needed), Arial 7 pt, lines >= 0.5 pt, no
titles, pdf with editable text plus a 600 dpi png.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from dataclasses import dataclass, asdict, field

MM = 1 / 25.4  # inches per mm
WIDTHS_MM = {"single column (89 mm)": 89.0, "1.5 columns (120 mm)": 120.0, "double column (183 mm)": 183.0}

# Okabe-Ito colour-blind-safe palette (recommended for Nature figures)
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]


@dataclass
class FigureSettings:
    width: str = "single column (89 mm)"
    font: str = "Arial"
    font_size: float = 7.0          # axis labels; tick labels and legends are 1 pt smaller
    line_width: float = 0.75        # data lines (pt)
    axes_width: float = 0.5         # axes and ticks (pt)
    tick_direction: str = "out"
    formats: str = "pdf, png"        # comma-separated: pdf, png, tiff, svg, eps
    dpi: int = 600                  # raster outputs

    def format_list(self) -> list[str]:
        out = [f.strip().lower().lstrip(".") for f in str(self.formats).replace(";", ",").split(",")]
        return [f for f in out if f in ("pdf", "png", "tiff", "tif", "svg", "eps", "jpg")] or ["pdf"]

    def width_in(self) -> float:
        return WIDTHS_MM.get(self.width, 89.0) * MM

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Panel:
    """One plot: draw(fig, ax) fills ax (and may add a colour bar to fig)."""

    name: str                  # file-name suffix, e.g. "logG_histogram"
    draw: object               # callable(fig, ax)
    aspect: float = 0.78       # exported height / width
    extra: dict = field(default_factory=dict)


def apply_style():
    """Screen style for the GUI."""
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.linewidth": 1.0,
        "pdf.fonttype": 42,          # editable text in exported PDFs
        "ps.fonttype": 42,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })


@contextmanager
def nature_style(fs: FigureSettings):
    import matplotlib as mpl

    small = max(5.0, fs.font_size - 1)
    rc = {
        "font.family": "sans-serif",
        "font.sans-serif": [fs.font, "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": fs.font_size,
        "axes.labelsize": fs.font_size,
        "axes.titlesize": fs.font_size,
        "xtick.labelsize": small,
        "ytick.labelsize": small,
        "legend.fontsize": small,
        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "axes.linewidth": fs.axes_width,
        "lines.linewidth": fs.line_width,
        "lines.markersize": 2.5,
        "patch.linewidth": fs.axes_width,
        "xtick.major.width": fs.axes_width,
        "ytick.major.width": fs.axes_width,
        "xtick.minor.width": fs.axes_width * 0.8,
        "ytick.minor.width": fs.axes_width * 0.8,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.direction": fs.tick_direction,
        "ytick.direction": fs.tick_direction,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        # subscripts / symbols (G$_0$, beta) in the same font as the text
        "mathtext.fontset": "custom",
        "mathtext.rm": fs.font,
        "mathtext.it": f"{fs.font}:italic",
        "mathtext.bf": f"{fs.font}:bold",
        "mathtext.sf": fs.font,
        "mathtext.default": "regular",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": fs.dpi,
        "figure.dpi": 150,
    }
    with mpl.rc_context(rc):
        yield


def export_panels(panels, folder: str, prefix: str, fs: FigureSettings | None = None) -> list[str]:
    """Write every panel as its own Nature-format figure; returns the written paths."""
    from matplotlib.figure import Figure

    fs = fs or FigureSettings()
    os.makedirs(folder, exist_ok=True)
    written = []
    with nature_style(fs):
        for p in panels:
            w = fs.width_in()
            fig = Figure(figsize=(w, w * p.aspect), layout="constrained")
            fig.export = True            # panels skip on-screen titles
            ax = fig.add_subplot(111)
            p.draw(fig, ax)
            base = os.path.join(folder, f"{prefix}_{p.name}")
            for fmt in fs.format_list():
                path = f"{base}.{fmt}"
                fig.savefig(path, dpi=fs.dpi, bbox_inches="tight", pad_inches=0.02)
                written.append(path)
    return written


def save_figure(fig, base_path: str, formats=("png",)) -> list[str]:
    """Save a whole (composite) figure as base_path.<fmt>."""
    root, _ = os.path.splitext(base_path)
    out = []
    for fmt in formats:
        p = f"{root}.{fmt.lower().lstrip('.')}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        out.append(p)
    return out


def is_export(ax) -> bool:
    return bool(getattr(ax.figure, "export", False))


def title(ax, text: str, **kw):
    """Panel title on screen only (Nature figures carry no panel titles)."""
    if not is_export(ax):
        ax.set_title(text, fontsize=kw.pop("fontsize", 9), **kw)


def light_palette(app):
    """Fusion style with an explicit light palette, independent of the OS theme."""
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor as C, QPalette as P

    app.setStyle("Fusion")
    p = P()
    p.setColor(P.Window, C(240, 240, 240))
    p.setColor(P.WindowText, Qt.black)
    p.setColor(P.Base, C(255, 255, 255))
    p.setColor(P.AlternateBase, C(233, 231, 227))
    p.setColor(P.ToolTipBase, Qt.white)
    p.setColor(P.ToolTipText, Qt.black)
    p.setColor(P.Text, Qt.black)
    p.setColor(P.Button, C(240, 240, 240))
    p.setColor(P.ButtonText, Qt.black)
    p.setColor(P.BrightText, Qt.red)
    p.setColor(P.Highlight, C(76, 163, 224))
    p.setColor(P.HighlightedText, Qt.black)
    p.setColor(P.Disabled, P.Text, C(150, 150, 150))
    p.setColor(P.Disabled, P.ButtonText, C(150, 150, 150))
    app.setPalette(p)
