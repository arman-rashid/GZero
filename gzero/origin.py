"""Sending the result plots to Origin (OriginLab) instead of saving matplotlib figures.

Every panel is drawn exactly as for the figure export (style.export_panels,
same size and Nature style). The finished matplotlib axes are then read back
-- lines, steps, fills, scatter, bars, heatmaps, contours, labels, legend,
limits, scales -- and rebuilt in Origin: one workbook (or matrix book) with
the plotted data and one graph per panel, with the same page size, axis box,
fonts and line widths. So the panels in figures.py need no Origin-specific code.

Needs Windows, a licensed Origin 2021 or newer and the originpro package
(py -m pip install originpro). Origin is driven through COM, so call this
from the GUI thread.
"""

from __future__ import annotations

import os
import re
import struct
import sys
import tempfile
from dataclasses import dataclass, field

import numpy as np

from .style import FigureSettings, nature_style

PAGE_PAD_IN = 0.02          # same padding as export_panels (pad_inches)
DASHES = {"--": 1, "dashed": 1, ":": 2, "dotted": 2, "-.": 3, "dashdot": 3}


def unavailable_reason() -> str | None:
    """None if Origin can be used, else a message saying why not."""
    if sys.platform != "win32":
        return "Origin output needs Windows (Origin runs only on Windows)."
    try:
        import originpro  # noqa: F401
    except ImportError:
        return ("The Python package 'originpro' is not installed.\n"
                "Install it with:  py -m pip install originpro\n"
                "(Origin 2021 or newer must be installed as well.)")
    return None


# Reading the matplotlib axes

@dataclass
class Series:
    kind: str                    # "line", "scatter", "column", "fill"
    x: np.ndarray
    y: np.ndarray
    label: str = ""              # legend text ("" = not in the legend)
    color: str = "#000000"
    width: float = 0.75          # line width (pt)
    dash: int = 0                # 0 solid, 1 dash, 2 dot, 3 dash-dot
    alpha: float = 1.0
    size: float = 3.0            # symbol size (pt)
    fill: str | None = None      # fill colour of a closed polygon
    gap: float = 0.0             # column gap (%)
    raw: tuple | None = None     # (x, y) before step expansion, kept in the sheet


@dataclass
class Heatmap:
    xc: np.ndarray               # cell centres
    yc: np.ndarray
    Z: np.ndarray                # (len(yc), len(xc))
    cmap: object
    vmin: float
    vmax: float
    cbar_label: str = ""
    cbar_box: tuple | None = None   # colour bar axes (left, top, width, height) in inches from page top-left
    cbar_label_pos: tuple | None = None   # colour bar label (left, top) in colour bar data coordinates


@dataclass
class Text:
    text: str
    x: float                     # left / top in data coordinates
    y: float
    color: str
    size: float


@dataclass
class PanelData:
    name: str
    page: tuple                  # page width, height (inches)
    box: tuple                   # axes left, top, width, height in % of the page
    xlabel: str
    ylabel: str
    xlim: tuple
    ylim: tuple
    xscale: str
    yscale: str
    xstep: float | None
    ystep: float | None
    xticklabels: list | None     # category labels (pos, text, x, y), e.g. violin plot, else None
    series: list = field(default_factory=list)
    heatmap: Heatmap | None = None
    texts: list = field(default_factory=list)
    legend: tuple | None = None  # (x, y of the centre in data, title)
    font_size: float = 7.0


def _hex(c) -> str:
    from matplotlib.colors import to_hex
    return to_hex(c, keep_alpha=False)


def _alpha(c, artist_alpha=None) -> float:
    """Opacity; face/edge colours of patches and collections already include the artist alpha."""
    from matplotlib.colors import to_rgba
    a = to_rgba(c)[3]
    return a * (artist_alpha if artist_alpha is not None else 1.0)


def _first(colors, default="k"):
    colors = np.atleast_2d(colors) if len(colors) else None
    return default if colors is None or len(colors) == 0 else colors[0]


def _to_data(ax, transform, pts):
    """Points given in an artist's transform -> data coordinates of ax."""
    pts = np.asarray(pts, float)
    if transform is ax.transData:
        return pts
    return ax.transData.inverted().transform(transform.transform(pts))


def _label(artist) -> str:
    lab = artist.get_label() or ""
    return "" if lab.startswith("_") else lab


def _steps(x, y, where):
    """matplotlib step drawing as explicit staircase points."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2:
        return x, y
    if where == "steps-mid":
        m = 0.5 * (x[1:] + x[:-1])
        xs = np.concatenate([[x[0]], np.repeat(m, 2), [x[-1]]])
        ys = np.repeat(y, 2)
    elif where == "steps-post":
        xs = np.repeat(x, 2)[1:]
        ys = np.repeat(y, 2)[:-1]
    else:                                   # "steps" / "steps-pre"
        xs = np.repeat(x, 2)[:-1]
        ys = np.repeat(y, 2)[1:]
    return xs, ys


def _closed(d):
    """Polygon with the first point repeated at the end (else Origin draws a line from the end to the axis)."""
    return d if len(d) < 2 or np.array_equal(d[0], d[-1]) else np.vstack([d, d[:1]])


def _segments(paths, transform, ax):
    """Paths -> one x, y pair with NaN between the pieces (Origin breaks lines at missing values)."""
    xs, ys = [], []
    for p in paths:
        for poly in p.to_polygons(closed_only=False):
            if len(poly) < 2:
                continue
            d = _to_data(ax, transform, poly)
            xs += [d[:, 0], [np.nan]]
            ys += [d[:, 1], [np.nan]]
    if not xs:
        return None, None
    return np.concatenate(xs[:-1]), np.concatenate(ys[:-1])


def _line_series(ax, ln) -> Series | None:
    if not ln.get_visible():
        return None
    xy = _to_data(ax, ln.get_transform(), np.column_stack([ln.get_xdata(orig=False), ln.get_ydata(orig=False)]))
    x, y = xy[:, 0], xy[:, 1]
    ls, marker = ln.get_linestyle(), ln.get_marker()
    has_line = ls not in ("None", "", " ", "none")
    has_marker = marker not in ("None", "", " ", None, "none")
    s = Series("line", x, y, _label(ln), _hex(ln.get_color()), ln.get_linewidth(), DASHES.get(ls, 0),
               _alpha(ln.get_color(), ln.get_alpha()))
    if has_marker:
        s.kind = "linesymbol" if has_line else "scatter"
        s.size = max(1.5, ln.get_markersize() * (0.5 if marker in (".", ",") else 1.0))
        s.color = _hex(ln.get_markerfacecolor() if not has_line else ln.get_color())
    ds = ln.get_drawstyle()
    if ds and ds.startswith("steps") and has_line:
        s.raw = (x, y)
        s.x, s.y = _steps(x, y, ds)
    return s


def _collection_series(ax, c) -> list:
    from matplotlib.collections import LineCollection, PathCollection, PolyCollection
    if not c.get_visible():
        return []
    if isinstance(c, PathCollection):                   # scatter
        off = _to_data(ax, c.get_offset_transform(), c.get_offsets())
        sizes = c.get_sizes()
        fc = _first(c.get_facecolor())
        return [Series("scatter", off[:, 0], off[:, 1], _label(c), _hex(fc),
                       size=float(np.sqrt(sizes[0])) if len(sizes) else 3.0, alpha=_alpha(fc))]
    if isinstance(c, PolyCollection):                   # fill_between, violin bodies
        out = []
        fcs, ecs = c.get_facecolor(), c.get_edgecolor()
        for i, p in enumerate(c.get_paths()):
            for poly in p.to_polygons(closed_only=False):
                d = _closed(_to_data(ax, c.get_transform(), poly))
                fc = fcs[i % len(fcs)] if len(fcs) else (0, 0, 0, 0)
                ec = ecs[i % len(ecs)] if len(ecs) else fc
                lw = c.get_linewidth()[i % len(c.get_linewidth())] if len(c.get_linewidth()) else 0
                out.append(Series("fill", d[:, 0], d[:, 1], _label(c) if i == 0 else "",
                                  _hex(ec if lw > 0 else fc), lw, alpha=_alpha(fc), fill=_hex(fc)))
        return out
    paths = c.get_paths()                               # LineCollection, contour sets
    if not paths:
        return []
    ecs = c.get_edgecolor()
    lws = c.get_linewidth() if len(np.atleast_1d(c.get_linewidth())) else [0.75]
    lws = np.atleast_1d(lws)
    if isinstance(c, LineCollection) or len(ecs) <= 1:
        groups = [(paths, ecs[0] if len(ecs) else "k", lws[0])]
    else:                                               # one series per contour level
        groups = [([p], ecs[i % len(ecs)], lws[i % len(lws)]) for i, p in enumerate(paths)]
    out = []
    for ps, ec, lw in groups:
        x, y = _segments(ps, c.get_transform(), ax)
        if x is not None:
            out.append(Series("line", x, y, "", _hex(ec), float(lw), alpha=_alpha(ec)))
    return out


def _patch_series(ax, p) -> Series | None:
    if not p.get_visible():
        return None
    path = p.get_path()
    polys = path.to_polygons(closed_only=False)
    if not polys:
        return None
    d = _to_data(ax, p.get_transform(), polys[0])
    lw = p.get_linewidth()
    if p.get_fill():
        d = _closed(d)
        fc = p.get_facecolor()
        return Series("fill", d[:, 0], d[:, 1], _label(p), _hex(p.get_edgecolor() if lw > 0 else fc), lw,
                      alpha=_alpha(fc), fill=_hex(fc))
    ec = p.get_edgecolor()
    return Series("line", d[:, 0], d[:, 1], _label(p), _hex(ec), lw, DASHES.get(p.get_linestyle(), 0),
                  _alpha(ec))


def _bar_series(ax, bc) -> Series | None:
    rects = [r for r in bc.patches if r.get_visible()]
    if not rects:
        return None
    x = np.array([r.get_x() + r.get_width() / 2 for r in rects])
    h = np.array([r.get_height() for r in rects])
    w = np.median([r.get_width() for r in rects])
    pitch = np.median(np.diff(x)) if len(x) > 1 else w
    fc = rects[0].get_facecolor()
    return Series("column", x, h, _label(bc), _hex(fc), alpha=_alpha(fc),
                  gap=float(np.clip(100 * (1 - w / pitch), 0, 100)) if pitch > 0 else 0.0)


def _ticks_step(ax, which):
    if getattr(ax, f"get_{which}scale")() != "linear":
        return None
    t = getattr(ax, f"get_{which}ticks")()
    lo, hi = sorted(getattr(ax, f"get_{which}lim")())
    t = t[(t >= lo - 1e-9 * abs(hi - lo)) & (t <= hi + 1e-9 * abs(hi - lo))]
    return float(np.median(np.diff(t))) if len(t) > 1 else None


def _category_labels(ax, renderer):
    """Text tick labels set with set_xticks(pos, labels) -> [(pos, text, x_left, y_top)], else None."""
    from matplotlib.ticker import FixedFormatter, FuncFormatter     # set_xticks(labels=) uses either
    if not isinstance(ax.xaxis.get_major_formatter(), (FixedFormatter, FuncFormatter)):
        return None
    inv = ax.transData.inverted()
    labs = []
    for t in ax.get_xticklabels():
        if t.get_text():
            bb = t.get_window_extent(renderer)
            x, y = inv.transform((bb.x0, bb.y1))
            labs.append((float(t.get_position()[0]), t.get_text(), float(x), float(y)))
    try:
        [float(lab[1].replace("−", "-")) for lab in labs]
        return None                                   # numbers: Origin's own tick labels do the same
    except ValueError:
        return labs


def read_axes(fig, ax, name: str, fs: FigureSettings) -> PanelData:
    """Turn a drawn matplotlib panel into PanelData (fig must be drawn already)."""
    from matplotlib.collections import Collection, QuadMesh
    from matplotlib.container import BarContainer
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    renderer = fig.canvas.get_renderer()
    tight = fig.get_tightbbox(renderer).padded(PAGE_PAD_IN)          # inches
    W, H = tight.width, tight.height
    fw, fh = fig.get_size_inches()

    def box(a):                       # axes -> (left, top, width, height) in inches from the page top-left
        p = a.get_position()
        return p.x0 * fw - tight.x0, tight.y1 - p.y1 * fh, p.width * fw, p.height * fh

    l, t, w, h = box(ax)
    pd = PanelData(name, (W, H), (100 * l / W, 100 * t / H, 100 * w / W, 100 * h / H),
                   ax.get_xlabel(), ax.get_ylabel(), ax.get_xlim(), ax.get_ylim(),
                   ax.get_xscale(), ax.get_yscale(), _ticks_step(ax, "x"), _ticks_step(ax, "y"),
                   _category_labels(ax, renderer), font_size=fs.font_size)

    in_bars = set()
    bars = {}
    for bc in ax.containers:
        if isinstance(bc, BarContainer):
            in_bars.update(id(r) for r in bc.patches)
            bars[id(bc.patches[0])] = bc

    fills_first, rest = [], []
    for a in sorted(ax.get_children(), key=lambda a: a.get_zorder()):   # stable: keeps draw order
        if a is ax.patch or a in ax.spines.values():
            continue
        found = []
        if isinstance(a, Line2D):
            found = [_line_series(ax, a)]
        elif isinstance(a, QuadMesh):
            pd.heatmap = _heatmap(ax, a, box, renderer)
        elif isinstance(a, Collection):
            found = _collection_series(ax, a)
        elif isinstance(a, Patch):
            if id(a) in bars:
                found = [_bar_series(ax, bars[id(a)])]
            elif id(a) not in in_bars:
                found = [_patch_series(ax, a)]
        for s in found:
            if s is None or len(s.x) == 0:
                continue
            # translucent fills (spans, violins) go below the lines, as they look in matplotlib
            (fills_first if s.kind == "fill" and s.alpha < 1 else rest).append(s)
    pd.series = fills_first + rest

    inv = ax.transData.inverted()
    for txt in ax.texts:
        if not txt.get_visible() or not txt.get_text():
            continue
        bb = txt.get_window_extent(renderer)
        x, y = inv.transform((bb.x0, bb.y1))
        pd.texts.append(Text(txt.get_text(), float(x), float(y), _hex(txt.get_color()), txt.get_fontsize()))

    leg = ax.get_legend()
    if leg is not None and leg.get_visible():
        bb = leg.get_window_extent(renderer)
        x, y = inv.transform(((bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2))
        pd.legend = (float(x), float(y), leg.get_title().get_text())
    return pd


def _heatmap(ax, mesh, box, renderer) -> Heatmap:
    coords = mesh.get_coordinates()                    # (ny+1, nx+1, 2) cell corners
    xe, ye = coords[0, :, 0], coords[:, 0, 1]
    Z = np.ma.filled(np.ma.asarray(mesh.get_array(), float), np.nan).reshape(len(ye) - 1, len(xe) - 1)
    cb = getattr(mesh, "colorbar", None)
    vmin, vmax = mesh.get_clim()
    hm = Heatmap(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), Z, mesh.cmap, float(vmin), float(vmax))
    if cb is not None:
        hm.cbar_label, hm.cbar_box = cb.ax.get_ylabel(), box(cb.ax)
        if hm.cbar_label:                              # left, top of the label in colour bar data coordinates
            bb = cb.ax.yaxis.label.get_window_extent(renderer)
            hm.cbar_label_pos = tuple(float(v) for v in cb.ax.transData.inverted().transform((bb.x0, bb.y1)))
    return hm


def draw_panels(panels, fs: FigureSettings) -> list:
    """Draw every panel like export_panels does and read it back as PanelData."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    out = []
    with nature_style(fs):
        for p in panels:
            w = fs.width_in()
            fig = Figure(figsize=(w, w * p.aspect), layout="constrained")
            FigureCanvasAgg(fig)
            fig.export = True
            ax = fig.add_subplot(111)
            p.draw(fig, ax)
            fig.canvas.draw()                          # runs the constrained layout
            out.append(read_axes(fig, ax, p.name, fs))
    return out


# matplotlib text -> Origin text

GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "Delta": "Δ", "epsilon": "ε", "kappa": "κ",
         "lambda": "λ", "mu": "μ", "pi": "π", "rho": "ρ", "sigma": "σ", "tau": "τ", "phi": "φ", "omega": "ω",
         "Omega": "Ω", "pm": "±", "times": "×", "cdot": "·", "circ": "°", "infty": "∞", "approx": "≈",
         "leq": "≤", "geq": "≥", "propto": "∝"}


def origin_text(s: str) -> str:
    """Mathtext ($G_0$, $\\beta$, nm$^{-1}$) -> Origin escapes (G\\-(0), β, nm\\+(-1))."""
    parts = s.split("$")
    out = []
    for i, part in enumerate(parts):
        if i % 2 == 0:
            out.append(part)
            continue
        part = re.sub(r"\\(mathrm|mathit|text)\{([^{}]*)\}", r"\2", part)
        part = re.sub(r"\\([A-Za-z]+)", lambda m: GREEK.get(m.group(1), m.group(1)), part)
        part = re.sub(r"_\{([^{}]*)\}", r"\\-(\1)", part)
        part = re.sub(r"\^\{([^{}]*)\}", r"\\+(\1)", part)
        part = re.sub(r"_(\S)", r"\\-(\1)", part)
        part = re.sub(r"\^(\S)", r"\\+(\1)", part)
        out.append(part.replace("{", "").replace("}", ""))
    return "".join(out)


# Writing to Origin

def _write_pal(cmap, folder) -> str:
    """Colormap as a RIFF .pal file (256 colours) that Origin can load, so heatmaps keep their colours."""
    rgb = (np.array([cmap(i / 255)[:3] for i in range(256)]) * 255).round().astype(int)
    data = struct.pack("<HH", 0x300, 256) + b"".join(struct.pack("BBBB", *c, 0) for c in rgb)
    body = b"PAL data" + struct.pack("<I", len(data)) + data
    name = re.sub(r"[^A-Za-z0-9_]+", "_", getattr(cmap, "name", "cmap"))
    path = os.path.join(folder, f"GZero_{name}.pal")
    with open(path, "wb") as fh:
        fh.write(b"RIFF" + struct.pack("<I", len(body)) + body)
    return path


def _nice_levels(vmin, vmax, n=5):
    """Major colour-scale levels at round numbers covering [vmin, vmax]."""
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        return [vmin, vmin + 1]
    raw = (vmax - vmin) / n
    mag = 10 ** np.floor(np.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=raw)
    lo, hi = np.floor(vmin / step) * step, np.ceil(vmax / step - 1e-9) * step
    return list(np.round(np.arange(lo, hi + step / 2, step), 12))


def _rgb(hexcol):
    return tuple(int(hexcol[i:i + 2], 16) for i in (1, 3, 5))


class OriginWriter:
    def __init__(self, op, fs: FigureSettings):
        self.op, self.fs = op, fs
        self.small = max(5.0, fs.font_size - 1)
        self.tmp = os.path.join(tempfile.gettempdir(), "gzero_origin")     # colour map .pal files
        os.makedirs(self.tmp, exist_ok=True)

    def lt(self, page, cmd):
        self.op.lt_exec(f"win -a {page.name}; {cmd}")

    def panel(self, pd: PanelData, title: str):
        op = self.op
        if pd.heatmap is not None:
            g = self._heatmap_graph(pd, title)
        else:
            g = op.new_graph(lname=title)
        gl = g[0]
        plots = self._series(gl, pd, title) if pd.series else []
        self._layout(g, gl, pd)
        self._legend(g, gl, pd, plots)        # last: Origin rebuilds the legend when plots or axes change
        return g

    def _series(self, gl, pd, title):
        op = self.op
        wks = op.new_sheet("w", lname=title)
        col = 0
        plots = []
        for k, s in enumerate(pd.series):
            name = s.label or f"{s.kind} {k + 1}"
            if s.raw is not None:                       # original data first (e.g. histogram bins)
                wks.from_list(col, s.raw[0], lname=origin_text(pd.xlabel), axis="X")
                wks.from_list(col + 1, s.raw[1], lname=origin_text(name), comments="data", axis="Y")
                col += 2
            wks.from_list(col, s.x, lname=origin_text(pd.xlabel), axis="X",
                          comments="step outline" if s.raw is not None else "")
            wks.from_list(col + 1, s.y, lname=origin_text(name), axis="Y",
                          comments="step outline" if s.raw is not None else "")
            ptype = {"scatter": "s", "linesymbol": "y", "column": "c"}.get(s.kind, "l")
            p = gl.add_plot(wks, col + 1, col, type=ptype)
            col += 2
            if p is None:
                continue
            plots.append((p, s))
            p.color = s.color
            cmds = [f"-w {max(s.width, 0.05) * 500:.0f}", f"-d {s.dash}"]
            if s.kind in ("scatter", "linesymbol"):
                p.symbol_kind = 2
                p.symbol_interior = 0
                cmds.append(f"-z {max(s.size, 0.5):.2f}")
            if s.kind == "fill":
                cmds += ["-pf 1", "-pfv 2", "-pfb color({},{},{})".format(*_rgb(self._blend(s.fill, s.alpha)))]
                if s.width <= 0:
                    p.color = self._blend(s.fill, s.alpha)
            if s.kind == "column":
                cmds += [f"-vg {s.gap:.0f}", "-pfb color({},{},{})".format(*_rgb(self._blend(s.color, s.alpha)))]
            p.set_cmd(*cmds)
            if s.alpha < 1 and s.kind in ("line", "scatter", "linesymbol"):
                p.transparency = int(round(100 * (1 - s.alpha)))
        return plots

    @staticmethod
    def _blend(hexcol, alpha):
        """Translucent colour on white as an opaque colour (Origin fills have no alpha here)."""
        if alpha >= 1:
            return hexcol
        c = np.array(_rgb(hexcol)) * alpha + 255 * (1 - alpha)
        return "#" + "".join(f"{int(round(v)):02x}" for v in c)

    def _legend(self, g, gl, pd, plots):
        first = 2 if pd.heatmap is not None else 1           # the heatmap itself is plot 1
        entries = [(first + i, s.label) for i, (_, s) in enumerate(plots) if s.label]
        if pd.legend is None or not entries:
            self.lt(g, "label -r legend;")
            return
        self.lt(g, "legend -r;")
        leg = gl.label("Legend")
        if leg is None:
            return
        lines = ([origin_text(pd.legend[2])] if pd.legend[2] else []) + \
                [f"\\l({i}) {origin_text(t)}" for i, t in entries]
        leg.text = "\n".join(lines)
        self.lt(g, f"legend.fsize = {self.small}; legend.font = font({self.fs.font}); legend.background = 0; "
                   f"legend.x = {pd.legend[0]!r}; legend.y = {pd.legend[1]!r};")

    def _heatmap_graph(self, pd, title):
        op = self.op
        hm = pd.heatmap
        ms = op.new_sheet("m", lname=title)
        ms.from_np(np.ascontiguousarray(hm.Z, dtype=float))
        ms.xymap = float(hm.xc[0]), float(hm.xc[-1]), float(hm.yc[0]), float(hm.yc[-1])
        # default template (same axis layout as the other graphs); Origin's colour scale object has no
        # scriptable label font, so the colour bar is a second layer (see _colorbar)
        g = op.new_graph(lname=title)
        self._colormap(g[0], ms, hm)
        return g

    def _colormap(self, gl, ms, hm):
        """Heatmap of matrix sheet ms in layer gl with the matplotlib colormap and limits."""
        p = gl.add_mplot(ms, 0, type=105)
        gl.rescale()                                     # initialises the colour map of the new plot
        majors = _nice_levels(hm.vmin, hm.vmax)
        p.zlevels = {"minors": max(1, 100 // max(1, len(majors) - 1) - 1), "levels": majors}
        p.colormap = _write_pal(hm.cmap, self.tmp)       # after the levels, which reset the colours
        # clipped values get the end colours and missing values (NaN) stay white, as in matplotlib
        for prop, c in (("colorbelow", _hex(hm.cmap(0.0))), ("colorabove", _hex(hm.cmap(1.0))), ("colorMiss", "#ffffff")):
            gl.obj.SetNumProp(f"plot1.cmap.{prop}", self.op.ocolor(c))
        return majors

    def _colorbar(self, g, pd, res):
        """Colour bar as a narrow second layer at the place of the matplotlib colour bar."""
        hm = pd.heatmap
        W, H = pd.page
        l, t, w, h = hm.cbar_box
        majors = _nice_levels(hm.vmin, hm.vmax)
        vals = np.linspace(majors[0], majors[-1], 256)
        ms = self.op.new_sheet("m", lname=f"{g.lname} colour bar", hidden=True)
        ms.from_np(np.column_stack([vals, vals]))
        ms.xymap = 0.0, 1.0, float(vals[0]), float(vals[-1])
        gl = g.add_layer(0)
        self._colormap(gl, ms, hm)
        gl.set_xlim(0, 1)
        gl.set_ylim(hm.vmin, hm.vmax, majors[1] - majors[0])
        fs = self.fs
        tick = 2 if fs.tick_direction == "out" else 1
        self.lt(g, f"page.active = 2; layer.unit = 1; layer.left = {100 * l / W:.3f}; layer.top = {100 * t / H:.3f}; "
                   f"layer.width = {100 * w / W:.3f}; layer.height = {100 * h / H:.3f}; "
                   "layer.x.showAxes = 3; layer.y.showAxes = 3; layer.x.showLabels = 0; layer.x.ticks = 0; "
                   f"layer.y.showLabels = 2; layer.y.ticks = 0; layer.y2.ticks = {tick}; layer.y.minorTicks = 0; "
                   f"layer.y2.label.pt = {self.small}; layer.y2.label.font = font({fs.font}); layer.x2.ticks = 0; "
                   f"layer.x.thickness = {fs.axes_width}; layer.y.thickness = {fs.axes_width}; "
                   f"layer.x2.thickness = {fs.axes_width}; layer.y2.thickness = {fs.axes_width}; layer.tickL = 3; "
                   "layer.factor = 1; label -r xb; label -r yl; label -r legend;")
        if hm.cbar_label and hm.cbar_label_pos:
            lab = gl.add_label(origin_text(hm.cbar_label), *hm.cbar_label_pos)
            if lab is not None:
                lab.set_float("fsize", fs.font_size)
                lab.set_int("rotate", 90)
                lab.set_float("x1", hm.cbar_label_pos[0])      # again: rotating moves the box
                lab.set_float("y1", hm.cbar_label_pos[1])
        self.lt(g, "page.active = 1;")

    def _layout(self, g, gl, pd):
        fs = self.fs
        W, H = pd.page
        self.lt(g, "double __res = page.resx;")
        res = self.op.po.LT_get_var("__res") or 600
        lt = [f"page.width = {W * res:.0f}; page.height = {H * res:.0f};",
              "layer.unit = 1; layer.left = {:.3f}; layer.top = {:.3f}; layer.width = {:.3f}; layer.height = {:.3f};"
              .format(*pd.box)]
        tick = 10 if fs.tick_direction == "out" else 5
        for a in ("x", "y"):
            lt.append(f"layer.{a}.showAxes = 3; layer.{a}.label.pt = {self.small}; layer.{a}.label.font = font({fs.font}); "
                      f"layer.{a}.label.rotate = 0; layer.{a}.thickness = {fs.axes_width}; "
                      f"layer.{a}2.thickness = {fs.axes_width}; layer.{a}.ticks = {tick};")
        lt.append("layer.tickL = 3;")
        for obj in ("xb", "yl"):
            lt.append(f"{obj}.fsize = {fs.font_size}; {obj}.font = font({fs.font});")
        self.lt(g, " ".join(lt))

        for a, scale, lim, step in (("x", pd.xscale, pd.xlim, pd.xstep), ("y", pd.yscale, pd.ylim, pd.ystep)):
            setattr(gl, f"{a}scale", "log10" if scale == "log" else "linear")
            lo, hi = lim
            if scale == "log":
                getattr(gl, f"set_{a}lim")(lo, hi, 1)
                self.lt(g, f"layer.{a}.label.numFormat = 2;")        # 10^n, as matplotlib
            else:
                getattr(gl, f"set_{a}lim")(lo, hi, step) if step else getattr(gl, f"set_{a}lim")(lo, hi)
                self.lt(g, f"layer.{a}.minorTicks = 0;")
        gl.axis("x").title = origin_text(pd.xlabel)
        gl.axis("y").title = origin_text(pd.ylabel)
        if pd.xticklabels:
            self._category_ticks(g, gl, pd)

        for t in pd.texts:
            lab = gl.add_label(origin_text(t.text), t.x, t.y)
            if lab is not None:
                lab.color = t.color
                lab.set_float("fsize", self.small)

        # Origin scales all text of a layer by layer.factor, which changes when the layer is resized
        self.lt(g, "layer.factor = 1;")
        if pd.heatmap is not None and pd.heatmap.cbar_box is not None:
            self._colorbar(g, pd, res)

    def _category_ticks(self, g, gl, pd):
        """Text tick labels (violin plot): ticks at the positions, the texts as labels where matplotlib has them."""
        pos = [lab[0] for lab in pd.xticklabels]
        step = np.median(np.diff(pos)) if len(pos) > 1 else 1
        gl.set_xlim(pd.xlim[0], pd.xlim[1], step)
        self.lt(g, f"layer.x.showLabels = 0; layer.x.firstTick = {pos[0]};")
        for _, text, x, y in pd.xticklabels:
            lab = gl.add_label(origin_text(text), x, y)
            if lab is not None:
                lab.set_float("fsize", self.small)


def send(groups, fs: FigureSettings | None = None, project_path: str | None = None, progress=None) -> list[str]:
    """Rebuild panels in Origin; returns the graph names.

    groups: [(folder, prefix, panels)] -- each group goes into its own Project
    Explorer folder, graphs are named <prefix>_<panel name>. Uses a running
    Origin if there is one (else starts it) and leaves it open and visible.
    With project_path the Origin project is saved there as well.
    progress(done, total, name) is called after every graph.
    """
    reason = unavailable_reason()
    if reason:
        raise RuntimeError(reason)
    import originpro as op

    fs = fs or FigureSettings()
    total = sum(len(p) for _, _, p in groups)
    names = []
    op.attach()
    try:
        op.set_show(True)
        w = OriginWriter(op, fs)
        for folder, prefix, panels in groups:
            op.lt_exec('pe_cd "/";')
            if folder:
                op.lt_exec(f'pe_mkdir "{_safe(folder)}" chk:=1 cd:=1;')
            for p in panels:
                g = w.panel(draw_panels([p], fs)[0], f"{prefix}_{p.name}")
                names.append(g.lname)
                if progress:
                    progress(len(names), total, g.lname)
        op.lt_exec('pe_cd "/";')
        if project_path:
            # @SDO: save even if a dialog is open in Origin (else Origin asks first and the call never returns)
            op.lt_exec("__gz_sdo = @SDO; @SDO = 1;")
            ok = op.save(project_path)
            op.lt_exec("@SDO = __gz_sdo;")
            if not ok:
                raise RuntimeError(f"Origin could not save the project to {project_path}")
    finally:
        op.detach()
    return names


def _safe(name: str) -> str:
    """Project Explorer folder name without characters Origin does not accept."""
    return re.sub(r'[\\/:*?"<>|&]+', " ", name).strip() or "GZero"
