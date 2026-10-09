"""Qt helpers: plot canvas, a settings form built from a dataclass, background worker."""

from __future__ import annotations

import dataclasses
import traceback

import numpy as np
from PySide6.QtCore import QObject, QRunnable, QThreadPool, Qt, Signal
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QMessageBox, QPushButton, QScrollArea, QSizePolicy, QSpinBox, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from .. import style


# Matplotlib canvas

class PlotPanel(QWidget):
    """Figure + navigation toolbar."""

    def __init__(self, parent=None, figsize=(8, 6)):
        super().__init__(parent)
        self.figure = Figure(figsize=figsize, layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas, 1)

    def clear(self):
        self.figure.clear()

    def draw(self):
        self.canvas.draw_idle()

    def save(self, path_without_ext: str, formats=("png",)):
        return style.save_figure(self.figure, path_without_ext, formats)

    def show_panels(self, panels, ncols: int = 2, width_ratios=None):
        """Draw a list of :class:`style.Panel` in a grid (the same panels are exported one per file)."""
        self.figure.clear()
        n = len(panels)
        if n == 0:
            self.draw()
            return
        ncols = max(1, min(ncols, n))
        nrows = int(np.ceil(n / ncols))
        gs = self.figure.add_gridspec(nrows, ncols, width_ratios=width_ratios if width_ratios and len(width_ratios) == ncols else None)
        for i, p in enumerate(panels):
            ax = self.figure.add_subplot(gs[i // ncols, i % ncols])
            p.draw(self.figure, ax)
        self.draw()


from ..figures import decimate as decimate_for_plot  # noqa: E402,F401


# Dataclass -> form

class FloatEdit(QLineEdit):
    """Line edit for floats of any magnitude (1e-10 ... 1e6); empty = None if optional."""

    def __init__(self, value, optional=False):
        super().__init__()
        self.optional = optional
        v = QDoubleValidator()
        v.setNotation(QDoubleValidator.ScientificNotation)
        self.setValidator(v)
        self.set_value(value)
        if optional:
            self.setPlaceholderText("auto / off")

    def set_value(self, value):
        self.setText("" if value is None else f"{value:g}")

    def value(self):
        t = self.text().strip()
        if not t:
            if self.optional:
                return None
            raise ValueError("A value is required.")
        return float(t)


class TupleEdit(QLineEdit):
    def __init__(self, value, cast=float):
        super().__init__(", ".join(f"{v:g}" for v in value))
        self.cast = cast

    def set_value(self, value):
        self.setText(", ".join(f"{v:g}" for v in value))

    def value(self):
        return tuple(self.cast(float(s)) for s in self.text().replace(";", ",").split(",") if s.strip())


class SettingsForm(QGroupBox):
    """Builds an editable form from a dataclass instance.

    labels  : field -> label text (fields not listed get a prettified name)
    choices : field -> list of allowed strings (combo box)
    tips    : field -> tooltip
    hidden  : fields not shown (kept at their current value)
    Nested dataclass fields get their own sub-group.
    """

    changed = Signal()

    def __init__(self, title: str, instance, labels=None, choices=None, tips=None, hidden=(),
                 parent=None):
        super().__init__(title, parent)
        self._cls = type(instance)
        self._instance = instance
        self.labels = labels or {}
        self.choices = choices or {}
        self.tips = tips or {}
        self.hidden = set(hidden)
        self.widgets = {}
        self.sub = {}
        lay = QFormLayout(self)
        lay.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        lay.setLabelAlignment(Qt.AlignRight)
        for f in dataclasses.fields(instance):
            if f.name in self.hidden:
                continue
            val = getattr(instance, f.name)
            label = self.labels.get(f.name, f.name.replace("_", " "))
            if dataclasses.is_dataclass(val):
                sub = SettingsForm(label, val, labels=self.labels.get(f"{f.name}.labels"),
                                   choices=self.choices.get(f"{f.name}.choices"))
                sub.changed.connect(self.changed)
                self.sub[f.name] = sub
                lay.addRow(sub)
                continue
            w = self._make_widget(f.name, val, f.type)
            if f.name in self.tips:
                w.setToolTip(self.tips[f.name])
            self.widgets[f.name] = w
            lay.addRow(label, w)

    @staticmethod
    def _parse_type(hint):
        """(base type, optional) from a dataclass annotation.

        With from __future__ import annotations the annotation is a string
        such as "float | None"; Python 3.9 cannot evaluate that, so parse it.
        """
        if not isinstance(hint, str):
            return hint, False
        parts = [p.strip() for p in hint.split("|")]
        optional = "None" in parts
        names = {"float": float, "int": int, "bool": bool, "str": str, "tuple": tuple}
        base = next((names[p] for p in parts if p in names), None)
        return base, optional

    def _make_widget(self, name, val, hint):
        base, optional = self._parse_type(hint)
        if name in self.choices:
            w = QComboBox()
            w.addItems(self.choices[name])
            if val is not None and str(val) in self.choices[name]:
                w.setCurrentText(str(val))
            w.currentIndexChanged.connect(self.changed)
            return w
        if base is bool or isinstance(val, bool):
            w = QCheckBox()
            w.setChecked(bool(val))
            w.toggled.connect(self.changed)
            return w
        if (base is int or isinstance(val, int)) and not optional:
            w = QSpinBox()
            w.setRange(-1_000_000_000, 1_000_000_000)
            w.setValue(int(val))
            w.valueChanged.connect(self.changed)
            return w
        if base is tuple or isinstance(val, tuple):
            cast = int if val and all(isinstance(v, int) for v in val) else float
            w = TupleEdit(val, cast)
            w.editingFinished.connect(self.changed)
            return w
        if base is str or isinstance(val, str):
            w = QLineEdit(str(val))
            w.editingFinished.connect(self.changed)
            return w
        if optional and base is int:
            w = FloatEdit(val, optional=True)
            w._is_int = True
            w.editingFinished.connect(self.changed)
            return w
        w = FloatEdit(val, optional=optional)
        w.editingFinished.connect(self.changed)
        return w

    def value(self):
        kwargs = {}
        for f in dataclasses.fields(self._instance):
            if f.name in self.sub:
                kwargs[f.name] = self.sub[f.name].value()
                continue
            w = self.widgets.get(f.name)
            if w is None:
                kwargs[f.name] = getattr(self._instance, f.name)
            elif isinstance(w, QComboBox):
                kwargs[f.name] = w.currentText()
            elif isinstance(w, QCheckBox):
                kwargs[f.name] = w.isChecked()
            elif isinstance(w, QSpinBox):
                kwargs[f.name] = w.value()
            elif isinstance(w, (FloatEdit, TupleEdit)):
                try:
                    v = w.value()
                except ValueError:
                    raise ValueError(f"Invalid value for '{self.labels.get(f.name, f.name)}'.")
                if getattr(w, "_is_int", False) and v is not None:
                    v = int(v)
                kwargs[f.name] = v
            else:
                kwargs[f.name] = w.text()
        return self._cls(**kwargs)

    def set_value(self, instance):
        self._instance = instance
        for f in dataclasses.fields(instance):
            val = getattr(instance, f.name)
            if f.name in self.sub:
                self.sub[f.name].set_value(val)
                continue
            w = self.widgets.get(f.name)
            if w is None:
                continue
            w.blockSignals(True)
            if isinstance(w, QComboBox):
                w.setCurrentText(str(val))
            elif isinstance(w, QCheckBox):
                w.setChecked(bool(val))
            elif isinstance(w, QSpinBox):
                w.setValue(int(val))
            elif isinstance(w, (FloatEdit, TupleEdit)):
                w.set_value(val)
            else:
                w.setText(str(val))
            w.blockSignals(False)

    def set_choices(self, name, items, current=None):
        w = self.widgets.get(name)
        if isinstance(w, QComboBox):
            w.blockSignals(True)
            w.clear()
            w.addItems(items)
            if current is not None:
                w.setCurrentText(str(current))
            w.blockSignals(False)


# Layout helpers

def scroll_panel(*widgets, width=330) -> QScrollArea:
    """Vertical stack of widgets in a scroll area (left-hand control column)."""
    inner = QWidget()
    lay = QVBoxLayout(inner)
    lay.setContentsMargins(4, 4, 4, 4)
    for w in widgets:
        if isinstance(w, str):
            lab = QLabel(w)
            lab.setWordWrap(True)
            lay.addWidget(lab)
        elif w is None:
            continue
        else:
            lay.addWidget(w)
    lay.addStretch(1)
    area = QScrollArea()
    area.setWidget(inner)
    area.setWidgetResizable(True)
    area.setMinimumWidth(width)
    area.setMaximumWidth(width + 140)
    return area


def button_row(*buttons) -> QWidget:
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    for b in buttons:
        lay.addWidget(b)
    return w


def make_button(text, slot, tip=None) -> QPushButton:
    b = QPushButton(text)
    b.clicked.connect(slot)
    if tip:
        b.setToolTip(tip)
    return b


def error_box(parent, title, exc):
    msg = str(exc) if not isinstance(exc, str) else exc
    QMessageBox.warning(parent, title, msg)


# Background work

class _Signals(QObject):
    finished = Signal(object)
    failed = Signal(str)


class Worker(QRunnable):
    """Run fn(*args) in the global thread pool; emit result or error text."""

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs
        self.signals = _Signals()

    def run(self):
        try:
            res = self.fn(*self.args, **self.kwargs)
        except (ValueError, KeyError) as exc:   # expected problems: message only
            self.signals.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - report everything else with a traceback
            self.signals.failed.emit(f"{exc}\n\n{traceback.format_exc(limit=4)}")
        else:
            self.signals.finished.emit(res)


def run_in_background(fn, *args, on_done=None, on_error=None, **kwargs) -> Worker:
    w = Worker(fn, *args, **kwargs)
    if on_done:
        w.signals.finished.connect(on_done)
    if on_error:
        w.signals.failed.connect(on_error)
    QThreadPool.globalInstance().start(w)
    return w
