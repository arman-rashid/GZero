"""Base class for the analysis tabs."""

from __future__ import annotations

import dataclasses
import os

from PySide6.QtWidgets import QApplication, QHBoxLayout, QWidget
from PySide6.QtCore import Qt

from .widgets import error_box, run_in_background


class AnalysisTab(QWidget):
    """A tab = left control column + right plot area.

    Subclasses fill self.forms (name -> SettingsForm) so settings are
    saved/loaded automatically, and implement export(folder).
    """

    title = "Tab"

    def __init__(self, project, main):
        super().__init__()
        self.project = project
        self.main = main
        self.forms = {}
        self.layout_ = QHBoxLayout(self)
        self.layout_.setContentsMargins(4, 4, 4, 4)
        self._workers = []

    # -- settings persistence ---------------------------------------------------
    def settings_dict(self) -> dict:
        out = {}
        for name, form in self.forms.items():
            try:
                out[name] = dataclasses.asdict(form.value())
            except ValueError:
                pass
        return out

    def apply_settings(self, d: dict):
        for name, form in self.forms.items():
            if name not in d:
                continue
            cls = type(form.value())
            try:
                form.set_value(_from_dict(cls, d[name]))
            except (TypeError, ValueError):
                pass

    # -- results ------------------------------------------------------------------
    def export(self, folder: str, prefix: str) -> list[str]:
        return []

    def export_panels(self) -> list:
        """Panels written as individual figures (default: the panels on screen)."""
        return getattr(self, "panels", lambda: [])()

    def figure_prefix(self, prefix: str) -> str:
        """Start of the figure file names (and Origin graph names)."""
        return prefix

    def export_figures(self, folder: str, prefix: str) -> list[str]:
        """Every panel as its own Nature-format figure in folder/figures."""
        from .. import style
        panels = self.export_panels()
        if not panels:
            return []
        return style.export_panels(panels, os.path.join(folder, "figures"), self.figure_prefix(prefix),
                                   self.main.figure_settings)

    def has_results(self) -> bool:
        return False

    # -- helpers ------------------------------------------------------------------
    def form_value(self, name):
        try:
            return self.forms[name].value()
        except ValueError as exc:
            error_box(self, "Invalid setting", exc)
            return None

    def status(self, msg: str):
        self.project.status.emit(msg)

    def background(self, fn, *args, on_done=None, busy="Working ..."):
        """Run fn in a worker thread with a busy cursor; on_done(result) in the GUI thread."""
        self.status(busy)
        QApplication.setOverrideCursor(Qt.WaitCursor)

        def done(res):
            QApplication.restoreOverrideCursor()
            self.status("Ready")
            if on_done:
                try:
                    on_done(res)
                except Exception as exc:  # noqa: BLE001
                    error_box(self, self.title, exc)

        def failed(msg):
            QApplication.restoreOverrideCursor()
            self.status("Error")
            error_box(self, self.title, msg)

        w = run_in_background(fn, *args, on_done=done, on_error=failed)
        self._workers.append(w)       # keep a reference until finished
        self._workers = self._workers[-8:]
        return w


def _from_dict(cls, d):
    kwargs = {}
    for f in dataclasses.fields(cls):
        if f.name not in d:
            continue
        v = d[f.name]
        default = f.default_factory() if f.default_factory is not dataclasses.MISSING else f.default
        if dataclasses.is_dataclass(default) and isinstance(v, dict):
            v = _from_dict(type(default), v)
        elif isinstance(default, tuple) and isinstance(v, list):
            v = tuple(v)
        kwargs[f.name] = v
    return cls(**kwargs)
