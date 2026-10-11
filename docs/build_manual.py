"""Makes docs/User_Manual.pdf out of docs/User_Manual.md with pandoc and LaTeX.

    py docs/build_manual.py

Needs pandoc (https://pandoc.org) and a LaTeX distribution with latexmk
(MiKTeX or TeX Live). The layout is in docs/latex/preamble.tex, the
conversion rules (formulas, table widths, title) in docs/latex/manual.lua.
The intermediate docs/latex/User_Manual.tex is kept, so the PDF can also be
built with any LaTeX editor.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LATEX = os.path.join(HERE, "latex")
MD = os.path.join(HERE, "User_Manual.md")
TEX = os.path.join(LATEX, "User_Manual.tex")
PDF = os.path.join(HERE, "User_Manual.pdf")


def run(cmd, cwd):
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def build():
    for tool in ("pandoc", "latexmk"):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not found: install pandoc and a LaTeX distribution (MiKTeX or TeX Live).")
    run(["pandoc", MD, "-f", "markdown+pipe_tables-implicit_figures", "-t", "latex", "-s", "--toc",
         "-o", TEX, "--lua-filter", os.path.join(LATEX, "manual.lua"),
         "-H", os.path.join(LATEX, "preamble.tex"),
         "-V", "documentclass=article", "-V", "fontsize=11pt", "-V", "papersize=a4",
         "-V", "date=", "-V", "colorlinks=true", "-V", "toc-title=Contents"], HERE)
    out = os.path.join(LATEX, "build")
    run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", f"-outdir={out}",
         os.path.basename(TEX)], LATEX)
    shutil.copyfile(os.path.join(out, "User_Manual.pdf"), PDF)
    print(f"wrote {PDF}")


if __name__ == "__main__":
    build()
