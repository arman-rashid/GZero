"""Makes the PDFs in docs/ out of their Markdown sources with pandoc and LaTeX:
User_Manual.md -> User_Manual.pdf and GZero_methods.md -> GZero_methods.pdf.

    py docs/build_manual.py

Needs pandoc (https://pandoc.org) and a LaTeX distribution with latexmk
(MiKTeX or TeX Live). The layout is in docs/latex/preamble.tex, the
conversion rules (formulas, table widths, title) in docs/latex/manual.lua.
The intermediate .tex files in docs/latex/ are kept, so the PDFs can also be
built with any LaTeX editor.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LATEX = os.path.join(HERE, "latex")
DOCUMENTS = [("User_Manual", True), ("GZero_methods", False)]   # (name, with contents)


def run(cmd, cwd):
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def build_one(name, toc):
    tex = os.path.join(LATEX, name + ".tex")
    cmd = ["pandoc", os.path.join(HERE, name + ".md"), "-f", "markdown+pipe_tables-implicit_figures",
           "-t", "latex", "-s", "-o", tex, "--lua-filter", os.path.join(LATEX, "manual.lua"),
           "-H", os.path.join(LATEX, "preamble.tex"), "-V", "documentclass=article", "-V", "fontsize=11pt",
           "-V", "papersize=a4", "-V", "date=", "-V", "colorlinks=true"]
    if toc:
        cmd += ["--toc", "-V", "toc-title=Contents"]
    run(cmd, HERE)
    out = os.path.join(LATEX, "build")
    run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", f"-outdir={out}",
         os.path.basename(tex)], LATEX)
    pdf = os.path.join(HERE, name + ".pdf")
    shutil.copyfile(os.path.join(out, name + ".pdf"), pdf)
    print(f"wrote {pdf}")


def build():
    for tool in ("pandoc", "latexmk"):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not found: install pandoc and a LaTeX distribution (MiKTeX or TeX Live).")
    for name, toc in DOCUMENTS:
        build_one(name, toc)


if __name__ == "__main__":
    build()
