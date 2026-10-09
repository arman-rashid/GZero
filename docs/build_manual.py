"""Makes docs/User_Manual.pdf out of User_Manual.md.

The ```latex blocks are drawn with matplotlib, Qt does the layout and prints
the pdf.

    py docs/build_manual.py
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(HERE, "User_Manual.md")
PDF = os.path.join(HERE, "User_Manual.pdf")
IMG = os.path.join(HERE, "img")


def render_formulas(md: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(IMG, exist_ok=True)
    count = 0

    def repl(m):
        nonlocal count
        count += 1
        tex = m.group(1).strip()
        # mathtext has no \qquad spacing inside \mathrm etc.; keep it simple
        tex = tex.replace(r"\qquad", r"\quad\quad")
        name = f"eq{count}.png"
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.text(0, 0, f"${tex}$", fontsize=13)
        fig.savefig(os.path.join(IMG, name), dpi=220, bbox_inches="tight", pad_inches=0.06,
                    transparent=False, facecolor="white")
        plt.close(fig)
        return f"![formula {count}](img/{name})"

    return re.sub(r"```latex\n(.*?)\n```", repl, md, flags=re.S)


def build():
    from PySide6.QtCore import QMarginsF, QSizeF, QUrl
    from PySide6.QtGui import QFont, QPageLayout, QPageSize, QTextDocument
    from PySide6.QtPrintSupport import QPrinter
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv)  # noqa: F841 (fonts need an app)
    md = open(MD, encoding="utf-8").read()
    md = render_formulas(md)

    doc = QTextDocument()
    doc.setBaseUrl(QUrl.fromLocalFile(HERE + os.sep))
    doc.setDefaultFont(QFont("Segoe UI", 10))
    doc.setDefaultStyleSheet("table { border-collapse: collapse; } td, th { padding: 3px; }")
    doc.setMarkdown(md)
    # formula images are rendered at 220 dpi: give them their physical size,
    # capped at the text width, so they stay sharp and never overflow the page
    from PySide6.QtGui import QImage
    html = doc.toHtml()

    def size_img(m):
        src = m.group(1)
        im = QImage(os.path.join(HERE, src))
        w = min(im.width() * 96 / 220 * 0.85, 600) if not im.isNull() else 600
        return f'<img src="{src}" width="{int(w)}"'

    html = re.sub(r'<img src="([^"]+)"', size_img, html)
    doc.setHtml(html)

    printer = QPrinter(QPrinter.HighResolution)
    printer.setOutputFormat(QPrinter.PdfFormat)
    printer.setOutputFileName(PDF)
    printer.setPageLayout(QPageLayout(QPageSize(QPageSize.A4), QPageLayout.Portrait,
                                      QMarginsF(18, 16, 18, 16), QPageLayout.Millimeter))
    printer.setDocName("GZero - User Manual")
    doc.setPageSize(QSizeF(printer.pageRect(QPrinter.Point).size()))
    doc.print_(printer)
    print(f"wrote {PDF}")


if __name__ == "__main__":
    build()
