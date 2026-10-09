"""Builds the standalone Windows version with PyInstaller.

    py -m pip install pyinstaller
    py build_exe.py

Everything ends up in dist/GZero/. Copy the whole folder, the exe
alone won't run.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    docs = os.path.join(HERE, "docs")
    args = [
        sys.executable, "-m", "PyInstaller",
        "--noconfirm", "--clean", "--windowed",
        "--name", "GZero",
        "--distpath", os.path.join(HERE, "dist"),
        "--workpath", os.path.join(HERE, "build"),
        "--specpath", os.path.join(HERE, "build"),
        "--paths", HERE,
        "--collect-submodules", "gzero",
        "--hidden-import", "igor2.packed",
        "--hidden-import", "igor2.binarywave",
        "--hidden-import", "sklearn.utils._typedefs",
        "--hidden-import", "sklearn.neighbors._partition_nodes",
        # figure export backends are loaded by file type at save time; PyInstaller cannot see them
        "--hidden-import", "matplotlib.backends.backend_pdf",
        "--hidden-import", "matplotlib.backends.backend_svg",
        "--hidden-import", "matplotlib.backends.backend_ps",
        "--hidden-import", "matplotlib.backends.backend_agg",
        "--exclude-module", "statsmodels",
        "--exclude-module", "tkinter",
        "--exclude-module", "PyQt5",
        "--exclude-module", "PyQt6",
    ]
    if os.path.isdir(docs):
        args += ["--add-data", f"{docs}{os.pathsep}docs"]
    try:                # optional: plotting in Origin (gzero/origin.py), bundled if installed
        import originpro  # noqa: F401
        args += ["--collect-all", "originpro", "--collect-all", "OriginExt"]
    except ImportError:
        print("originpro not installed: the exe will not be able to plot in Origin")
    args.append(os.path.join(HERE, "GZero.py"))
    print(" ".join(args))
    subprocess.check_call(args, cwd=HERE)
    # also place the manual next to the exe so it is easy to find
    out = os.path.join(HERE, "dist", "GZero")
    if os.path.isdir(docs):
        shutil.copytree(docs, os.path.join(out, "docs"), dirs_exist_ok=True)
    print(f"\nBuilt: {os.path.join(out, 'GZero.exe')}")


if __name__ == "__main__":
    main()
