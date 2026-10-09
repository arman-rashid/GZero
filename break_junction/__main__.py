"""python -m break_junction starts the GUI (running this file directly also works)."""

import os
import sys

if not __package__:
    # run as a plain script: make the package importable and use absolute imports
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from break_junction.gui.app import main
else:
    from .gui.app import main

sys.exit(main())
