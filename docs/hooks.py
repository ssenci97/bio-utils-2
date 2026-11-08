"""MkDocs hooks: render docs/assets/operations.png from docs/scripts/operations_graph.py before each build."""

from __future__ import annotations

import importlib.util
import logging
from pathlib import Path
from typing import Any

LOG = logging.getLogger("mkdocs.hooks.operations_graph")
SCRIPT = Path(__file__).resolve().parent / "scripts" / "operations_graph.py"


def on_pre_build(config: Any) -> None:
    """Render the operations graph; keep the committed PNG when graphviz is missing.

    The PNG is rewritten only when its bytes change, so 'mkdocs serve' does not rebuild in a loop.

    Args:
        config: MkDocs configuration (unused).
    """
    spec = importlib.util.spec_from_file_location("operations_graph", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        module.render()
    except (ImportError, OSError) as error:   # graphviz package or dot executable missing
        LOG.info("operations graph not rendered (%s); keeping the committed PNG.", error)   # info: --strict fails on warnings
