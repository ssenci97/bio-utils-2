"""Type aliases shared across the package."""

from __future__ import annotations

import os
import sys
from collections.abc import Iterable, Sequence
from typing import Any

if sys.version_info >= (3, 11):
    from typing import Self
else:  # Python 3.10
    from typing_extensions import Self

__all__ = ["ChainPair", "ChainRef", "ChainRefs", "Moieties", "PathLike", "Record", "Self"]

PathLike = str | os.PathLike[str]
ChainRef = str
ChainRefs = Iterable[str]
ChainPair = tuple[str, str]
Moieties = Sequence[str]
Record = dict[str, Any]
