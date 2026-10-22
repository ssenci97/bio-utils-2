"""Second test layer: type hints and docstrings, checked statically and at runtime.

- test_every_function_is_annotated: needs nothing but the standard library.
- ruff (docstring style, annotation presence), pydoclint (docstrings match signatures) and mypy (static types): each
  is skipped when its tool is not installed (pip install -e ".[dev]").
- test_runtime_types: re-runs the offline suite under typeguard, which checks every annotated argument and return
  value of synced_structure at runtime.
"""

import ast
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "synced_structure"
OFFLINE_TESTS = ["tests/test_smoke.py", "tests/test_regressions.py", "tests/test_docking.py"]


def _run(command):
    """Run a command from the repository root and return (exit code, combined output)."""
    process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    return process.returncode, process.stdout + process.stderr


def _module_command(module):
    if importlib.util.find_spec(module) is None:
        pytest.skip(f"{module} is not installed")
    return [sys.executable, "-m", module]


# annotations, without external tools

def _functions(tree):
    """Yield (function node, is_method, is_static) for every function in a module tree."""
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    static = any(getattr(d, "id", None) == "staticmethod" for d in item.decorator_list)
                    yield item, True, static
        elif isinstance(node, (ast.Module, ast.FunctionDef)):
            for item in ast.iter_child_nodes(node):
                if isinstance(item, ast.FunctionDef):
                    yield item, False, False


def test_every_function_is_annotated_and_documented():
    problems = []
    for path in sorted(PACKAGE.glob("*.py")):
        for node, is_method, is_static in _functions(ast.parse(path.read_text())):
            arguments = node.args.posonlyargs + node.args.args + node.args.kwonlyargs
            if is_method and not is_static:
                arguments = arguments[1:]
            arguments += [arg for arg in (node.args.vararg, node.args.kwarg) if arg is not None]
            where = f"{path.name}:{node.lineno} {node.name}"
            problems += [f"{where}: argument '{arg.arg}' has no type" for arg in arguments if arg.annotation is None]
            if node.returns is None:
                problems.append(f"{where}: no return type")
            if ast.get_docstring(node) is None:
                problems.append(f"{where}: no docstring")
    assert not problems, "\n".join(problems)


# static checks

def test_ruff_docstrings_and_annotations():
    code, output = _run([*_module_command("ruff"), "check", "src/synced_structure"])
    assert code == 0, output


def test_pydoclint_docstrings_match_signatures():
    executable = shutil.which("pydoclint") or str(Path(sys.executable).parent / "pydoclint")
    if not Path(executable).is_file():
        pytest.skip("pydoclint is not installed")
    code, output = _run([executable, "src/synced_structure"])
    assert code == 0, output


def test_mypy_static_types():
    code, output = _run(_module_command("mypy"))
    assert code == 0, output


# runtime checks

def test_runtime_types():
    if importlib.util.find_spec("typeguard") is None:
        pytest.skip("typeguard is not installed")
    code, output = _run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                         "--typeguard-packages=synced_structure", *OFFLINE_TESTS])
    assert code == 0, output[-6000:]
