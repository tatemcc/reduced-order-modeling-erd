"""Static docstring audit for ERD capstone modules.

Usage:
    python scripts/check_docstrings.py

The check enforces:
- module docstrings,
- class docstrings,
- public function/method docstrings,
- dataclass docstrings that include an ``Attributes:`` section.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List


ROOT = Path(__file__).resolve().parents[1]
TARGET_DIRS = [
    ROOT / "erd_fipy" / "erd_fipy",
    ROOT / "erd_fipy" / "scripts",
    ROOT / "model" / "erd_model",
    ROOT / "model" / "scripts",
    ROOT / "control" / "erd_control",
    ROOT / "control" / "scripts",
]


@dataclass
class Finding:
    path: Path
    line: int
    message: str



def _is_dataclass(class_node: ast.ClassDef) -> bool:
    for dec in class_node.decorator_list:
        if isinstance(dec, ast.Name) and dec.id == "dataclass":
            return True
        if isinstance(dec, ast.Call) and isinstance(dec.func, ast.Name) and dec.func.id == "dataclass":
            return True
    return False



def _iter_py_files(dirs: Iterable[Path]) -> Iterable[Path]:
    for d in dirs:
        if not d.exists():
            continue
        for path in sorted(d.rglob("*.py")):
            yield path



def collect_findings(path: Path) -> List[Finding]:
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    findings: List[Finding] = []

    if ast.get_docstring(tree) is None:
        findings.append(Finding(path, 1, "missing module docstring"))

    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            cdoc = ast.get_docstring(node)
            if cdoc is None:
                findings.append(Finding(path, node.lineno, f"class '{node.name}' missing docstring"))
            if _is_dataclass(node) and (cdoc is None or "Attributes:" not in cdoc):
                findings.append(
                    Finding(path, node.lineno, f"dataclass '{node.name}' docstring missing 'Attributes:' section")
                )

            for child in node.body:
                if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    continue
                if child.name.startswith("__"):
                    continue
                if ast.get_docstring(child) is None:
                    findings.append(
                        Finding(path, child.lineno, f"method '{node.name}.{child.name}' missing docstring")
                    )

        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("_"):
                continue
            if ast.get_docstring(node) is None:
                findings.append(Finding(path, node.lineno, f"function '{node.name}' missing docstring"))

    return findings



def main() -> int:
    all_findings: List[Finding] = []
    for path in _iter_py_files(TARGET_DIRS):
        try:
            all_findings.extend(collect_findings(path))
        except SyntaxError as exc:
            all_findings.append(Finding(path, exc.lineno or 1, f"syntax error during parse: {exc.msg}"))

    if all_findings:
        for f in all_findings:
            print(f"{f.path}:{f.line}: {f.message}")
        return 1

    print("Docstring audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
