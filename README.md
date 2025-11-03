# ERD-Capstone

**Electrodeless Ring Discharge (ERD) Capstone Project**

This repository implements the full ERD modeling pipeline:

1. **`erd_fipy`** – real-time plasma and field simulation (FiPy-based data generator)
2. **`model`** – reduced-order modeling / PySINDy analysis
3. **`control`** – control-loop and plant-model integration

All subprojects share a single reproducible environment managed by [**uv**](https://docs.astral.sh/uv/).

---

## 🧩  Project layout

```
erd-capstone/
├── pyproject.toml          # uv workspace & dev tools
├── uv.lock                 # lockfile (commit this)
├── README.md
├── .python-version         # "3.10"
├── erd_fipy/               # FiPy-based data generator
│   ├── pyproject.toml
│   └── erd_fipy/
│       ├── __init__.py
│       ├── config.py
│       ├── pdes.py
│       ├── ...
│       └── scripts/
│           └── run_demo.py
├── model/                  # PySINDy / ROM
│   ├── pyproject.toml
│   └── model/
│       └── ...
└── control/                # closed-loop control
    ├── pyproject.toml
    └── control/
        └── ...
```

---

## ⚙️  Environment setup (Python 3.10 + uv)

### 0. Prerequisites

* Python ≥ 3.10 installed (or let uv install it)
* [uv](https://docs.astral.sh/uv/) installed (`pip install uv` or system package)

### 1. Clone the repository

```bash
git clone git@github.com:<org-or-user>/erd-capstone.git
cd erd-capstone
```

### 2. Pin Python version

```bash
uv python install 3.10          # one time per machine
```

### 3. Create the environment and install dependencies

```bash
uv sync --python 3.10 --group dev
```

This command:

* creates a virtual environment at `.venv/`
* installs all workspace packages (`erd_fipy`, `model`, `control`)
* installs shared dev tools (black, ruff, pytest, mypy, etc.)

### 4. Activate the environment (optional)

You can prefix commands with `uv run`, or activate `.venv` directly:

```bash
source .venv/bin/activate
```

---

## 🔬  FiPy setup

By default, FiPy is installed from PyPI and works out of the box.
However, the FiPy repo itself comes with some helpful examples, so it may be useful to clone it separately and install in editable mode inside this environment:

```bash
git clone https://github.com/usnistgov/fipy.git <directory of your choice>
uv run pip uninstall -y fipy          # remove the PyPI wheel
uv run pip install -e <directory of your choice>     # use your local checkout
```

This overrides the packaged version locally without affecting others.
CI and teammates will continue to use the stable PyPI release.

---

## ▶️  Running & testing

### Run the demo simulation

```bash
uv run python erd_fipy/scripts/run_demo.py
```

Outputs appear under `outputs/demo_run.h5`.

### Run tests (once you add them)

```bash
uv run pytest -q
```

### Lint and format

```bash
uv run ruff check . --fix
uv run black .
```

---

## 🧑‍💻  Development workflow (Git + uv)

1. **Stay up-to-date**

   ```bash
   git pull origin main
   uv sync --frozen --group dev
   ```

2. **Create a feature branch**

   ```bash
   git checkout -b feat/<short-description>
   ```

3. **Develop**

   * Work in VS Code (select `.venv/bin/python` or kernel `erd-capstone`)
   * Run/lint/format/test using `uv run ...`

4. **Commit & push**

   ```bash
   git add -A
   git commit -m "Implement new plasma BC in erd_fipy"
   git push -u origin feat/<short-description>
   ```

5. **Open a Pull Request**

   * Target `main`
   * CI (ruff + black + pytest) runs automatically
   * Review → merge

6. **After merge**

   ```bash
   git checkout main
   git pull
   uv sync --frozen --group dev
   ```

---

## 🧹  Pre-commit hooks

Install once:

```bash
uv run pre-commit install
```

Then every `git commit` automatically runs:

* ruff (lint)
* black (format)
* whitespace fixes

Run manually anytime:

```bash
uv run pre-commit run --all-files
```

---

## 📂  Data management

* **Do not commit** large `.h5` output files.
  Use scratch storage (e.g., `/scratch/$USER/erd-capstone/outputs/<run-id>.h5`).
* Commit only small test fixtures under `model/tests/fixtures/`.
* Document dataset schema (HDF5 group / dataset names) in `erd_fipy/README.md`.

---

## 🧠  Quick reference

| Task                 | Command                                       |
| -------------------- | --------------------------------------------- |
| Install everything   | `uv sync --group dev`                         |
| Run code             | `uv run python ...`                           |
| Run tests            | `uv run pytest -q`                            |
| Lint / format        | `uv run ruff check . --fix && uv run black .` |
| Update deps          | `uv add <pkg>` → `uv lock` → commit `uv.lock` |
| Recreate environment | `rm -rf .venv && uv sync`                     |

---

## 📜  License

This project is released under the [MIT License](LICENSE).
