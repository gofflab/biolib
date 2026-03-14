# Contributing to biolib

## Development Setup

1. Clone the repository:

   ```bash
   git clone https://github.com/gofflab/biolib.git
   cd biolib
   ```

2. Create a virtual environment and install in editable mode with dev dependencies:

   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev]"
   ```

3. Install pre-commit hooks:

   ```bash
   pip install pre-commit
   pre-commit install
   ```

## Running Tests

```bash
pytest
```

With coverage report:

```bash
pytest --cov=src --cov-report=html
open htmlcov/index.html
```

## Code Style

This project uses [ruff](https://docs.astral.sh/ruff/) for linting and formatting.

Check for issues:

```bash
ruff check src/
```

Auto-fix issues:

```bash
ruff check --fix src/
```

Format code:

```bash
ruff format src/
```

## Branch Naming

- Features: `feature/<short-description>`
- Bug fixes: `fix/<short-description>`
- Automated branches: `claude/<description>-<id>`

## Commit Messages

Use clear, imperative commit messages:

- `Add GTFlib support for GFF3 format`
- `Fix off-by-one error in intervallib.overlap()`
- `Upgrade seqlib to Python 3.12`

## Adding a New Module

1. Create the module in `src/seqlib/` or `src/qpcr/`
2. Add it to `__all__` in the corresponding `__init__.py`
3. Add smoke tests in `tests/test_seqlib.py` or `tests/test_qpcr.py`
4. Document it in `README.md` module table
5. Note the addition in `CHANGELOG.md` under `[Unreleased]`

## Dependency Notes

- **pygr**: Legacy genome database library — unmaintained and Python 2 only.
  `seqlib.genomelib` and `seqlib.pygrlib` depend on it and are non-functional
  in Python 3. Do not add new code using `pygr`.

- **rasmus**: Legacy utility library — not Python 3 compatible.
  All `rasmus` references have been replaced with local implementations or removed.

- **rpy2**: Optional dependency for R integration. Required by `qpcr.qpcrAnalysis`
  for ddCt analysis. Not required for pure-Python functionality.
