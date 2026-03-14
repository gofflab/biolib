# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `pyproject.toml` with modern setuptools packaging configuration
- `requirements.txt` with pinned dependency ranges
- `tests/` directory with smoke tests for `qpcr` and `seqlib` modules
- GitHub Actions CI workflow for linting and testing
- `.pre-commit-config.yaml` with ruff and pre-commit-hooks
- `CHANGELOG.md` and `CONTRIBUTING.md`
- `ruff`, `pytest`, and `black` configuration in `pyproject.toml`
- `__version__` attribute to both `qpcr` and `seqlib` packages
- `__all__` export list to `seqlib/__init__.py`

### Changed
- Upgraded entire codebase from Python 2 to Python 3.12
- Replaced `seqlib/__init__.py` SHRiMP pipeline stub with proper package docstring and exports
- Expanded `qpcr/__init__.py` to expose all submodules (`abi`, `MinerMethod`, `qpcrAnalysis`, `util`)
- Removed dead `rasmus` library imports from `seqlib/util.py` (were already silently failing)
- Wrapped legacy `pygr` imports in `genomelib.py` and `pygrlib.py` with `try/except ImportError`
- Replaced `import sequencelib` with relative import in `genomelib.py`

### Deprecated
- `seqlib.genomelib` — requires the unmaintained `pygr` library; use `pysam` or `pybedtools` instead
- `seqlib.pygrlib` — experimental scratch file depending on `pygr`; not suitable for production use

## [0.2.0] — Python 3.12 upgrade

### Changed
- Full Python 2 → Python 3.12 migration across all modules
- Updated `print` statements to `print()` functions
- Modernised `dict.keys()`/`values()`/`items()` usage
- Fixed exception syntax (`except X as e`)
- Updated `urllib`/`urllib2` imports for Python 3
- Fixed integer division and string handling throughout

## [0.1.0] — Initial release

- Personal compbio utility library for sequence analysis and qPCR
