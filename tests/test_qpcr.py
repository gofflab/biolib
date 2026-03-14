"""Smoke tests for the qpcr package."""

import pytest


def test_abi_import():
    from qpcr import abi
    assert abi is not None


def test_miner_import():
    from qpcr import MinerMethod
    assert MinerMethod is not None


def test_qpcr_analysis_import():
    from qpcr import qpcrAnalysis
    assert qpcrAnalysis is not None


def test_util_import():
    from qpcr import util
    assert util is not None


def test_package_version():
    import qpcr
    assert hasattr(qpcr, "__version__")
    assert qpcr.__version__ == "0.2.0"
