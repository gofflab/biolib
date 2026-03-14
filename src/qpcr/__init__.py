"""
qpcr — Quantitative PCR data analysis utilities.

Modules:
    abi           ABI instrument file parsing and data loading
    qpcrAnalysis  ddCt analysis and qPCR workflows (requires rpy2)
    MinerMethod   Miner method for PCR efficiency estimation
    util          Utility functions for qPCR data processing
"""

__version__ = "0.2.0"

from . import abi
from . import MinerMethod
from . import qpcrAnalysis
from . import util

__all__ = ["abi", "MinerMethod", "qpcrAnalysis", "util"]
