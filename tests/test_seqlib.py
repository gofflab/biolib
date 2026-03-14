"""Smoke tests for the seqlib package."""

import pytest


def test_package_version():
    import seqlib
    assert hasattr(seqlib, "__version__")
    assert seqlib.__version__ == "0.2.0"


def test_stats_import():
    from seqlib import stats
    assert stats is not None


def test_util_import():
    from seqlib import util
    assert util is not None


def test_algorithms_import():
    from seqlib import algorithms
    assert algorithms is not None


def test_prob_import():
    from seqlib import prob
    assert prob is not None


def test_gtflib_import():
    from seqlib import GTFlib
    assert GTFlib is not None


def test_intervallib_import():
    from seqlib import intervallib
    assert intervallib is not None


def test_jensen_shannon_import():
    from seqlib import JensenShannon
    assert JensenShannon is not None


def test_seqstats_import():
    from seqlib import seqstats
    assert seqstats is not None


def test_mysam_import():
    from seqlib import mySam
    assert mySam is not None


def test_misc_import():
    from seqlib import misc
    assert misc is not None


def test_converters_import():
    from seqlib import converters
    assert converters is not None


def test_clustering_import():
    from seqlib import clustering
    assert clustering is not None


def test_blockIt_import():
    from seqlib import blockIt
    assert blockIt is not None


def test_continuous_data_import():
    from seqlib import continuousData
    assert continuousData is not None


def test_alignment_import():
    from seqlib import Alignment
    assert Alignment is not None


def test_chip_import():
    from seqlib import Chip
    assert Chip is not None


def test_lsflib_import():
    from seqlib import LSFlib
    assert LSFlib is not None


def test_qctools_import():
    from seqlib import QCtools
    assert QCtools is not None


def test_ripdiff_import():
    from seqlib import RIPDiff
    assert RIPDiff is not None


def test_bowtie_import():
    from seqlib import bowtie
    assert bowtie is not None


def test_bwa_import():
    from seqlib import bwa
    assert bwa is not None
