"""
seqlib — Computational biology sequence analysis utilities.

This package provides tools for:
- Sequence manipulation and analysis
- Genomic interval operations
- SAM/BAM file processing
- GTF/GFF annotation parsing
- Statistical analysis of sequencing data
- Alignment tool wrappers (Bowtie, BWA)
- ChIP-seq and RIP-seq analysis

Note: Some legacy modules (genomelib, pygrlib) require the unmaintained
'pygr' library and must be imported explicitly if needed.
"""

__version__ = "0.2.0"

__all__ = [
    "algorithms",
    "Alignment",
    "blockIt",
    "bowtie",
    "bwa",
    "Chip",
    "clustering",
    "continuousData",
    "converters",
    "GTFlib",
    "intervallib",
    "JensenShannon",
    "LSFlib",
    "misc",
    "mySam",
    "prob",
    "QCtools",
    "RIPDiff",
    "seqlib",
    "seqstats",
    "stats",
    "util",
]
