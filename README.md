# biolib

Personal computational biology utility library for sequence analysis and qPCR data
processing, built for Python 3.12+.

## Installation

```bash
pip install -e ".[dev]"
```

### Requirements

- Python >= 3.12
- numpy >= 1.26
- scipy >= 1.12
- pysam >= 0.22
- rpy2 >= 3.5 (required for R-based qPCR analysis and enrichment functions)

## Modules

### `seqlib` — Sequence Analysis Utilities

A broad collection of bioinformatics tools for next-generation sequencing analysis.

| Module                  | Description                                      |
|-------------------------|--------------------------------------------------|
| `seqlib.stats`          | Statistical functions for genomic data           |
| `seqlib.util`           | General-purpose utility functions                |
| `seqlib.seqlib`         | Core sequence manipulation                       |
| `seqlib.seqstats`       | Sequence-level statistics                        |
| `seqlib.intervallib`    | Genomic interval operations                      |
| `seqlib.mySam`          | SAM/BAM file handling                            |
| `seqlib.GTFlib`         | GTF/GFF annotation parsing                       |
| `seqlib.algorithms`     | Common bioinformatics algorithms                 |
| `seqlib.prob`           | Probability distributions                        |
| `seqlib.JensenShannon`  | Jensen-Shannon divergence                        |
| `seqlib.Alignment`      | Sequence alignment utilities                     |
| `seqlib.Chip`           | ChIP-seq analysis tools                          |
| `seqlib.clustering`     | Clustering algorithms                            |
| `seqlib.converters`     | Format conversion utilities                      |
| `seqlib.bowtie`         | Bowtie aligner wrappers                          |
| `seqlib.bwa`            | BWA aligner wrappers                             |
| `seqlib.LSFlib`         | LSF cluster job submission                       |
| `seqlib.QCtools`        | Quality control tools                            |
| `seqlib.RIPDiff`        | RIP-seq differential analysis                    |
| `seqlib.continuousData` | Continuous data representation and operations    |
| `seqlib.blockIt`        | Block-based data iteration                       |
| `seqlib.misc`           | Miscellaneous helper functions                   |

### `qpcr` — qPCR Analysis

Tools for quantitative PCR data processing and analysis.

| Module               | Description                                  |
|----------------------|----------------------------------------------|
| `qpcr.abi`           | ABI instrument file parsing                  |
| `qpcr.qpcrAnalysis`  | ddCt analysis and qPCR workflows             |
| `qpcr.MinerMethod`   | Miner method for PCR efficiency estimation   |
| `qpcr.util`          | Utility functions for qPCR data              |

## Usage Examples

### Parse a GTF annotation file

```python
from seqlib import GTFlib

gtf = GTFlib.GTFReader("annotation.gtf")
for gene in gtf:
    print(gene.gene_id, gene.chrom, gene.start, gene.end)
```

### Compute Jensen-Shannon divergence

```python
from seqlib.JensenShannon import JS_divergence

p = [0.25, 0.25, 0.25, 0.25]
q = [0.50, 0.50, 0.00, 0.00]
divergence = JS_divergence(p, q)
print(divergence)
```

### Work with genomic intervals

```python
from seqlib import intervallib

interval = intervallib.Interval("chr1", 1000, 2000, strand="+")
print(interval.length())
```

### Load ABI qPCR results

```python
from qpcr import abi

data = abi.parseABIResults("results.txt", "cycleData.txt")
```

### Run ddCt qPCR analysis

```python
from qpcr import qpcrAnalysis

results = qpcrAnalysis.ddCtAnalysis(
    data_file="results.txt",
    endogenous_control="GapDH",
    reference_sample="control"
)
```

## Development

### Setup

```bash
git clone https://github.com/gofflab/biolib.git
cd biolib
pip install -e ".[dev]"
```

### Running Tests

```bash
pytest
```

With coverage:

```bash
pytest --cov=src --cov-report=html
```

### Linting and Formatting

```bash
# Check for issues
ruff check src/

# Auto-fix issues
ruff check --fix src/

# Format code
ruff format src/
```

### Pre-commit Hooks

```bash
pip install pre-commit
pre-commit install
```

## Project Structure

```
biolib/
├── src/
│   ├── qpcr/           # qPCR analysis modules
│   └── seqlib/         # Sequence analysis modules
├── tests/              # Test suite
├── pyproject.toml      # Package configuration
└── requirements.txt    # Pinned dependencies
```

## License

MIT
