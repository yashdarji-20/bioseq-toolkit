# 🧬 BioSeq Toolkit

> A professional DNA, RNA, and Protein sequence analysis toolkit built entirely in Python.

[![CI](https://github.com/yashdarji-20/bioseq-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/yashdarji-20/bioseq-toolkit/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-264%20passing-brightgreen.svg)](tests/)
[![Coverage](https://img.shields.io/badge/coverage-99%25-brightgreen.svg)](tests/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Live Demo](https://img.shields.io/badge/live%20demo-streamlit-red.svg)](https://bioseq-toolkit.streamlit.app/)

**🔗 [Try it live: bioseq-toolkit.streamlit.app](https://bioseq-toolkit.streamlit.app/)** — no install required.

---

## Overview

**BioSeq Toolkit** is a command-line tool, Python library, and web application for
bioinformatics sequence analysis — validation, composition metrics, transcription,
translation, gene finding, sequence alignment, variant detection, and PCR primer
design — implemented from scratch in pure Python with a clean, modular, fully-tested
architecture.

It is built to demonstrate professional software-engineering practice applied to
scientific code: strict separation of concerns, comprehensive unit testing, static
type checking, custom exception hierarchies, continuous integration, and two
independent interfaces sharing one tested core.

## Motivation

Most bioinformatics scripting is quick and disposable. This project takes the
opposite approach: treat a small scientific toolkit as production software.
Every function is validated, type-hinted, documented, and tested; the code is
organised into clear layers (`core`, `io`, `cli`, `utils`); and the whole thing
installs and runs as a proper command-line tool and a deployed web app. The goal
is a codebase that reads the way professional software is built, not the way a
homework script is.

## Features

### Sequence fundamentals
- **Validation** — check whether a sequence is well-formed DNA, RNA, or protein
- **Composition** — GC content, AT content, per-base nucleotide frequency
- **Transformation** — complement, reverse complement, DNA→RNA transcription
- **Translation** — DNA→protein with configurable stop-codon handling
- **Statistics** — length, molecular weight, melting temperature (Wallace rule)

### Analysis
- **ORF Finder** — six-frame open reading frame detection across both strands
- **Sequence Alignment** — global alignment via the Needleman-Wunsch dynamic
  programming algorithm, with traceback and percent identity
- **Motif Finder** — pattern search supporting IUPAC ambiguity codes
  (`N`, `R`, `Y`, `W`, `S`, …), including overlapping occurrences
- **Restriction Enzyme Search** — recognition-site detection and in-silico
  digestion for common enzymes (EcoRI, BamHI, HindIII, NotI, PstI, SmaI, …)
- **SNP Detection** — variant calling between two sequences, with
  transition/transversion classification and Ti/Tv ratio; automatically aligns
  sequences of differing length so indels are not mistaken for substitutions
- **Primer Analysis** — evaluates PCR primers against standard design criteria
  (length, GC content, Tm, 3′ GC clamp, self-complementarity, homopolymer runs)
  and reports actionable warnings

### Interfaces and infrastructure
- **FASTA I/O** — read and write multi-record FASTA files with multi-line
  sequence support
- **Command-line interface** — a `bioseq` command with `gc`, `reverse`,
  `translate`, and `stats` subcommands
- **Web application** — a [live Streamlit interface](https://bioseq-toolkit.streamlit.app/)
  covering every analysis feature
- **Continuous integration** — GitHub Actions runs the full test suite, linting,
  formatting, and type checks across Python 3.10–3.13 on every push

## Installation

Requires Python 3.10 or newer.

```bash
git clone https://github.com/yashdarji-20/bioseq-toolkit.git
cd bioseq-toolkit

python -m venv .venv
source .venv/bin/activate         # Windows: .venv\Scripts\activate

pip install -e ".[dev]"           # package plus development tools
```

## Usage

### Web application

The easiest way to try the toolkit is the **[live web app](https://bioseq-toolkit.streamlit.app/)** —
no installation needed. It provides tabs for sequence composition, ORF detection,
pairwise comparison, motif and restriction-site search, and primer evaluation.

The app imports the same `core` functions as the CLI; no analysis logic is
duplicated between the two interfaces.

To run it locally:

```bash
pip install -e ".[app]"
streamlit run app/streamlit_app.py
```

### Command line

Given a FASTA file `data/sample.fasta`:

```
>seq1 test sequence
ATGGCTTAAGGGCCC
>seq2 another
ATGAAACGTACGT
```

```bash
$ bioseq gc data/sample.fasta
seq1 test sequence    GC=60.00%
seq2 another          GC=38.46%

$ bioseq stats data/sample.fasta
seq1 test sequence    length=15    GC=60.00%    A=3 C=4 G=5 T=3
seq2 another          length=13    GC=38.46%    A=5 C=2 G=3 T=3

$ bioseq translate data/sample.fasta
seq1 test sequence    MA
seq2 another          MKRT

$ bioseq reverse data/sample.fasta
seq1 test sequence    GGGCCCTTAAGCCAT
seq2 another          ACGTACGTTTCAT
```

### As a library

```python
from bioseq_toolkit.core.composition import gc_content
from bioseq_toolkit.core.translation import translate

gc_content("ATGGGCCC")                    # 75.0
translate("ATGAAACGT")                    # "MKR"
```

**Finding open reading frames:**

```python
from bioseq_toolkit.core.orf import find_orfs

orfs = find_orfs(genome, min_length=75)
for orf in orfs:
    print(orf.strand, orf.frame, orf.start, orf.protein)
```

**Aligning two sequences:**

```python
from bioseq_toolkit.core.alignment import align

result = align("ATCG", "ACG")
print(result.formatted())
# ATCG
# | ||
# A-CG
print(result.score, result.identity)
```

**Detecting variants:**

```python
from bioseq_toolkit.core.mutation import compare_sequences

report = compare_sequences("ATCGATCGAT", "GTCGCTCGAT")
print(report.total_variants)     # 2
print(report.transitions)        # 1
print(report.transversions)      # 1
print(report.ti_tv_ratio)        # 1.0
```

**Evaluating a PCR primer:**

```python
from bioseq_toolkit.core.primer import analyze_primer

print(analyze_primer("ACGTTGCATTGACCTGAAGC").report())
# Primer: ACGTTGCATTGACCTGAAGC
# Verdict: SUITABLE
#   Length:            20 bases
#   GC content:        50.0%
#   Melting temp:      60.0 C
#   GC clamp (3' end): yes
```

## Project Structure

```
bioseq-toolkit/
├── src/bioseq_toolkit/
│   ├── core/
│   │   ├── validation.py    # DNA / RNA / protein validation
│   │   ├── composition.py   # GC, AT, nucleotide frequency
│   │   ├── transforms.py    # complement, reverse complement, transcription
│   │   ├── translation.py   # DNA → protein
│   │   ├── stats.py         # molecular weight, Tm, summary statistics
│   │   ├── orf.py           # six-frame ORF detection
│   │   ├── alignment.py     # Needleman-Wunsch global alignment
│   │   ├── motif.py         # motif and restriction-site search
│   │   ├── mutation.py      # SNP detection and classification
│   │   └── primer.py        # PCR primer analysis
│   ├── io/fasta.py          # FASTA reading and writing
│   ├── cli/main.py          # argparse command-line interface
│   └── utils/               # constants, exceptions, logging
├── app/streamlit_app.py     # Streamlit web interface
├── tests/                   # unit tests mirroring the package layout
├── data/                    # sample FASTA files
├── .github/workflows/ci.yml # continuous integration
├── pyproject.toml           # packaging, dependencies, tool configuration
├── LICENSE
└── README.md
```

The architecture keeps a strict separation of concerns: `core` contains pure
logic with no knowledge of files, terminals, or browsers. Both the CLI and the
web app are thin adapters over it, which is why an entirely new interface could
be added without touching the library. Higher-level features are built by
*composing* lower-level ones — variant detection calls the alignment module,
primer analysis calls the composition and stats modules — rather than
duplicating logic.

## Biological Background

**GC content** is the percentage of bases that are guanine or cytosine. Because
G–C pairs are held by three hydrogen bonds versus two for A–T, GC-rich DNA is
more thermally stable — a fact used in primer design and organism classification.

**Transcription** copies DNA into messenger RNA (every `T` becomes `U`), and
**translation** reads that RNA in three-base *codons*, each specifying one amino
acid, until a stop codon is reached. Together these form the *central dogma* of
molecular biology: DNA → RNA → protein.

**Open reading frames** are stretches beginning with a start codon (`ATG`) and
ending at an in-frame stop codon — candidate protein-coding genes. Because genes
may lie in any of three reading frames on either strand, gene finding requires a
six-frame search.

**Sequence alignment** measures similarity between sequences by inserting gaps
to maximise matching. The Needleman-Wunsch algorithm solves this optimally with
dynamic programming in O(m×n) time, where a naive search would be exponential.

**SNPs** are single-base differences between sequences. Substitutions are
classified as *transitions* (within purines or within pyrimidines) or
*transversions* (across the two classes); transitions occur roughly twice as
often in nature, making the Ti/Tv ratio a standard quality metric.

## Testing

The project has **264 unit tests** with **99% code coverage**:

```bash
pytest
```

Tests cover happy paths, boundary conditions, error handling, and mathematical
invariants — for example, that the reverse complement of a reverse complement
returns the original sequence, that stripping gaps from an alignment recovers
the input sequences, and that restriction digest fragments rejoin into the
original DNA.

## Code Quality

Quality is enforced by tooling rather than convention, all configured in
`pyproject.toml` and run automatically by CI on every push:

```bash
ruff check src tests    # linting
black --check src tests # formatting
mypy src                # static type checking
```

## Roadmap

- ✅ **Version 1.0** — validation, composition, transforms, translation,
  statistics, FASTA I/O, CLI
- ✅ **Version 2.0** — ORF finder, sequence alignment, motif and restriction
  search, SNP detection, primer analysis
- ✅ **Version 3.0** — [Streamlit web app](https://bioseq-toolkit.streamlit.app/),
  GitHub Actions continuous integration
- ⏳ **Planned** — Docker support, PyPI packaging

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for details.