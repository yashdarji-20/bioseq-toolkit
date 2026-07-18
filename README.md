# 🧬 BioSeq Toolkit

> A professional DNA, RNA, and Protein sequence analysis toolkit built entirely in Python.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-264%20passing-brightgreen.svg)](tests/)
[![Coverage](https://img.shields.io/badge/coverage-99%25-brightgreen.svg)](tests/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Typed: mypy](https://img.shields.io/badge/typed-mypy-blue.svg)](https://mypy-lang.org/)

---

## Overview

**BioSeq Toolkit** is a command-line tool and Python library for bioinformatics
sequence analysis — validation, composition metrics, transcription, translation,
gene finding, sequence alignment, variant detection, and PCR primer design —
implemented from scratch in pure Python with a clean, modular, fully-tested
architecture.

It is built to demonstrate professional software-engineering practice applied to
scientific code: strict separation of concerns, comprehensive unit testing,
static type checking, custom exception hierarchies, and a real command-line
interface.

## Motivation

Most bioinformatics scripting is quick and disposable. This project takes the
opposite approach: treat a small scientific toolkit as production software.
Every function is validated, type-hinted, documented, and tested; the code is
organised into clear layers (`core`, `io`, `cli`, `utils`); and the whole thing
installs and runs as a proper command-line tool. The goal is a codebase that
reads the way professional software is built, not the way a homework script is.

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

### Infrastructure
- **FASTA I/O** — read and write multi-record FASTA files with multi-line
  sequence support
- **Command-line interface** — a `bioseq` command with `gc`, `reverse`,
  `translate`, and `stats` subcommands

## Installation

Requires Python 3.10 or newer.

```bash
git clone https://github.com/yashdarji-20/bioseq-toolkit.git
cd bioseq-toolkit

python -m venv .venv
source .venv/bin/activate         # Windows: .venv\Scripts\activate

pip install -e ".[dev]"           # installs the package plus dev tools
```

## Usage

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
from bioseq_toolkit.core.orf import find_orfs
from bioseq_toolkit.core.alignment import align
from bioseq_toolkit.core.mutation import compare_sequences
from bioseq_toolkit.core.primer import analyze_primer

gc_content("ATGGGCCC")                    # 75.0
translate("ATGAAACGT")                    # "MKR"
```

**Finding open reading frames:**

```python
orfs = find_orfs(genome, min_length=75)
for orf in orfs:
    print(orf.strand, orf.frame, orf.start, orf.protein)
```

**Aligning two sequences:**

```python
result = align("ATCG", "ACG")
print(result.formatted())
# ATCG
# | ||
# A-CG
print(result.score, result.identity)
```

**Detecting variants:**

```python
report = compare_sequences("ATCGATCGAT", "GTCGCTCGAT")
print(report.total_variants)     # 2
print(report.transitions)        # 1
print(report.transversions)      # 1
print(report.ti_tv_ratio)        # 1.0
```

**Evaluating a PCR primer:**

```python
print(analyze_primer("ACGTTGCATTGACCTGAAGC").report())
# Primer: ACGTTGCATTGACCTGAAGC
# Verdict: SUITABLE
#   Length:            20 bases
#   GC content:        50.0%
#   Melting temp:      60.0 C
#   GC clamp (3' end): yes
```

## Project Structure