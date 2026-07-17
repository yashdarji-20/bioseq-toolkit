# 🧬 BioSeq Toolkit

> A professional DNA, RNA, and Protein sequence analysis toolkit built entirely in Python.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-110%20passing-brightgreen.svg)](tests/)
[![Coverage](https://img.shields.io/badge/coverage-98%25-brightgreen.svg)](tests/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## Overview

**BioSeq Toolkit** is a command-line tool and Python library for common
bioinformatics sequence-analysis tasks — sequence validation, GC content,
transcription, translation, FASTA parsing, and more — implemented from scratch
in pure Python with a clean, modular, fully-tested architecture.

It is designed to demonstrate professional software-engineering practices:
separation of concerns, comprehensive unit testing, type safety, custom
exceptions, and a proper command-line interface.

## Motivation

Most bioinformatics scripting is quick and disposable. This project takes the
opposite approach: treat a small scientific toolkit as production software.
Every function is validated, type-hinted, documented, and tested; the code is
organized into clear layers (`core`, `io`, `cli`, `utils`); and the whole thing
installs and runs as a real command-line tool. The goal is a codebase that
reads the way professional software is built, not the way a homework script is.

## Features

**Validation** — check whether a sequence is well-formed DNA, RNA, or protein.

**Composition** — GC content, AT content, and per-base nucleotide frequency.

**Transformation** — complement, reverse complement, and DNA→RNA transcription.

**Translation** — DNA→protein translation with configurable stop-codon handling.

**Statistics** — length, molecular weight, and melting temperature (Wallace rule),
bundled into a typed result.

**FASTA I/O** — read and write multi-record FASTA files with multi-line sequence
support.

**Command-line interface** — a `bioseq` command with `gc`, `reverse`,
`translate`, and `stats` subcommands.

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