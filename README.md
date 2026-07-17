The design keeps a strict separation of concerns: `core` contains pure logic
with no knowledge of files or the command line, so the same functions power
both the CLI and (in a later version) a web interface.

## Biological Background

**GC content** is the percentage of bases that are guanine or cytosine. Because
G–C pairs are held by three hydrogen bonds versus two for A–T, GC-rich DNA is
more thermally stable — a fact used in primer design and organism classification.

**Transcription** copies DNA into messenger RNA (every `T` becomes `U`), and
**translation** reads that RNA in three-base *codons*, each specifying one amino
acid, until a stop codon is reached. Together these form the *central dogma* of
molecular biology: DNA → RNA → protein.

**Melting temperature (Tm)** estimates the temperature at which a DNA duplex
separates into single strands — essential for designing PCR primers.

## Testing

The project has **110 unit tests** with **98% code coverage**. Run them with:

```bash
pytest
```

Tests cover happy paths, edge cases (empty input, whitespace, wrong alphabet),
error handling, and mathematical invariants (e.g. the reverse complement of a
reverse complement returns the original sequence).

## Code Quality

Code quality is enforced with three tools, all configured in `pyproject.toml`:

```bash
ruff check src tests    # linting
black src tests         # formatting
mypy src                # static type checking
```

## Future Improvements

- **Version 2:** ORF finder, restriction-enzyme search, motif finder, SNP
  detection, sequence alignment, primer analysis.
- **Version 3:** Streamlit web application, Docker support, PyPI packaging,
  and GitHub Actions for continuous integration.

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for details.