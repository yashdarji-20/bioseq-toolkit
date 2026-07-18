"""Streamlit web interface for BioSeq Toolkit.

Like the CLI, this module is a thin presentation layer: it collects input,
calls the tested ``core`` functions, and renders the results. It contains no
sequence-analysis logic of its own, which is why an entirely new interface
could be added without touching the library.

Run locally with::

    streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path

import streamlit as st

from bioseq_toolkit.core.alignment import align
from bioseq_toolkit.core.composition import gc_content, nucleotide_frequency
from bioseq_toolkit.core.motif import (
    RESTRICTION_ENZYMES,
    digest,
    find_motif,
    find_restriction_sites,
)
from bioseq_toolkit.core.mutation import compare_sequences
from bioseq_toolkit.core.orf import find_orfs
from bioseq_toolkit.core.primer import analyze_primer
from bioseq_toolkit.core.stats import melting_temperature, molecular_weight
from bioseq_toolkit.core.transforms import reverse_complement, transcribe
from bioseq_toolkit.core.translation import translate
from bioseq_toolkit.io.fasta import read_fasta
from bioseq_toolkit.utils.exceptions import BioSeqError

EXAMPLE_SEQUENCE = "ATGGCTAAAGGGCCCTTTAAACGTACGTAGCTAGCTAGCATCGATCGTAA"


def _read_uploaded_fasta(uploaded_file: io.BytesIO) -> str:
    """Return the first sequence from an uploaded FASTA file."""
    text = uploaded_file.getvalue().decode("utf-8")
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "upload.fasta"
        path.write_text(text, encoding="utf-8")
        records = read_fasta(path)
    return records[0].sequence


def _sequence_input(key_prefix: str) -> str:
    """Render the shared sequence input controls and return the sequence.

    Args:
        key_prefix: A unique prefix for this tab's widget keys. Streamlit
            identifies widgets by type and parameters, so identical widgets
            rendered in different tabs need distinct keys to avoid collisions.

    Returns:
        The sequence the user supplied, or an empty string.
    """
    method = st.radio(
        "Provide a sequence by",
        ["Typing or pasting", "Uploading a FASTA file"],
        horizontal=True,
        key=f"{key_prefix}_method",
    )

    if method == "Uploading a FASTA file":
        uploaded = st.file_uploader(
            "FASTA file", type=["fasta", "fa", "txt"], key=f"{key_prefix}_upload"
        )
        if uploaded is None:
            return ""
        try:
            return _read_uploaded_fasta(uploaded)
        except BioSeqError as exc:
            st.error(f"Could not read that file: {exc}")
            return ""

    return st.text_area(
        "DNA sequence",
        value=EXAMPLE_SEQUENCE,
        height=120,
        key=f"{key_prefix}_text",
    )


def render_analysis_tab() -> None:
    """Composition, statistics, and transformations for one sequence."""
    st.subheader("Sequence analysis")
    sequence = _sequence_input("analysis")

    if not sequence:
        st.info("Enter a sequence above to see its analysis.")
        return

    try:
        columns = st.columns(4)
        columns[0].metric("Length", f"{len(sequence.strip())} bp")
        columns[1].metric("GC content", f"{gc_content(sequence):.1f}%")
        columns[2].metric("Melting temp", f"{melting_temperature(sequence):.0f} C")
        columns[3].metric("Mol. weight", f"{molecular_weight(sequence):,.0f} Da")

        st.bar_chart(nucleotide_frequency(sequence))

        st.markdown("**Transformations**")
        st.code(f"Reverse complement:  {reverse_complement(sequence)}", language=None)
        st.code(f"Transcribed to RNA:  {transcribe(sequence)}", language=None)
        st.code(f"Translated protein:  {translate(sequence)}", language=None)
    except BioSeqError as exc:
        st.error(f"{exc}")


def render_orf_tab() -> None:
    """Open reading frame detection."""
    st.subheader("Open reading frames")
    sequence = _sequence_input("orf")
    min_length = st.slider(
        "Minimum ORF length (bases)", 6, 300, 30, step=3, key="orf_min_length"
    )

    if not sequence:
        st.info("Enter a sequence above to search for ORFs.")
        return

    try:
        orfs = find_orfs(sequence, min_length=min_length)
    except BioSeqError as exc:
        st.error(f"{exc}")
        return

    if not orfs:
        st.warning("No open reading frames found at this minimum length.")
        return

    st.success(f"Found {len(orfs)} open reading frame(s).")
    st.dataframe(
        [
            {
                "Strand": orf.strand,
                "Frame": orf.frame,
                "Start": orf.start,
                "End": orf.end,
                "Length": orf.length,
                "Protein": orf.protein,
            }
            for orf in orfs
        ],
        use_container_width=True,
    )


def render_alignment_tab() -> None:
    """Pairwise alignment and variant comparison."""
    st.subheader("Compare two sequences")

    left, right = st.columns(2)
    sequence_a = left.text_area(
        "Reference", value="ATCGATCGATCG", height=100, key="compare_reference"
    )
    sequence_b = right.text_area(
        "Sample", value="ATCGATTGATCG", height=100, key="compare_sample"
    )

    if not sequence_a or not sequence_b:
        st.info("Enter both sequences to compare them.")
        return

    try:
        alignment = align(sequence_a, sequence_b)
        report = compare_sequences(sequence_a, sequence_b)
    except BioSeqError as exc:
        st.error(f"{exc}")
        return

    columns = st.columns(4)
    columns[0].metric("Score", alignment.score)
    columns[1].metric("Identity", f"{alignment.identity:.1f}%")
    columns[2].metric("Variants", report.total_variants)
    columns[3].metric("Ti / Tv", f"{report.ti_tv_ratio:.2f}")

    st.markdown("**Alignment**")
    st.code(alignment.formatted(), language=None)

    if report.variants:
        st.markdown("**Variants**")
        st.dataframe(
            [
                {
                    "Position": variant.position,
                    "Reference": variant.reference,
                    "Sample": variant.sample,
                    "Type": variant.mutation_type.value,
                }
                for variant in report.variants
            ],
            use_container_width=True,
        )


def render_search_tab() -> None:
    """Motif and restriction-site searching."""
    st.subheader("Search a sequence")
    sequence = _sequence_input("search")

    if not sequence:
        st.info("Enter a sequence above to search it.")
        return

    st.markdown("**Motif search**")
    motif = st.text_input(
        "Motif (IUPAC codes allowed, e.g. GGNTTA)", value="ATG", key="search_motif"
    )
    if motif:
        try:
            matches = find_motif(sequence, motif)
            if matches:
                positions = ", ".join(str(m.start) for m in matches)
                st.success(
                    f"Found {len(matches)} occurrence(s) at positions: {positions}"
                )
            else:
                st.warning("No occurrences of that motif.")
        except (BioSeqError, ValueError) as exc:
            st.error(f"{exc}")

    st.markdown("**Restriction digest**")
    enzyme = st.selectbox("Enzyme", sorted(RESTRICTION_ENZYMES), key="search_enzyme")
    try:
        sites = find_restriction_sites(sequence, enzyme)
        fragments = digest(sequence, enzyme)
        if sites:
            st.success(
                f"{enzyme} cuts at {len(sites)} site(s), producing "
                f"{len(fragments)} fragment(s)."
            )
            st.dataframe(
                [
                    {"Fragment": i, "Length": len(f), "Sequence": f}
                    for i, f in enumerate(fragments, start=1)
                ],
                use_container_width=True,
            )
        else:
            st.warning(f"{enzyme} does not cut this sequence.")
    except (BioSeqError, ValueError) as exc:
        st.error(f"{exc}")


def render_primer_tab() -> None:
    """PCR primer evaluation."""
    st.subheader("Evaluate a PCR primer")
    primer = st.text_input(
        "Primer sequence", value="ACGTTGCATTGACCTGAAGC", key="primer_sequence"
    )

    if not primer:
        st.info("Enter a primer to evaluate it.")
        return

    try:
        analysis = analyze_primer(primer)
    except BioSeqError as exc:
        st.error(f"{exc}")
        return

    if analysis.is_suitable:
        st.success("This primer meets every design criterion.")
    else:
        st.warning("This primer has problems (listed below).")

    columns = st.columns(4)
    columns[0].metric("Length", f"{analysis.length} bp")
    columns[1].metric("GC content", f"{analysis.gc_percent:.1f}%")
    columns[2].metric("Melting temp", f"{analysis.melting_temp:.0f} C")
    columns[3].metric("GC clamp", "Yes" if analysis.has_gc_clamp else "No")

    for warning in analysis.warnings:
        st.write(f"- {warning}")


def main() -> None:
    """Compose the page."""
    st.set_page_config(page_title="BioSeq Toolkit", page_icon="🧬", layout="wide")

    st.title("🧬 BioSeq Toolkit")
    st.caption(
        "DNA sequence analysis: composition, gene finding, alignment, "
        "restriction mapping, and primer design."
    )

    analysis, orfs, comparison, search, primer = st.tabs(
        ["Analysis", "ORFs", "Compare", "Search", "Primer"]
    )

    with analysis:
        render_analysis_tab()
    with orfs:
        render_orf_tab()
    with comparison:
        render_alignment_tab()
    with search:
        render_search_tab()
    with primer:
        render_primer_tab()


if __name__ == "__main__":
    main()
