"""Codon Optimiser - Streamlit web app.

Reverse-translates a protein sequence into a DNA sequence optimised for a
target organism and expression level, avoiding user-specified restriction
sites and undesirable motifs.

Run locally:  streamlit run app.py
"""

from __future__ import annotations

import csv
import io
import json
import sys
from pathlib import Path

import streamlit as st

# make the codonopt package importable when running from the app/ folder
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from codonopt import (  # noqa: E402
    PRESET_LABELS,
    OptimizationConfig,
    cai,
    gc_content,
    get_table,
    load_preset_tables,
    n_end_flag,
    optimize,
    parse_custom_table,
    validate_table,
)
from codonopt.codon_tables import table_gc3  # noqa: E402
from codonopt.translate import back_translate, clean_protein  # noqa: E402

st.set_page_config(page_title="Codon Optimiser", page_icon="🧬", layout="wide")

# ------------------------------------------------------------------ helpers

COMMON_RESTRICTION_SITES = [
    "GAATTC",  # EcoRI
    "GGATCC",  # BamHI
    "AAGCTT",  # HindIII
    "CTGCAG",  # PstI
    "GTCGAC",  # SalI
    "GCGGCCGC",  # NotI
    "ACCGGT",  # AgeI
    "TCTAGA",  # XbaI
    "ACTAGT",  # SpeI
    "GGTACC",  # KpnI
    "GAGCTC",  # SacI
    "CCGCGG",  # MluI... (SacII)
    "GCTAGC",  # NheI
    "CATATG",  # NdeI
    "CCCATGG",  # NcoI
]


def to_fasta(name: str, seq: str, width: int = 60) -> str:
    lines = [f">{name}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i:i + width])
    return "\n".join(lines)


# ------------------------------------------------------------------- sidebar

st.title("Codon Optimiser")
st.markdown(
    "Reverse-translate a **protein sequence** into a DNA sequence optimised for "
    "your target organism and desired expression level."
)

with st.sidebar:
    st.header("Parameters")

    organism = st.selectbox(
        "Target organism",
        options=list(PRESET_LABELS.values()),
        index=0,
        help="Codon usage tables derived from the Kazusa CUTG database.",
    )
    use_custom_table = st.checkbox("Use a custom codon table instead")

    custom_table_text = ""
    if use_custom_table:
        custom_table_text = st.text_area(
            "Paste custom codon table (Kazusa CUTG format)",
            height=260,
            placeholder=(
                "UUU 26.1(170666) UCU 23.5(153557) UAU 18.8(122728) UGU 8.1( 52903)\n"
                "UUC 18.4(120510) UCC 14.2( 92923) UAC 14.8( 96596) UGC 4.8( 31095)\n"
                "UUA 26.2(170884) UCA 18.7(122028) UAA 1.1( 6913) UGA 0.7( 4447)\n"
                "UUG 27.2(177573) UCG 8.6( 55951) UAG 0.5( 3312) UGG 10.4( 67789)\n"
                "\n...(all 64 codons, as shown in the codon usage database)..."
            ),
            help=(
                "Copy the table straight from the codon usage database "
                "(www.kazusa.or.jp/codon/): four codon columns per block, "
                "each entry as CODON frequency(count). All 64 codons are required."
            ),
        )

    expression = st.radio(
        "Desired expression level",
        options=["high", "medium", "low"],
        index=0,
        horizontal=True,
        help=(
            "high = maximise codon adaptation (CAI); medium = moderate CAI; "
            "low = deliberately enrich rare codons (e.g. to slow translation "
            "for folding-sensitive proteins)."
        ),
    )

    st.subheader("Restriction sites to avoid")
    selected_sites = st.multiselect("Common sites", COMMON_RESTRICTION_SITES)
    custom_sites = st.text_input(
        "Additional sites (comma-separated, IUPAC codes allowed)",
        placeholder="e.g. GAATTC, GGTCTC",
    )

    st.subheader("Sequence hygiene checks")
    use_builtin_motifs = st.toggle("Screen built-in motifs", value=True,
                                   help="PolyA signals, TATA boxes, RBS, Chi sites, cryptic splice motifs.")
    avoid_terminators = st.toggle("Avoid rho-independent terminators", value=True,
                                  help="Relevant for bacterial expression.")
    user_motifs = st.text_input(
        "Custom motifs to avoid (comma-separated, IUPAC)",
        placeholder="e.g. GCTGGTGG, AATAAA",
    )

    seed = st.number_input("Random seed (reproducibility)", min_value=0, value=42, step=1)

# ---------------------------------------------------------------- main input

st.subheader("Protein sequence")
seq_name = st.text_input(
    "Sequence name",
    value="my_sequence",
    help="Used in the downloaded FASTA header and the CSV report.",
)
protein_input = st.text_area(
    "Paste a protein sequence (raw or FASTA)",
    height=160,
    placeholder="MKVLAAALLA...\nor\n>my_protein\nMKVLAAALLA...",
)
uploaded_fasta = st.file_uploader("…or upload a protein FASTA", type=["fa", "fasta", "txt"])
if uploaded_fasta is not None:
    protein_input = uploaded_fasta.read().decode("utf-8", errors="replace")

if st.button("Optimise", type="primary"):
    # ---- resolve codon table
    table = None
    if use_custom_table:
        if not custom_table_text.strip():
            st.error("Custom table selected but empty.")
            st.stop()
        try:
            table = parse_custom_table(custom_table_text)
            problems = validate_table(table)
            if problems:
                st.error("Custom table problems: " + "; ".join(problems))
                st.stop()
        except (ValueError, json.JSONDecodeError) as exc:
            st.error(f"Could not parse custom codon table: {exc}")
            st.stop()
        table_label = "custom table"
    else:
        table = get_table(organism)
        table_label = organism

    # ---- resolve sites / motifs
    sites = list(selected_sites)
    sites += [s.strip().upper() for s in custom_sites.split(",") if s.strip()]
    motifs = [m.strip().upper() for m in user_motifs.split(",") if m.strip()]

    # ---- validate protein
    try:
        protein = clean_protein(protein_input)
    except ValueError as exc:
        st.error(f"Invalid protein sequence: {exc}")
        st.stop()

    # ---- run optimisation
    cfg = OptimizationConfig(
        expression_level=expression,
        avoid_sites=sites,
        avoid_motifs=motifs,
        use_builtin_motifs=use_builtin_motifs,
        avoid_terminators=avoid_terminators,
        seed=int(seed),
    )
    with st.spinner("Optimising codons (simulated annealing)…"):
        result = optimize(protein, table, cfg)

    # ---- outputs
    st.success(
        f"Optimised {len(protein)} residues → {len(result.dna)} nt for {table_label} "
        f"({expression} expression)."
    )

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("CAI", f"{result.cai:.3f}")
    m2.metric("GC content", f"{result.gc * 100:.1f}%")
    m3.metric("GC3", f"{result.gc3 * 100:.1f}%")
    m4.metric("5' hairpin score", f"{result.hairpin:.2f}")

    st.subheader("Optimised DNA sequence")
    st.code(result.dna, language=None)

    fasta_buf = io.StringIO(to_fasta(seq_name or "optimized_sequence", result.dna))
    st.download_button("Download FASTA", fasta_buf.getvalue(),
                       file_name="optimized_sequence.fasta", mime="text/plain")
    st.download_button("Download raw DNA (.txt)", result.dna,
                       file_name="optimized_sequence.txt", mime="text/plain")

    # ---- CSV design report
    report_buf = io.StringIO()
    writer = csv.writer(report_buf)
    writer.writerow(["field", "value"])
    writer.writerow(["sequence_name", seq_name])
    writer.writerow(["protein_sequence", protein])
    writer.writerow(["target_organism", table_label])
    writer.writerow(["expression_level", expression])
    writer.writerow(["restriction_sites_to_avoid", "; ".join(sites) if sites else "none"])
    writer.writerow(["screen_builtin_motifs", "yes" if use_builtin_motifs else "no"])
    writer.writerow(["avoid_rho_independent_terminators",
                     "yes" if avoid_terminators else "no"])
    writer.writerow(["custom_motifs_to_avoid", "; ".join(motifs) if motifs else "none"])
    writer.writerow(["random_seed", int(seed)])
    writer.writerow(["optimised_dna_sequence", result.dna])
    writer.writerow(["CAI", f"{result.cai:.4f}"])
    writer.writerow(["GC_content", f"{result.gc:.4f}"])
    writer.writerow(["GC3", f"{result.gc3:.4f}"])
    writer.writerow(["5prime_hairpin_score", f"{result.hairpin:.4f}"])
    st.download_button("Download report", report_buf.getvalue(),
                       file_name=f"{(seq_name or 'sequence')}_report.csv",
                       mime="text/csv")

    # ---- QC report
    st.subheader("Quality report")
    n_sites = sum(len(v) for v in result.sites_found.values())
    n_motifs = sum(len(v) for v in result.motifs_found.values())
    col_a, col_b = st.columns(2)

    with col_a:
        st.markdown("**Restriction sites avoided**")
        if sites:
            if n_sites == 0:
                st.markdown("✅ None of the specified sites occur (either strand).")
            else:
                for name, pos in result.sites_found.items():
                    if pos:
                        st.markdown(f"⚠️ {name}: {len(pos)} occurrence(s)")
        else:
            st.markdown("No sites specified.")

        st.markdown("**Undesirable motifs**")
        if n_motifs == 0:
            st.markdown("✅ No built-in or custom motifs found (either strand).")
        else:
            for name, pos in result.motifs_found.items():
                if pos:
                    st.markdown(f"⚠️ {name}: {len(pos)} occurrence(s)")

        if result.homopolymers:
            st.markdown(f"⚠️ Homopolymer runs ≥5: {len(result.homopolymers)}")
        else:
            st.markdown("✅ No homopolymer runs ≥5.")

        if avoid_terminators:
            if result.terminators:
                st.markdown(f"⚠️ Possible rho-independent terminators: {len(result.terminators)}")
            else:
                st.markdown("✅ No rho-independent terminator motifs.")

    with col_b:
        st.markdown("**Back-translation check**")
        st.markdown("✅ DNA translates exactly to the input protein."
                    if back_translate(result.dna) == protein else "❌ mismatch")
        n_end = n_end_flag(protein)
        if n_end:
            st.markdown(f"ℹ️ {n_end}")
        st.markdown("**Codon usage table (target organism)**")
        st.caption(f"Usage-weighted GC3 of the table: {table_gc3(table) * 100:.1f}%")
        with st.expander("Show table"):
            st.dataframe({aa: codons for aa, codons in table.items()},
                         use_container_width=True)

    with st.expander("What do these metrics mean?"):
        st.markdown(
            """
- **CAI (Codon Adaptation Index)**: 1.0 = every codon is the most-used synonym in the target organism.
  High-expression designs target CAI ≈ 1; medium/low designs intentionally use less-adapted codons.
- **GC content / GC3**: fraction of G/C over all positions / third codon positions. Extreme values can
  hurt synthesis, transcription, and mRNA stability.
- **5' hairpin score**: heuristic stability of structure in the first 50 nt; lower is better for
  translation initiation.
- **Motif screening**: polyA signals, TATA boxes, Shine-Dalgarno motifs, Chi sites, cryptic splice
  motifs, and (for bacteria) rho-independent terminators are scanned on both strands.
            """
        )

else:
    st.info("Paste a protein sequence and press **Optimise**.")
