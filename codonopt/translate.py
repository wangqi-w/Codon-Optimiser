"""Protein sequence validation and back-translation checks."""

from __future__ import annotations

from .codon_tables import STANDARD_CODE, CodonTable, most_frequent_codon

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_protein(seq: str) -> str:
    """Uppercase, strip whitespace/FASTA header, and validate a protein sequence.

    Returns the cleaned sequence. Raises ValueError with the position of the
    first invalid character (ambiguous residues such as X/B/Z are rejected
    because they have no unique codon).
    """
    if seq is None:
        raise ValueError("no protein sequence provided")
    lines = seq.strip().splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    cleaned = "".join(lines).upper().replace(" ", "").replace("\t", "")
    if not cleaned:
        raise ValueError("protein sequence is empty")
    for i, ch in enumerate(cleaned):
        if ch == "*" and i != len(cleaned) - 1:
            raise ValueError(f"internal stop character '*' at position {i + 1}")
        if ch != "*" and ch not in VALID_AA:
            raise ValueError(
                f"invalid or ambiguous residue {ch!r} at position {i + 1} "
                "(X/B/Z/J and lowercase are not supported)"
            )
    if len(cleaned) < 2:
        raise ValueError("protein sequence must be at least 2 residues long")
    return cleaned


def back_translate(dna: str) -> str:
    """Translate a DNA string (first frame) into protein; '*' for stop."""
    protein = []
    for i in range(0, len(dna) - len(dna) % 3, 3):
        protein.append(STANDARD_CODE.get(dna[i:i + 3].upper(), "X"))
    return "".join(protein)


def verify_translation(dna: str, protein: str) -> None:
    """Raise ValueError unless dna translates exactly to protein (no extra stop)."""
    translated = back_translate(dna)
    if translated != protein:
        # find first mismatch for a helpful message
        for i, (a, b) in enumerate(zip(translated, protein)):
            if a != b:
                raise ValueError(
                    f"back-translation mismatch at residue {i + 1}: "
                    f"DNA gives {a}, input was {b}"
                )
        raise ValueError(
            f"length mismatch: DNA encodes {len(translated)} residues, "
            f"input has {len(protein)}"
        )


def choose_stop_codon(table: CodonTable) -> str:
    """Most frequently used stop codon in the target organism."""
    return most_frequent_codon(table, "*")
