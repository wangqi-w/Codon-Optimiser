"""Codon usage tables: presets (Kazusa CUTG via python_codon_tables) and custom input.

A codon table is a dict mapping one-letter amino acid codes (including "*" for
stop) to {codon: relative_frequency} with frequencies summing to ~1 per amino
acid. Codons use T (DNA alphabet).
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
from typing import Dict, Union

# Standard genetic code: codon -> amino acid (DNA alphabet).
STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "AGT": "S", "AGC": "S",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY*"

CodonTable = Dict[str, Dict[str, float]]

# Friendly display names for the bundled preset tables.
PRESET_LABELS = {
    "e_coli_316407": "Escherichia coli",
    "b_subtilis_1423": "Bacillus subtilis",
    "s_cerevisiae_4932": "Saccharomyces cerevisiae",
    "c_elegans_6239": "Caenorhabditis elegans",
    "d_melanogaster_7227": "Drosophila melanogaster",
    "h_sapiens_9606": "Homo sapiens",
    "m_musculus_10090": "Mus musculus",
    "m_musculus_domesticus_10092": "Mus musculus domesticus",
    "g_gallus_9031": "Gallus gallus",
}

_DATA_FILE = os.path.join(os.path.dirname(__file__), "data", "codon_tables.json")
_PRESET_CACHE = None


def load_preset_tables() -> Dict[str, CodonTable]:
    """Load all bundled preset tables (Kazusa CUTG-derived), cached."""
    global _PRESET_CACHE
    if _PRESET_CACHE is None:
        with open(_DATA_FILE) as fh:
            _PRESET_CACHE = json.load(fh)
    return _PRESET_CACHE


def get_table(name: str) -> CodonTable:
    """Look up a preset table by key or friendly label."""
    tables = load_preset_tables()
    if name in tables:
        return tables[name]
    for key, label in PRESET_LABELS.items():
        if label.lower() == name.strip().lower():
            return tables[key]
    raise KeyError(f"Unknown codon table: {name!r}. Available: {sorted(tables)}")


def validate_table(table: CodonTable, tol: float = 0.05) -> list[str]:
    """Validate a codon table; return a list of problems (empty = valid)."""
    problems: list[str] = []
    for aa in AMINO_ACIDS:
        if aa not in table or not table[aa]:
            problems.append(f"missing amino acid {aa!r}")
            continue
        codons = table[aa]
        expected = [c for c, a in STANDARD_CODE.items() if a == aa]
        for codon in codons:
            codon = codon.upper().replace("U", "T")
            if codon not in expected:
                problems.append(f"{codon} is not a codon for {aa}")
        total = sum(codons.values())
        if abs(total - 1.0) > tol:
            problems.append(f"frequencies for {aa} sum to {total:.4f} (expected ~1)")
        if any(v < 0 for v in codons.values()):
            problems.append(f"negative frequency for {aa}")
    extra = set(table) - set(AMINO_ACIDS)
    if extra:
        problems.append(f"unexpected keys: {sorted(extra)}")
    return problems


def _normalise(table: CodonTable) -> CodonTable:
    """Uppercase codons, convert U to T, and renormalise each amino acid to 1."""
    out: CodonTable = {}
    for aa, codons in table.items():
        aa = aa.strip()
        if aa == "U":  # some tables use U for stop
            aa = "*"
        norm: Dict[str, float] = {}
        for codon, freq in codons.items():
            codon = codon.strip().upper().replace("U", "T")
            if len(codon) != 3 or set(codon) - set("ACGT"):
                raise ValueError(f"invalid codon {codon!r}")
            norm[codon] = float(freq)
        total = sum(norm.values())
        if total <= 0:
            raise ValueError(f"frequencies for {aa} sum to {total}")
        out[aa] = {c: f / total for c, f in norm.items()}
    return out


# Kazusa CUTG web-format token: e.g. "UUU 26.1(170666)" (freq per 1000 + count).
# Both RNA (U) and DNA (T) spellings of the codons are accepted.
_KAZUSA_TOKEN = re.compile(r"([ACGUT]{3})\s+([\d.]+)\s*(?:\((\d+)\))?")


def parse_kazusa_table(text: str) -> CodonTable:
    """Parse a codon usage table in the Kazusa CUTG web-display format.

    Expected layout (blank lines and column grouping are ignored)::

        UUU 26.1(170666) UCU 23.5(153557) UAU 18.8(122728) UGU 8.1( 52903)
        UUC 18.4(120510) UCC 14.2( 92923) UAC 14.8( 96596) UGC 4.8( 31095)
        ...

    Codons are converted from U to T and frequencies are renormalised to sum
    to 1 per amino acid. Raises ValueError if no valid tokens are found.
    """
    table: CodonTable = {}
    n_tokens = 0
    for codon, freq, _count in _KAZUSA_TOKEN.findall(text):
        n_tokens += 1
        codon = codon.replace("U", "T")
        aa = STANDARD_CODE[codon]
        table.setdefault(aa, {})[codon] = float(freq)
    if n_tokens < 64:
        raise ValueError(
            f"expected 64 codon entries in Kazusa CUTG format "
            f"(e.g. 'UUU 26.1(170666)'), found {n_tokens}"
        )
    return _normalise(table)


def parse_custom_table(text: str) -> CodonTable:
    """Parse a custom codon table.

    Accepted formats (tried in order):
      1. Kazusa CUTG web format: ``UUU 26.1(170666) UCU 23.5(153557) ...``
      2. JSON: {"A": {"GCT": 0.4, ...}, ...}
      3. CSV/TSV: rows of ``codon, frequency`` or ``amino_acid, codon, frequency``.
    """
    text = text.strip()
    if not text:
        raise ValueError("empty codon table")
    if _KAZUSA_TOKEN.search(text):
        return parse_kazusa_table(text)
    if text[0] in "{[":
        raw = json.loads(text)
        return _normalise(raw)

    # sniff delimiter
    sample = text[:2048]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t ")
    except csv.Error:
        dialect = csv.excel
    rows = [r for r in csv.reader(io.StringIO(text), dialect) if r and r[0].strip()]
    if not rows:
        raise ValueError("no rows found in codon table")

    # drop a header row if present
    first = [c.strip().lower() for c in rows[0]]
    if any(h in ("codon", "frequency", "freq", "aa", "amino_acid") for h in first):
        rows = rows[1:]
    if not rows:
        raise ValueError("codon table has no data rows")

    table: CodonTable = {}
    if len(rows[0]) == 2:
        for codon, freq in rows:
            codon = codon.strip().upper().replace("U", "T")
            if codon not in STANDARD_CODE:
                raise ValueError(f"unknown codon {codon!r}")
            aa = STANDARD_CODE[codon]
            table.setdefault(aa, {})[codon] = float(freq)
    elif len(rows[0]) == 3:
        for aa, codon, freq in rows:
            aa = aa.strip().upper()
            if aa == "U":
                aa = "*"
            if aa not in AMINO_ACIDS:
                raise ValueError(f"unknown amino acid {aa!r}")
            codon = codon.strip().upper().replace("U", "T")
            table.setdefault(aa, {})[codon] = float(freq)
    else:
        raise ValueError("rows must have 2 columns (codon, freq) or 3 (aa, codon, freq)")
    return _normalise(table)


def most_frequent_codon(table: CodonTable, aa: str) -> str:
    """Return the most frequent codon for an amino acid."""
    codons = table[aa]
    return max(codons, key=codons.get)


def table_gc3(table: CodonTable) -> float:
    """Usage-weighted mean GC content at third codon positions."""
    num = den = 0.0
    for aa, codons in table.items():
        if aa == "*":
            continue
        for codon, freq in codons.items():
            num += freq * (codon[2] in "GC")
            den += freq
    return num / den if den else 0.0
