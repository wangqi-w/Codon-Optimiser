"""Sequence metrics: CAI, GC, motif/restriction scans, 5' hairpin heuristic, N-end rule."""

from __future__ import annotations

import math
import re
from functools import lru_cache
from typing import Dict, List, Optional

from .codon_tables import CodonTable

# ---------------------------------------------------------------- IUPAC / DNA

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@lru_cache(maxsize=4096)
def _rc_cached(dna: str) -> str:
    return dna.translate(COMPLEMENT)[::-1]


def reverse_complement(dna: str) -> str:
    return _rc_cached(dna)


@lru_cache(maxsize=512)
def iupac_to_regex(pattern: str) -> re.Pattern:
    """Compile an IUPAC ambiguity-code motif into a regex (case-insensitive).

    Regex character classes like ``[AT]`` and quantifiers like ``{10}`` are
    passed through unchanged.
    """
    pattern = pattern.strip().upper().replace("U", "T")
    out = []
    for ch in pattern:
        if ch in IUPAC:
            out.append(IUPAC[ch])
        elif ch.isalpha():
            raise ValueError(f"invalid IUPAC character {ch!r} in motif {pattern!r}")
        else:  # regex syntax: [ ] { } etc.
            out.append(ch)
    return re.compile("".join(out))


# ------------------------------------------------------------------ CAI / GC

def _codon_weights(table: CodonTable, floor: Optional[float] = None) -> Dict[str, float]:
    """Relative adaptiveness (codon freq / max freq for its amino acid)."""
    weights: Dict[str, float] = {}
    for aa, codons in table.items():
        max_f = max(codons.values())
        for codon, freq in codons.items():
            w = freq / max_f if max_f > 0 else 0.0
            if floor is not None:
                w = max(w, floor)
            weights[codon] = w
    return weights


def cai(dna: str, table: CodonTable, floor: Optional[float] = None) -> float:
    """Codon Adaptation Index (Sharp & Li 1987).

    With ``floor=None`` (default) codons absent from the table / zero-frequency
    codons give a true CAI of 0. Pass a small floor (e.g. 1e-3) for a smoothed
    score usable inside an optimiser.
    """
    dna = dna.upper().replace("U", "T")
    weights = _codon_weights(table, floor=floor)
    log_sum, n = 0.0, 0
    for i in range(0, len(dna) - len(dna) % 3, 3):
        codon = dna[i:i + 3]
        if codon == "ATG":  # methionine: single codon, weight 1
            n += 1
            log_sum += 0.0
            continue
        w = weights.get(codon, 0.0)
        if w <= 0:
            return 0.0
        log_sum += math.log(w)
        n += 1
    if n == 0:
        return 0.0
    return math.exp(log_sum / n)


def gc_content(dna: str) -> float:
    dna = dna.upper().replace("U", "T")
    if not dna:
        return 0.0
    return sum(ch in "GC" for ch in dna) / len(dna)


def gc3_content(dna: str) -> float:
    """GC content at third codon positions only."""
    dna = dna.upper().replace("U", "T")
    third = dna[2::3]
    if not third:
        return 0.0
    return sum(ch in "GC" for ch in third) / len(third)


def window_gc_extremes(dna: str, window: int = 20, lo: float = 0.25, hi: float = 0.75) -> List[tuple]:
    """Sliding windows whose GC fraction falls outside [lo, hi]."""
    dna = dna.upper()
    bad = []
    for i in range(0, max(1, len(dna) - window + 1)):
        seg = dna[i:i + window]
        gc = sum(ch in "GC" for ch in seg) / len(seg)
        if gc < lo or gc > hi:
            bad.append((i, round(gc, 3)))
    return bad


# ------------------------------------------------------------- motif / sites

def find_motifs(dna: str, patterns: Dict[str, str]) -> Dict[str, List[int]]:
    """Scan DNA (both strands) for IUPAC motifs.

    ``patterns`` maps a motif name to an IUPAC pattern string. Returns
    {name: [match positions on the given strand]}.
    """
    dna = dna.upper().replace("U", "T")
    rc = reverse_complement(dna)
    hits: Dict[str, List[int]] = {}
    for name, pattern in patterns.items():
        rx = iupac_to_regex(pattern)
        fwd = [m.start() for m in rx.finditer(dna)]
        rev = [m.start() for m in rx.finditer(rc)]
        hits[name] = fwd + rev
    return hits


def find_restriction_sites(dna: str, sites: List[str]) -> Dict[str, List[int]]:
    """Locate restriction sites (IUPAC allowed) on both strands."""
    patterns = {s.upper(): s for s in sites}
    return find_motifs(dna, patterns)


def homopolymer_runs(dna: str, min_run: int = 5) -> List[tuple]:
    """Runs of a single nucleotide of length >= min_run: [(start, base, length)]."""
    dna = dna.upper()
    runs = []
    i = 0
    while i < len(dna):
        j = i
        while j < len(dna) and dna[j] == dna[i]:
            j += 1
        if j - i >= min_run:
            runs.append((i, dna[i], j - i))
        i = j
    return runs


def terminator_scan(dna: str, min_stem: int = 4, min_polyu: int = 4) -> List[int]:
    """Heuristic rho-independent terminator positions (bacteria).

    Finds GC-rich inverted repeats (stem >= min_stem, loop 3-8 nt) followed
    within 15 nt by a run of >= min_polyu T's on the coding strand. Scans
    upstream of poly-T runs only, since the poly-U tail is mandatory.
    """
    dna = dna.upper()
    hits = []
    n = len(dna)
    for m in re.finditer(f"T{{{min_polyu},}}", dna):
        tail_start = m.start()
        region_start = max(0, tail_start - 15 - 8)
        for i in range(region_start, tail_start - 2 * min_stem - 3 + 1):
            if i < 0:
                continue
            for stem in (8, 7, 6, 5, 4):
                if stem < min_stem:
                    break
                stem_seq = dna[i:i + stem]
                rc_stem = reverse_complement(stem_seq)
                for loop in range(3, 9):
                    j = i + stem + loop
                    k = j + stem
                    if k > tail_start + 15 or k > n:
                        break
                    if dna[j:k] == rc_stem:
                        gc_pairs = sum(1 for a in stem_seq if a in "GC")
                        if gc_pairs >= 3:
                            hits.append(i)
                    # keep scanning loops; longest stem at this i wins via max below
            if hits and hits[-1] == i:
                break  # already recorded a terminator starting at i
    return sorted(set(hits))


# ------------------------------------------------------- 5' hairpin heuristic

def hairpin_score(dna: str, window: int = 50) -> float:
    """Heuristic 5'-end mRNA secondary-structure score in [0, 1].

    Detects perfect inverted-repeat stem-loops (stem >= 4 contiguous
    Watson-Crick pairs, loop 3-8 nt) within the first ``window`` nt and scores
    by the longest stem found: (max_stem - 3) / 5, capped at 1. Higher = more
    stable structure near the start codon = worse for translation initiation.
    Lightweight approximation (no thermodynamic model).
    """
    seg = dna.upper().replace("U", "T")[:window]
    n = len(seg)
    if n < 11:  # 4 + 3 + 4 minimum hairpin
        return 0.0
    max_stem = 0
    for i in range(n - 2 * 4 - 3):
        for stem in range(8, 3, -1):
            found = False
            for loop in range(3, 9):
                j = i + stem + loop
                k = j + stem
                if k > n:
                    break
                if seg[i:i + stem] == reverse_complement(seg[j:k]):
                    # require at least one GC pair so AT-rich stems don't dominate
                    if any(a in "GC" for a in seg[i:i + stem]):
                        max_stem = max(max_stem, stem)
                        found = True
                        break
            if found:
                break
    return min(1.0, max(0, max_stem - 3) / 5)


# ------------------------------------------------------------------ N-end rule

# Destabilising N-terminal residues after Met removal (Varshavsky N-end rule).
N_END_PRIMARY = {"R", "K", "H", "L", "F", "Y", "W", "I"}  # primary destabilising
N_END_SECONDARY = {"D", "E", "C"}  # secondary (via N-terminal modification)
# Residues permitting initiator-Met removal (small side chains).
MET_REMOVAL_OK = set("GASCTPV")


def n_end_flag(protein: str) -> Optional[str]:
    """Return a human-readable N-end-rule note, or None if the N-terminus is fine."""
    if not protein:
        return None
    # residue 1 of the mature protein: Met is removed if residue 2 is small
    body = protein[1:] if protein[0] == "M" else protein
    if not body:
        return None
    first = body[0]
    if first in N_END_PRIMARY:
        return (
            f"N-terminal {first} is a primary destabilising residue under the "
            "N-end rule; consider adding a stabilising N-terminal residue "
            "(e.g. A, G, S, M) if protein half-life matters."
        )
    if first in N_END_SECONDARY:
        return (
            f"N-terminal {first} is a secondary destabilising residue under the "
            "N-end rule (can be converted to a primary destabilising residue)."
        )
    return None
