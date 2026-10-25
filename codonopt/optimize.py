"""Simulated-annealing codon optimisation.

The optimiser freely chooses synonymous codons per residue and minimises a
penalty objective combining:
  * distance from the CAI target set by the desired expression level,
  * restriction sites (user list, both strands),
  * undesirable motifs (built-in library + user motifs, both strands),
  * GC-window extremes and homopolymer runs,
  * stable 5'-end mRNA structure,
  * bacterial rho-independent terminators (optional).
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .codon_tables import CodonTable, most_frequent_codon
from .metrics import (
    cai,
    find_motifs,
    gc_content,
    hairpin_score,
    homopolymer_runs,
    reverse_complement,
    terminator_scan,
    window_gc_extremes,
)
from .translate import clean_protein, verify_translation

# Built-in undesirable motifs (IUPAC, coding strand; scanned on both strands).
BUILTIN_MOTIFS: Dict[str, str] = {
    "polyA_signal_AATAAA": "AATAAA",
    "polyA_signal_ATTAAA": "ATTAAA",
    "TATA_box": "TATA[AT]A",
    "Shine_Dalgarno_RBS": "AGGAGG",
    "Chi_site": "GCTGGTGG",
    "splice_donor_like": "[AC]AGGT[AG]AGT",
    "splice_acceptor_like": "[TC]{10}[AC]AG",
}

# CAI targets per desired expression level (quadratic penalty around target).
EXPRESSION_CAI_TARGETS = {"high": 1.0, "medium": 0.75, "low": 0.45}


@dataclass
class OptimizationConfig:
    expression_level: str = "high"            # high | medium | low
    avoid_sites: List[str] = field(default_factory=list)   # restriction sites (IUPAC)
    avoid_motifs: List[str] = field(default_factory=list)  # extra user motifs (IUPAC)
    use_builtin_motifs: bool = True
    avoid_terminators: bool = True            # bacterial rho-independent terminators
    hairpin_window: int = 50
    gc_window: int = 20
    gc_lo: float = 0.25
    gc_hi: float = 0.75
    min_homopolymer: int = 5
    # penalty weights
    w_site: float = 100.0
    w_motif: float = 30.0
    w_gc: float = 5.0
    w_homopolymer: float = 20.0
    w_hairpin: float = 5.0
    w_terminator: float = 30.0
    w_cai: float = 50.0
    # annealing
    n_iterations: int = 6000
    t_start: float = 1.0
    t_end: float = 0.01
    seed: Optional[int] = None


@dataclass
class OptimizationResult:
    dna: str
    protein: str
    cai: float
    gc: float
    gc3: float
    hairpin: float
    sites_found: Dict[str, List[int]]
    motifs_found: Dict[str, List[int]]
    homopolymers: List[tuple]
    terminators: List[int]
    objective: float
    expression_level: str


def _motif_patterns(cfg: OptimizationConfig) -> Dict[str, str]:
    patterns: Dict[str, str] = {}
    if cfg.use_builtin_motifs:
        patterns.update(BUILTIN_MOTIFS)
    for i, m in enumerate(cfg.avoid_motifs):
        patterns[f"user_motif_{i + 1}_{m.upper()}"] = m
    return patterns


def _site_patterns(cfg: OptimizationConfig) -> Dict[str, str]:
    return {s.upper(): s for s in cfg.avoid_sites}


def _gc_window_penalty(dna: str, cfg: OptimizationConfig) -> Tuple[float, int]:
    """Soft GC-window penalty: proportional to how far windows exceed bounds."""
    dna = dna.upper()
    total, n_bad = 0.0, 0
    for i in range(0, max(1, len(dna) - cfg.gc_window + 1)):
        seg = dna[i:i + cfg.gc_window]
        gc = (seg.count("G") + seg.count("C")) / len(seg)
        excess = max(0.0, gc - cfg.gc_hi, cfg.gc_lo - gc)
        if excess > 0:
            total += excess
            n_bad += 1
    return cfg.w_gc * total, n_bad


def _objective(
    dna: str,
    protein: str,
    table: CodonTable,
    cfg: OptimizationConfig,
    log_weights: List[List[float]],
    aa_indices: List[int],
) -> Tuple[float, Dict[str, float]]:
    """Penalty objective (lower is better) plus a component breakdown."""
    target = EXPRESSION_CAI_TARGETS[cfg.expression_level]
    cai_val = cai(dna, table, floor=1e-4)
    if target >= 0.99:
        # pure maximisation: linear penalty is better conditioned near CAI ~ 1
        cai_pen = cfg.w_cai * (1.0 - cai_val)
    else:
        cai_pen = cfg.w_cai * (cai_val - target) ** 2

    sites = find_motifs(dna, _site_patterns(cfg)) if cfg.avoid_sites else {}
    n_sites = sum(len(v) for v in sites.values())
    site_pen = cfg.w_site * n_sites

    motifs = find_motifs(dna, _motif_patterns(cfg))
    n_motifs = sum(len(v) for v in motifs.values())
    motif_pen = cfg.w_motif * n_motifs

    gc_pen, n_bad_gc = _gc_window_penalty(dna, cfg)

    runs = homopolymer_runs(dna, cfg.min_homopolymer)
    homo_pen = cfg.w_homopolymer * sum(r[2] - cfg.min_homopolymer + 1 for r in runs)

    hair = hairpin_score(dna, cfg.hairpin_window)
    hair_pen = cfg.w_hairpin * hair

    term_pen = 0.0
    if cfg.avoid_terminators:
        term_pen = cfg.w_terminator * len(terminator_scan(dna))

    total = cai_pen + site_pen + motif_pen + gc_pen + homo_pen + hair_pen + term_pen
    parts = {
        "cai_penalty": cai_pen, "cai": cai_val,
        "sites": n_sites, "motifs": n_motifs,
        "bad_gc_windows": n_bad_gc, "homopolymers": len(runs),
        "hairpin": hair, "terminators": term_pen / cfg.w_terminator if cfg.avoid_terminators else 0,
    }
    return total, parts


def _initial_dna(protein: str, table: CodonTable, rng: random.Random,
                 expression_level: str) -> Tuple[str, List[List[int]], List[int]]:
    """Build a starting DNA and index structures for fast mutation.

    ``log_weights[i]`` = list of log(frequency) for codons of residue i's amino
    acid; ``aa_indices[i]`` = which codon index is currently used at residue i.
    """
    codons_by_aa: Dict[str, List[str]] = {}
    freqs_by_aa: Dict[str, List[float]] = {}
    for aa, codons in table.items():
        items = sorted(codons.items())
        codons_by_aa[aa] = [c for c, _ in items]
        freqs_by_aa[aa] = [f for _, f in items]

    log_weights, aa_indices = [], []
    out = []
    for i, aa in enumerate(protein):
        codons = codons_by_aa[aa]
        freqs = freqs_by_aa[aa]
        log_weights.append([math.log(max(f, 1e-6)) for f in freqs])
        if aa == "M" and i == 0:
            idx = codons.index("ATG") if "ATG" in codons else 0
        elif aa == "*":
            idx = freqs.index(max(freqs))
        elif expression_level == "high":
            idx = freqs.index(max(freqs))
        else:
            idx = rng.choices(range(len(codons)), weights=freqs)[0]
        aa_indices.append(idx)
        out.append(codons[idx])
    return "".join(out), log_weights, aa_indices


def _codon_lists(table: CodonTable) -> Dict[str, List[str]]:
    return {aa: sorted(codons) for aa, codons in table.items()}


def optimize(
    protein: str,
    table: CodonTable,
    config: Optional[OptimizationConfig] = None,
) -> OptimizationResult:
    """Reverse-translate a protein into an optimised coding sequence."""
    cfg = config or OptimizationConfig()
    if cfg.expression_level not in EXPRESSION_CAI_TARGETS:
        raise ValueError(f"expression_level must be one of {sorted(EXPRESSION_CAI_TARGETS)}")
    protein = clean_protein(protein)
    rng = random.Random(cfg.seed)

    codon_lists = _codon_lists(table)
    dna, log_weights, aa_indices = _initial_dna(protein, table, rng, cfg.expression_level)

    # residue -> amino acid lookup for fast mutation
    aas = list(protein)

    def dna_from_indices() -> str:
        return "".join(codon_lists[aas[i]][aa_indices[i]] for i in range(len(aas)))

    current, parts = _objective(dna, protein, table, cfg, log_weights, aa_indices)
    best_dna, best_score, best_parts = dna, current, parts

    for it in range(cfg.n_iterations):
        t = cfg.t_start * (cfg.t_end / cfg.t_start) ** (it / max(1, cfg.n_iterations - 1))
        # propose: mutate one mutable residue to a random synonymous codon
        i = rng.randrange(len(aas))
        aa = aas[i]
        choices = codon_lists[aa]
        if len(choices) < 2 or aa == "*":
            continue
        old_idx = aa_indices[i]
        new_idx = rng.randrange(len(choices))
        if new_idx == old_idx:
            continue
        aa_indices[i] = new_idx
        candidate = dna_from_indices()
        score, cand_parts = _objective(candidate, protein, table, cfg, log_weights, aa_indices)
        delta = score - current
        if delta <= 0 or rng.random() < math.exp(-delta / t):
            current, parts = score, cand_parts
            if score < best_score:
                best_dna, best_score, best_parts = candidate, score, cand_parts
        else:
            aa_indices[i] = old_idx  # revert

    final_dna = best_dna
    verify_translation(final_dna, protein)
    sites = find_motifs(final_dna, _site_patterns(cfg)) if cfg.avoid_sites else {}
    motifs = find_motifs(final_dna, _motif_patterns(cfg))

    return OptimizationResult(
        dna=final_dna,
        protein=protein,
        cai=cai(final_dna, table),
        gc=gc_content(final_dna),
        gc3=_gc3(final_dna),
        hairpin=hairpin_score(final_dna, cfg.hairpin_window),
        sites_found=sites,
        motifs_found=motifs,
        homopolymers=homopolymer_runs(final_dna, cfg.min_homopolymer),
        terminators=terminator_scan(final_dna) if cfg.avoid_terminators else [],
        objective=best_score,
        expression_level=cfg.expression_level,
    )


def _gc3(dna: str) -> float:
    third = dna[2::3]
    return sum(ch in "GC" for ch in third) / len(third) if third else 0.0
