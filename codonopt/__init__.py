"""codonopt: reverse-translate proteins into expression-tuned DNA sequences.

Modules
-------
codon_tables : preset + custom codon usage tables
translate    : protein validation and back-translation checks
metrics      : CAI, GC, motif/restriction scans, 5' hairpin heuristic, N-end rule
optimize     : simulated-annealing codon optimisation
"""

from .codon_tables import (
    PRESET_LABELS,
    get_table,
    load_preset_tables,
    parse_custom_table,
    validate_table,
)
from .metrics import (
    cai,
    find_motifs,
    find_restriction_sites,
    gc3_content,
    gc_content,
    hairpin_score,
    homopolymer_runs,
    n_end_flag,
    reverse_complement,
    terminator_scan,
)
from .optimize import OptimizationConfig, optimize
from .translate import back_translate, clean_protein

__version__ = "1.0.0"

__all__ = [
    "PRESET_LABELS",
    "OptimizationConfig",
    "back_translate",
    "cai",
    "clean_protein",
    "find_motifs",
    "find_restriction_sites",
    "gc3_content",
    "gc_content",
    "get_table",
    "hairpin_score",
    "homopolymer_runs",
    "load_preset_tables",
    "n_end_flag",
    "optimize",
    "parse_custom_table",
    "reverse_complement",
    "terminator_scan",
    "validate_table",
]
