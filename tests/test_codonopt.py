"""Tests for the codonopt engine."""

import json

import pytest

from codonopt import (
    PRESET_LABELS,
    back_translate,
    cai,
    clean_protein,
    find_motifs,
    find_restriction_sites,
    gc_content,
    get_table,
    hairpin_score,
    homopolymer_runs,
    load_preset_tables,
    n_end_flag,
    optimize,
    parse_custom_table,
    reverse_complement,
    terminator_scan,
    validate_table,
    OptimizationConfig,
)
from codonopt.codon_tables import table_gc3
from codonopt.translate import verify_translation

ECOLI = get_table("e_coli_316407")
GFP = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGV"
    "QCFSRYPDHMKQHDFFSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNIL"
    "GHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQS"
    "ALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


# ------------------------------------------------------------------- tables

def test_all_preset_tables_valid():
    tables = load_preset_tables()
    assert len(tables) == len(PRESET_LABELS)
    for name, table in tables.items():
        assert validate_table(table) == [], f"{name} invalid"


def test_get_table_by_label():
    assert get_table("Escherichia coli") is ECOLI
    with pytest.raises(KeyError):
        get_table("no such organism")


def test_custom_table_kazusa_format():
    kazusa = """
UUU 26.1(170666) UCU 23.5(153557) UAU 18.8(122728) UGU 8.1( 52903)
UUC 18.4(120510) UCC 14.2( 92923) UAC 14.8( 96596) UGC 4.8( 31095)
UUA 26.2(170884) UCA 18.7(122028) UAA 1.1( 6913) UGA 0.7( 4447)
UUG 27.2(177573) UCG 8.6( 55951) UAG 0.5( 3312) UGG 10.4( 67789)

CUU 12.3( 80076) CCU 13.5( 88263) CAU 13.6( 89007) CGU 6.4( 41791)
CUC 5.4( 35545) CCC 6.8( 44309) CAC 7.8( 50785) CGC 2.6( 16993)
CUA 13.4( 87619) CCA 18.3(119641) CAA 27.3(178251) CGA 3.0( 19562)
CUG 10.5( 68494) CCG 5.3( 34597) CAG 12.1( 79121) CGG 1.7( 11351)

AUU 30.1(196893) ACU 20.3(132522) AAU 35.7(233124) AGU 14.2( 92466)
AUC 17.2(112176) ACC 12.7( 83207) AAC 24.8(162199) AGC 9.8( 63726)
AUA 17.8(116254) ACA 17.8(116084) AAA 41.9(273618) AGA 21.3(139081)
AUG 20.9(136805) ACG 8.0( 52045) AAG 30.8(201361) AGG 9.2( 60289)

GUU 22.1(144243) GCU 21.2(138358) GAU 37.6(245641) GGU 23.9(156109)
GUC 11.8( 76947) GCC 12.6( 82357) GAC 20.2(132048) GGC 9.8( 63903)
GUA 11.8( 76927) GCA 16.2(105910) GAA 45.6(297944) GGA 10.9( 71216)
GUG 10.8( 70337) GCG 6.2( 40358) GAG 19.2(125717) GGG 6.0( 39359)
"""
    table = parse_custom_table(kazusa)
    assert validate_table(table) == []
    assert set(table["M"]) == {"ATG"}
    assert abs(sum(table["L"].values()) - 1.0) < 1e-9
    # U-to-T conversion: TTT present, UTT absent
    assert "TTT" in table["F"] and "UTT" not in table["F"]


def test_custom_table_kazusa_incomplete_rejected():
    with pytest.raises(ValueError, match="64 codon entries"):
        parse_custom_table("UUU 26.1(170666) UUC 18.4(120510)")


def test_custom_table_json():
    raw = json.dumps(ECOLI)
    table = parse_custom_table(raw)
    assert validate_table(table) == []


def test_custom_table_csv_two_columns():
    rows = []
    for aa, codons in ECOLI.items():
        for codon, freq in codons.items():
            rows.append(f"{codon},{freq}")
    table = parse_custom_table("\n".join(rows))
    assert validate_table(table) == []


def test_custom_table_handles_uracil_and_renorm():
    text = "codon,freq\nAUG,2\nGCU,1\nGCG,1\nUAA,1\n"
    table = parse_custom_table(text)
    assert abs(sum(table["A"].values()) - 1.0) < 1e-9
    assert "ATG" in table["M"] and "AUG" not in table["M"]


def test_custom_table_rejects_bad_codon():
    with pytest.raises(ValueError):
        parse_custom_table("codon,freq\nZZZ,1\n")


# ---------------------------------------------------------------- translate

def test_clean_protein_strips_fasta_and_whitespace():
    assert clean_protein(">x\nMKV\n FLV") == "MKVFLV"


def test_clean_protein_rejects_ambiguous():
    with pytest.raises(ValueError, match="position 3"):
        clean_protein("MKXVL")


def test_clean_protein_rejects_internal_stop():
    with pytest.raises(ValueError, match="internal stop"):
        clean_protein("MK*VL")


def test_back_translate_roundtrip():
    assert back_translate("ATGAAATTCTAA") == "MKF*"


def test_verify_translation_detects_mismatch():
    with pytest.raises(ValueError):
        verify_translation("ATGAAATTCTAA", "MKL*")


# ------------------------------------------------------------------ metrics

def test_cai_perfect_and_poor():
    # all most-frequent codons -> CAI 1
    best = "".join(max(c, key=c.get) for a, c in ECOLI.items() if a != "*")
    assert cai(best, ECOLI) == pytest.approx(1.0)
    # all least-frequent codons -> low CAI
    worst = "".join(min(c, key=c.get) for a, c in ECOLI.items() if a != "*")
    assert cai(worst, ECOLI) < 0.5


def test_cai_zero_for_absent_codon():
    # AGG is near-zero in E. coli table for Arg in many CUTG builds; use a
    # synthetic table to guarantee a zero-frequency codon
    table = {"R": {"CGT": 1.0, "AGG": 0.0}}
    assert cai("CGT", table) == pytest.approx(1.0)
    assert cai("AGG", table) == 0.0


def test_gc_and_gc3():
    assert gc_content("GGCC") == 1.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ATGC") == pytest.approx(0.5)


def test_reverse_complement():
    assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI site is palindromic
    assert reverse_complement("ATGC") == "GCAT"


def test_restriction_site_both_strands():
    # non-palindromic site present only on the reverse strand must be found
    site = "GACT"  # rc = AGTC
    dna = reverse_complement(site) + "AAAA"
    hits = find_restriction_sites(dna, [site])
    assert len(hits["GACT"]) == 1
    # and a clean sequence has none
    assert find_restriction_sites("AAAAAAAAAAAAAAA", [site]) == {"GACT": []}


def test_motif_iupac_and_user():
    hits = find_motifs("AAATAAATATAAA", {"TATA_box": "TATA[AT]A"})
    assert hits["TATA_box"]


def test_homopolymer_runs():
    assert homopolymer_runs("AAAAAAGGGTT", 5) == [(0, "A", 6)]
    assert homopolymer_runs("AATATA", 5) == []


def test_hairpin_score_detects_stem_loop():
    stem_loop = "ACGT" + "TTTTTT" + "ACGT"  # perfect 4-bp stem, 6-nt loop
    assert hairpin_score(stem_loop * 2, window=50) > 0
    assert hairpin_score("AAAAAAAAAAAA", window=50) == 0.0


def test_terminator_scan_finds_classic_terminator():
    # GC-rich inverted repeat followed by poly-U (T's in DNA)
    seq = "AAAAAAGCGCGCGCGCGCTTTTTTTTTTAAAAA"
    assert terminator_scan(seq)


def test_n_end_flag():
    assert n_end_flag("MR") is not None  # mature N-term Arg
    assert n_end_flag("MAGST") is None
    assert n_end_flag("MD") is not None  # secondary destabilising


# ---------------------------------------------------------------- optimizer

def test_optimize_high_expression_ecoli():
    cfg = OptimizationConfig(
        expression_level="high", avoid_sites=["GAATTC", "GCGC"],
        seed=42,
    )
    res = optimize(GFP, ECOLI, cfg)
    verify_translation(res.dna, res.protein)
    assert res.cai > 0.9
    assert res.sites_found == {"GAATTC": [], "GCGC": []}
    assert res.gc > 0.3 and res.gc < 0.75


def test_optimize_expression_levels_order():
    scores = {}
    for level in ("high", "medium", "low"):
        cfg = OptimizationConfig(expression_level=level, seed=7, n_iterations=1500)
        res = optimize(GFP, ECOLI, cfg)
        verify_translation(res.dna, res.protein)
        scores[level] = res.cai
    assert scores["high"] > scores["medium"] > scores["low"]


def test_optimize_motif_avoidance():
    cfg = OptimizationConfig(
        expression_level="high", use_builtin_motifs=True,
        avoid_motifs=["GCTGGTGG"], seed=3, n_iterations=1500,
    )
    res = optimize(GFP, ECOLI, cfg)
    assert all(len(v) == 0 for v in res.motifs_found.values())


def test_optimize_custom_table():
    table = parse_custom_table(json.dumps(get_table("h_sapiens_9606")))
    cfg = OptimizationConfig(expression_level="high", seed=1, n_iterations=1000)
    res = optimize("MKVL", table, cfg)
    verify_translation(res.dna, "MKVL")


def test_optimize_short_peptide_and_seed_reproducibility():
    cfg = OptimizationConfig(expression_level="medium", seed=11, n_iterations=500)
    r1 = optimize("MKV", ECOLI, cfg)
    r2 = optimize("MKV", ECOLI, cfg)
    assert r1.dna == r2.dna


def test_optimize_rejects_bad_expression_level():
    with pytest.raises(ValueError):
        optimize("MKV", ECOLI, OptimizationConfig(expression_level="max"))


def test_optimize_start_codon_fixed():
    cfg = OptimizationConfig(expression_level="low", seed=5, n_iterations=500)
    res = optimize("MKVL", ECOLI, cfg)
    assert res.dna.startswith("ATG")


def test_table_gc3_sane():
    g = table_gc3(ECOLI)
    assert 0.2 < g < 1.0
