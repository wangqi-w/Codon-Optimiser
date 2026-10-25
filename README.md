# Codon Optimiser

A browser-based web app that reverse-translates any **protein sequence** into a
**DNA sequence** optimised for a target organism and desired expression level.

## Features

- **Target organism**: 9 built-in codon usage tables (Kazusa CUTG-derived):
  *E. coli*, *B. subtilis*, *S. cerevisiae*, *C. elegans*, *D. melanogaster*,
  human, mouse (x2), chicken — plus a **custom codon table** pasted in the
  codon usage database (Kazusa CUTG) web format, e.g.
  `UUU 26.1(170666) UCU 23.5(153557) ...` (all 64 codons; U or T spelling).
- **Desired expression level**: high (maximise CAI), medium (moderate CAI), or
  low (deliberate rare-codon enrichment, e.g. for folding-sensitive proteins).
- **Restriction sites to avoid**: pick from common sites or specify any
  IUPAC-coded site; scanned on both strands.
- **Undesirable motifs**: built-in screen for premature polyA signals
  (AATAAA/ATTAAA), TATA boxes, Shine-Dalgarno/RBS motifs, Chi sites, cryptic
  splice donor/acceptor motifs, and rho-independent terminators (bacteria);
  custom motifs also supported.
- **Additional factors**: 5'-end mRNA hairpin minimisation (first 50 nt),
  sliding-window GC balance, homopolymer-run avoidance, and an N-end-rule
  note for the mature N-terminus.
- **Outputs**: optimised DNA (FASTA + text download), a CSV design report
  (sequence name, parameters, optimised sequence, CAI / GC / GC3 / hairpin
  metrics), a quality report, and a guaranteed back-translation check (DNA
  translates exactly to your input protein).

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app expects the `codonopt` package folder next to `app.py` (as in this
repository).

## Deploy to Streamlit Community Cloud (free)

1. Push this repository to GitHub.
2. Sign up at <https://share.streamlit.io> with a free account.
3. Click **New app**, select your repo, set the main file to `app.py`, and deploy.

## Project layout

```
codonopt/           optimisation engine (pure Python)
  codon_tables.py   preset + custom codon usage tables
  translate.py      protein validation and back-translation checks
  metrics.py        CAI, GC, motif/restriction scans, hairpin heuristic, N-end rule
  optimize.py       simulated-annealing optimiser
  data/             bundled Kazusa CUTG codon tables (JSON)
tests/              pytest suite
app/app.py          Streamlit web app
```

## Design rationale

See `report_codon_optimisation_research.md` for the literature summary of the
sequence properties that drive expression (CAI, 5' mRNA structure, GC
composition, regulatory motifs, restriction sites, N-end rule) and how the
optimiser weights them.

## Limitations

- The 5' hairpin score is a lightweight stem-loop heuristic, not a
  thermodynamic folding model; verify final candidates with an RNA-folding
  tool if 5' structure is critical.
- Expression-level tuning is CAI-targeted; it does not model codon-pair bias
  or codon harmonisation along the message.
- Ambiguous protein residues (X, B, Z, J) are rejected.
