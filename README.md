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
