"""Writes Ki.Bio's genetic-codes.txt from Biopython's copy of NCBI's gc.prt.

Usage (from the repository root): python3 -I tools/genetic_codes.py src/main/resources/io/kixi/bio/genetic-codes.txt
"""
import sys, Bio
from Bio.Data import CodonTable as C
out = sys.argv[1]
B = "TCAG"
codons = [a + b + c for a in B for b in B for c in B]
lines = [
    "# NCBI genetic codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), from gc.prt",
    f"# by way of Biopython {Bio.__version__} (Bio.Data.CodonTable), written by tools/genetic_codes.py. Do not edit by hand.",
    "# id<TAB>name<TAB>amino acids<TAB>starts<TAB>stops; 64 codons each, bases in TCAG order",
    "# (TTT TTC TTA TTG TCT ...). Amino acids: '*' where the codon only stops; starts: 'M' for a start",
    "# codon; stops: '*' for a stop codon (some codes read a stop codon as an amino acid too).",
]
for i in sorted(C.unambiguous_dna_by_id):
    t = C.unambiguous_dna_by_id[i]
    aas = "".join(t.forward_table.get(c, "*") for c in codons)
    starts = "".join("M" if c in t.start_codons else "-" for c in codons)
    stops = "".join("*" if c in t.stop_codons else "-" for c in codons)
    name = "; ".join(n for n in t.names if n and not n.startswith("SGC"))
    lines.append(f"{i}\t{name}\t{aas}\t{starts}\t{stops}")
open(out, "w").write("\n".join(lines) + "\n")
