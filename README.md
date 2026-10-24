# Ki.Bio-JVM

The `bio` type bag: DNA, RNA and protein sequences for KD, KS, KQ and RevDB. A
type bag is a library that adds value types, functions and constants to every Ki
layer at once (the Ki type bags design, in Ki.Docs). Ki.Bio depends on Ki.Core
and Ki.Sci (for formulas and masses).

```kd
oligo name="M13F" seq=.dna(GTAAAACGACGGCCAGT)
plasmid name="pToy" seq=.dna(ds circular "TTCATGGCCTAAGAA")
peptide name="ub-N" seq=.protein("MQIFVKTLTGK")
```

```ks
use ki.bio.*
let primer = .dna("GTAAAACGACGGCCAGT")
say primer.reverseComplement()        // ACTGGCCGTCGTTTTAC
say "GAYGG" in primer                 // true: IUPAC patterns, Y is C or T
say primer.tm()                       // about 53.58°C (IDT's defaults: 50 mM Na+, 250 nM)
say primer.tm(50mM, 250nM, 2mM, 0.2mM) // about 60.39°C in a PCR buffer (Na+, oligo, Mg2+, dNTPs)
say .dna(ds circular "TTCAAAAAAAGAA").find("GAATTC")  // [11]: the site wraps the origin
say translate("ATGGCCTAA")            // MA*
say .protein("MQIFVKTLTGK").mass      // 1265.581Da
```

```kq
from o in oligos where o.seq == .dna(gtaaaacgacggccagt) select o.name
from o in oligos from p in plasmids where o.seq in p.seq or o.seq.reverseComplement() in p.seq select o.name, p.name
```

## Types

| Type | Literal | Orders by |
|---|---|---|
| `Sequence` | none: any of the three below, for a schema that takes them all | |
| `DNA` | `.dna("…")`, `.dna(ds "…")`, `.dna(ds circular "…")` | letters, then strandedness and topology |
| `RNA` | `.rna("…")` (also `ds`, `circular`) | letters, then strandedness and topology |
| `Protein` | `.protein("…")` | letters |

- **Letters as written**, uppercased, 5' to 3' (N- to C-terminus): no canonical
  rotation or strand, so `==` compares the writing and `sameMolecule` asks the
  chemical question (the same circle from another origin, or the other strand of
  a duplex). A literal is letters only (either case, no spaces); the functions
  read looser text: `dna(text)` drops whitespace, and `seq(text)` also drops FASTA headers and
  numbering and works out the kind.
- **Alphabets**: DNA `ACGT`, RNA `ACGU`, both with the IUPAC codes
  `RYSWKMBDHVN`; protein the 20 amino acids, `U` (selenocysteine), `O`
  (pyrrolysine), `B`, `Z`, `J`, `X` and `*` (a stop).
- **Positions** are 1-based and inclusive, as in GenBank. `find`, `count` and
  `in` take IUPAC patterns (`GGATYC`), count overlaps, and run across the origin
  of a circular sequence; so does `slice`.
- **Stored** packed: 2 bits a base for plain sequences, 4 with IUPAC codes, 5 a
  residue for proteins (a 10 kb plasmid takes 2.5 KB). Indexable, so a primer is
  found by its letters; for long sequences index `digest` (GA4GH refget
  `sha512t24u`, `SQ.…`).

Members of every sequence: `text`, `kind`, `length`, `formula`, `mass`,
`monoisotopicMass`, `digest`, `find`, `count`, `contains`, `slice`, `reverse`,
`sameMolecule`. Of DNA and RNA: `ds`, `circular`, `gc` (percent of the bases known
to be G/C or A/T, as Biopython's `gc_fraction`; nil when none is), `complement`,
`reverseComplement`, `codons`, `translate(frame, code)`, `orfs(minLength, code)`,
`orf`; DNA adds `transcribe` and `tm`, RNA `reverseTranscribe`. Of a protein:
`pI`, `charge(pH)`, `extinction280`, `extinction280Reduced`.

Functions: `dna(x)`, `rna(x)`, `protein(x)`, `seq(x, kind?)`, `revcomp`,
`transcribe`, `translate(x, frame?, code?)`, `gc`, `tm(x, Na+?, oligo?, Mg2+?,
dNTPs?)`, `pI`, `a260(absorbance, kind)` (µg/mL: 50 for dsDNA, 33 ssDNA, 40
RNA), `hamming`, `levenshtein`. Generic names (`length`, `count`, `find`,
`reverse`, ...) are members only, so `use ki.bio.*` binds nothing a script is
likely to define itself.

## Methods

- **Melting temperature** as IDT's OligoAnalyzer computes it: nearest-neighbor
  thermodynamics (Allawi and SantaLucia 1997; SantaLucia 1998) with the salt
  corrections of Owczarzy et al. 2004 (Na+) and 2008 (Mg2+, with dNTPs binding
  it), at IDT's defaults of 50 mM Na+, 250 nM oligo and no Mg2+ or dNTPs. The
  oligo is taken to be in excess over its target, as IDT assumes (k is the oligo
  concentration). Plain A, C, G and T, linear, at every length.
- **Isoelectric point** as ExPASy's Compute pI/Mw finds it (Bjellqvist pKs,
  terminal residues included), by bisection to 0.0001.
- **Extinction coefficient at 280 nm** (Pace et al. 1995): 5500 per Trp, 1490
  per Tyr, 125 per cystine.
- **Masses** from Ki.Sci formulas with the CIAAW 2021 standard atomic weights:
  DNA and RNA with a 5' hydroxyl, as synthesized; a circular sequence has no
  ends; a protein adds water; ambiguous letters have no mass.
- **Genetic codes**: the 27 NCBI tables, by number (`translate(x, 1, 2)` reads
  frame 1 with the vertebrate mitochondrial code); an ambiguous codon
  translates when every reading of it agrees (`TTY` is F), as `B`, `Z` or `J`
  when the readings are D or N, E or Q, I or L, and is `X` otherwise.

Not yet: alignment, oligo modifications (5' phosphate, dyes), RNA duplex Tm,
ε260, and lengths in bp, nt and aa units.

## Build

A composite build with the sibling checkouts `../Ki.Core-JVM`, `../Ki.Sci-JVM`,
and, for the tests only, `../Ki.KD-JVM`, `../Ki.KS-JVM` and `../RevDB-JVM`.
Kotlin 2.3.0, Kotest 6.0.7, JDK 21.

```
./gradlew test
```

Tm, pI, translation and masses are checked against Biopython 1.88 (the same
methods: `MeltingTemp.Tm_NN` with `DNA_NN3` and salt correction 7,
`IsoelectricPoint`), Tm within 0.001 °C over 1,179 cases. `tools/reference.py`
writes the reference values and `tools/genetic_codes.py` the genetic code table;
both regenerate identically. The corpus round trip puts every sample value
through KD, KS and RevDB.

MIT license. The genetic codes are NCBI's, by way of Biopython.
