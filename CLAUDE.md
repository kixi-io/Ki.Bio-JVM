# Ki.Bio-JVM

The `bio` type bag (Ki type bags design, sections 2, 5.5 and 11; the design is in Ki.Docs and the RevDB project docs). Kotlin 2.3.0 with `explicitApi()`, Kotest 6.0.7, JDK 21.

## Composite build
Includes `../Ki.Core-JVM` and `../Ki.Sci-JVM` (the main dependencies) and `../Ki.KD-JVM`, `../Ki.KS-JVM`, `../RevDB-JVM` (tests only). Ask before adding a dependency.

## Rules
- The bag is `BioBag` (a class with a public no-argument constructor, listed in `META-INF/services/io.kixi.types.TypeBag`); nothing registers it. Its manifest, `META-INF/ki/bags/io.kixi.bio.BioBag.properties`, is what discovery reads instead of creating the bag (type bags design, 13): it says `preload=true`, and lists the bag's types, literals, functions, constants and members, so regenerate it with `TypeBags.manifestOf(bag, preload = true)` whenever a name changes (`BagCorpus` fails until it matches). Tests reach it with `TypeBags.byName("bio")`.
- Stored encodings and sort keys (`Packing.kt`) never change for a released version: bump `BioBag.VERSION` and keep decoding older versions.
- `Sequence` (Kotlin `BioSequence`) has no literal and no samples; its members are the union of the three types' members, dispatched to the value's own type, so a `type=Sequence` attribute can be queried with any of them. A member added to `DNA`, `RNA` or `Protein` appears there by itself.
- Tm follows IDT and pI follows ExPASy (design decisions 15 and 16). `BioReferenceTest` compares them, translation and masses with Biopython: regenerate `src/test/resources/io/kixi/bio/biopython-reference.tsv` with `tools/reference.py`, never edit it. Biopython clamps its pI search to 4.05 to 12; the test allows for that.
- The genetic codes come from NCBI by way of Biopython (`tools/genetic_codes.py` writes `genetic-codes.txt`); never by hand.
- Masses come from Ki.Sci formulas (`Masses.kt` holds residue formulas, not weights); never add a mass table.
- Letters are kept as written (decision 13): no rotation or strand canonicalization in a value, its text or its sort key.
- Function and constant names are global across bags (and KQ's built-ins win in KQ): check the other bags before adding one. Generic names (`length`, `count`, `find`, `reverse`, ...) are members, never bag functions (decision 20).
- `BagCorpus.kt` (the corpus round trip) is shared with Ki.Math-JVM's and Ki.Sci-JVM's tests; keep the copies alike.

## Working agreements
- Commit straight to `master`; never create a branch unless Dan asks. Code and docs in separate commits, author `dleuck <dan@kixi.io>`. Don't publish or change the version without asking.
