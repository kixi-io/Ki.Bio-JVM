package io.kixi.bio

import io.kixi.sci.Element
import io.kixi.sci.Elements
import io.kixi.sci.Formula

/**
 * Sequence formulas, so masses come from Ki.Sci's standard atomic weights
 * (CIAAW 2021, design decision 10) and monoisotopic masses from its isotopes.
 *
 * A residue is its monomer less one water, as it sits in the chain. A linear
 * strand of DNA or RNA has a 5' hydroxyl, as an oligo is synthesized (IDT's
 * convention): its residues less HPO3 plus H2O, one phosphate fewer than
 * residues. A circular strand is its residues alone. A protein is its residues
 * plus one water. An ambiguous letter (N, R, X, ...) has no formula.
 */
internal object Masses {

    // Element order of the count arrays.
    private val SYMBOLS = arrayOf("C", "H", "N", "O", "P", "S", "Se")
    private val ELEMENTS: Array<Element> by lazy { Array(SYMBOLS.size) { Elements.bySymbol(SYMBOLS[it])!! } }

    private fun f(c: Int, h: Int, n: Int, o: Int, p: Int = 0, s: Int = 0, se: Int = 0) = intArrayOf(c, h, n, o, p, s, se)

    /** Deoxynucleotide residues: dNMP less H2O (dA C10H12N5O5P, 313.21 Da). */
    private val DNA = mapOf(
        'A' to f(10, 12, 5, 5, 1), 'C' to f(9, 12, 3, 6, 1), 'G' to f(10, 12, 5, 6, 1), 'T' to f(10, 13, 2, 7, 1),
    )

    /** Ribonucleotide residues: NMP less H2O (A C10H12N5O6P, 329.21 Da). */
    private val RNA = mapOf(
        'A' to f(10, 12, 5, 6, 1), 'C' to f(9, 12, 3, 7, 1), 'G' to f(10, 12, 5, 7, 1), 'U' to f(9, 11, 2, 8, 1),
    )

    /**
     * Threose nucleotide residues, **derived and unverified** ([warnUnverified]).
     *
     * TNA's sugar is α-L-threofuranose (C4H8O4) where DNA has 2-deoxyribose
     * (C5H10O4), so every residue is the deoxynucleotide residue less CH2:
     * tA C9H10N5O5P, 299.18 Da. The backbone is a 3'-2' phosphodiester rather
     * than 5'-3', but it is still one phosphate per internucleotide bond, so the
     * end correction in [strand] is DNA's.
     *
     * No published residue table was used. Check against one before relying on a
     * TNA mass.
     */
    private val TNA = mapOf(
        'A' to f(9, 10, 5, 5, 1), 'C' to f(8, 10, 3, 6, 1), 'G' to f(9, 10, 5, 6, 1), 'T' to f(9, 11, 2, 7, 1),
    )

    /**
     * Locked nucleotide residues, **derived and unverified** ([warnUnverified]).
     *
     * An LNA sugar is ribose closed by a 2'-O,4'-C methylene bridge: the bridge
     * adds CH2 and the two new bonds take an H from the 2' oxygen and one from
     * C4', so the sugar is ribose plus one carbon and every residue is the
     * ribonucleotide residue plus C. LNA is made with thymine rather than uracil,
     * so T is the 5-methyl residue (uracil plus CH2) with the bridge:
     * lA C11H12N5O6P, 341.22 Da; lT C11H13N2O8P, 332.20 Da.
     *
     * No published residue table was used.
     */
    private val LNA = mapOf(
        'A' to f(11, 12, 5, 6, 1), 'C' to f(10, 12, 3, 7, 1), 'G' to f(11, 12, 5, 7, 1), 'T' to f(11, 13, 2, 8, 1),
    )

    /**
     * Peptide nucleic acid residues, derived but **reproducing the values PNA mass
     * calculators use** (A 275.27, C 251.25, G 291.27, T 266.26 Da), so these are
     * the firmest of the three derived tables. Still flagged ([warnUnverified]),
     * since nothing here is checked against a reference file the way DNA, RNA and
     * protein masses are.
     *
     * A PNA monomer is N-(2-aminoethyl)glycine (C4H10N2O2) carrying its base on
     * the secondary nitrogen through a methylenecarbonyl linker, so the residue is
     * the backbone plus the base's acetic acid less two waters, one for the amide
     * to the base and one for the peptide bond in the chain: for adenine,
     * C4H10N2O2 + C7H7N5O2 - 2 H2O = C11H13N7O2.
     *
     * The backbone has no phosphate and no charge, and a strand ends in a free
     * amine and a free acid, so [pna] adds one water as [protein] does.
     */
    private val PNA = mapOf(
        'A' to f(11, 13, 7, 2), 'C' to f(10, 13, 5, 3), 'G' to f(11, 13, 7, 3), 'T' to f(11, 14, 4, 4),
    )

    /** Amino acid residues: the amino acid less H2O (Ala C3H5NO, 71.08 Da). */
    private val PROTEIN = mapOf(
        'A' to f(3, 5, 1, 1), 'R' to f(6, 12, 4, 1), 'N' to f(4, 6, 2, 2), 'D' to f(4, 5, 1, 3),
        'C' to f(3, 5, 1, 1, s = 1), 'E' to f(5, 7, 1, 3), 'Q' to f(5, 8, 2, 2), 'G' to f(2, 3, 1, 1),
        'H' to f(6, 7, 3, 1), 'I' to f(6, 11, 1, 1), 'L' to f(6, 11, 1, 1), 'K' to f(6, 12, 2, 1),
        'M' to f(5, 9, 1, 1, s = 1), 'F' to f(9, 9, 1, 1), 'P' to f(5, 7, 1, 1), 'S' to f(3, 5, 1, 2),
        'T' to f(4, 7, 1, 2), 'W' to f(11, 10, 2, 1), 'Y' to f(9, 9, 1, 2), 'V' to f(5, 9, 1, 1),
        'U' to f(3, 5, 1, 1, se = 1), 'O' to f(12, 19, 3, 2),
    )

    /**
     * One strand of nucleic acid: null when a letter is ambiguous. A linear strand
     * has one phosphate fewer than residues, whichever backbone carries them.
     */
    fun strand(letters: String, backbone: Backbone, circular: Boolean): IntArray? {
        val table = when (backbone) {
            Backbone.DNA -> DNA
            Backbone.RNA -> RNA
            Backbone.TNA -> TNA
            Backbone.LNA -> LNA
        }
        val sum = IntArray(SYMBOLS.size)
        for (ch in letters) add(sum, table[ch] ?: return null)
        if (!circular) { sum[4] -= 1; sum[3] -= 2; sum[1] += 1 }     // less HPO3, plus H2O: P-1, O-2, H+1
        return sum
    }

    /** A PNA strand: its residues plus one water, as a peptide. Null when a letter is ambiguous. */
    fun pna(letters: String): IntArray? {
        val sum = IntArray(SYMBOLS.size)
        for (ch in letters) add(sum, PNA[ch] ?: return null)
        sum[1] += 2; sum[3] += 1                                        // plus H2O
        return sum
    }

    private val warned = java.util.concurrent.ConcurrentHashMap<String, Boolean>()

    /**
     * Warns once per kind per JVM that a mass comes from a derived residue table
     * rather than a checked reference (TNA, PNA and LNA; see their tables above).
     * The line goes through Ki.Core's [io.kixi.warn], so RevNexus shows it in its
     * orange Messages pane and a plain JVM prints it to stdout.
     */
    fun warnUnverified(kind: String) {
        if (warned.putIfAbsent(kind, true) != null) return
        io.kixi.warn(
            "$kind masses are derived from the backbone's chemistry and are not verified against a " +
                "reference table, unlike DNA, RNA and protein masses. Treat ${Letters.article(kind)} $kind " +
                "mass, monoisotopic mass or formula as provisional and check it before relying on it."
        )
    }

    /** A protein: null when a letter has no formula (B Z J X) or a stop is inside it. A final stop is ignored. */
    fun protein(letters: String): IntArray? {
        val body = letters.removeSuffix("*")
        if (body.isEmpty()) return null                                  // a stop alone has no residues
        val sum = IntArray(SYMBOLS.size)
        for (ch in body) add(sum, PROTEIN[ch] ?: return null)
        sum[1] += 2; sum[3] += 1                                        // plus H2O
        return sum
    }

    fun plus(a: IntArray, b: IntArray): IntArray = IntArray(a.size) { a[it] + b[it] }

    fun formula(counts: IntArray): Formula {
        val m = LinkedHashMap<Element, Int>()
        for (i in counts.indices) if (counts[i] > 0) m[ELEMENTS[i]] = counts[i]
        return Formula.of(m)
    }

    private fun add(sum: IntArray, r: IntArray) { for (i in sum.indices) sum[i] = Math.addExact(sum[i], r[i]) }
}

/**
 * What a nucleic acid's residues hang from: the sugar and the linkage. It picks the
 * residue table in [Masses.strand], and `RNA` is also what makes a sequence read
 * and write U where the others read T.
 *
 * `TNA` and `LNA` masses are derived rather than checked; see the tables.
 */
internal enum class Backbone { DNA, RNA, TNA, LNA }
