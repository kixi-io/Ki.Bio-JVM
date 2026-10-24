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

    /** Amino acid residues: the amino acid less H2O (Ala C3H5NO, 71.08 Da). */
    private val PROTEIN = mapOf(
        'A' to f(3, 5, 1, 1), 'R' to f(6, 12, 4, 1), 'N' to f(4, 6, 2, 2), 'D' to f(4, 5, 1, 3),
        'C' to f(3, 5, 1, 1, s = 1), 'E' to f(5, 7, 1, 3), 'Q' to f(5, 8, 2, 2), 'G' to f(2, 3, 1, 1),
        'H' to f(6, 7, 3, 1), 'I' to f(6, 11, 1, 1), 'L' to f(6, 11, 1, 1), 'K' to f(6, 12, 2, 1),
        'M' to f(5, 9, 1, 1, s = 1), 'F' to f(9, 9, 1, 1), 'P' to f(5, 7, 1, 1), 'S' to f(3, 5, 1, 2),
        'T' to f(4, 7, 1, 2), 'W' to f(11, 10, 2, 1), 'Y' to f(9, 9, 1, 2), 'V' to f(5, 9, 1, 1),
        'U' to f(3, 5, 1, 1, se = 1), 'O' to f(12, 19, 3, 2),
    )

    /** One strand of nucleic acid: null when a letter is ambiguous. */
    fun strand(letters: String, rna: Boolean, circular: Boolean): IntArray? {
        val table = if (rna) RNA else DNA
        val sum = IntArray(SYMBOLS.size)
        for (ch in letters) add(sum, table[ch] ?: return null)
        if (!circular) { sum[4] -= 1; sum[3] -= 2; sum[1] += 1 }     // less HPO3, plus H2O: P-1, O-2, H+1
        return sum
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
