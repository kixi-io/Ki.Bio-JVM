package io.kixi.bio

/**
 * The letters of DNA, RNA and protein, and the IUPAC codes.
 *
 * A nucleotide letter is a set of bases, kept as a 4-bit mask (A 1, C 2, G 4,
 * T or U 8): `R` is A or G (5), `N` any (15). A sequence letter matches a pattern
 * letter when its bases are among the pattern's, so the pattern `N` matches
 * anything and the sequence letter `N` matches only `N`.
 */
internal object Letters {

    const val DNA: String = "ACGTRYSWKMBDHVN"
    const val RNA: String = "ACGURYSWKMBDHVN"

    /** The 20, selenocysteine (U), pyrrolysine (O), the ambiguity codes B Z J X, and a stop. */
    const val PROTEIN: String = "ACDEFGHIKLMNPQRSTVWYUOBZJX*"

    private val MASK = IntArray(128).also { m ->
        for ((c, v) in listOf('A' to 1, 'C' to 2, 'G' to 4, 'T' to 8, 'U' to 8, 'R' to 5, 'Y' to 10, 'S' to 6,
            'W' to 9, 'K' to 12, 'M' to 3, 'B' to 14, 'D' to 13, 'H' to 11, 'V' to 7, 'N' to 15)) m[c.code] = v
    }

    /** The DNA letter for each mask, 1 to 15 (0 is no letter). */
    private const val DNA_OF_MASK = "?ACMGRSVTWYHKDBN"

    fun mask(c: Char): Int = if (c.code < 128) MASK[c.code] else 0

    fun dnaLetter(mask: Int): Char = DNA_OF_MASK[mask]

    fun rnaLetter(mask: Int): Char = if (mask == 8) 'U' else DNA_OF_MASK[mask]

    /** The complement of a base set: A and T swap, C and G swap. */
    fun complementMask(m: Int): Int = ((m and 1) shl 3) or ((m and 8) ushr 3) or ((m and 2) shl 1) or ((m and 4) ushr 1)

    /** True for A, C, G, T and U. */
    fun isPlainBase(c: Char): Boolean = c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U'

    /** A nucleotide [s] of a sequence matches the pattern letter [p] when its bases are among [p]'s. */
    fun nucleotideMatches(s: Char, p: Char): Boolean {
        val ms = mask(s)
        val mp = mask(p)
        return ms != 0 && mp != 0 && (ms and mp.inv()) == 0
    }

    /** A residue [s] matches the pattern letter [p]: itself, X for any residue, B for D or N, Z for E or Q, J for I or L. */
    fun residueMatches(s: Char, p: Char): Boolean = s == p || when (p) {
        'X' -> s != '*'
        'B' -> s == 'D' || s == 'N'
        'Z' -> s == 'E' || s == 'Q'
        'J' -> s == 'I' || s == 'L'
        else -> false
    }

    /**
     * The letters of [text] uppercased, whitespace dropped, checked against
     * [alphabet]; [what] names the kind in errors ("DNA").
     */
    fun letters(text: String, alphabet: String, what: String): String {
        val sb = StringBuilder(text.length)
        for (ch in text) if (!ch.isWhitespace()) sb.append(ch.uppercaseChar())
        return checked(sb.toString(), alphabet, what)
    }

    fun checked(letters: String, alphabet: String, what: String): String {
        if (letters.isEmpty()) throw BioException("An empty $what sequence")
        var bad: LinkedHashSet<Char>? = null
        for (ch in letters) if (alphabet.indexOf(ch) < 0) (bad ?: LinkedHashSet<Char>().also { bad = it }).add(ch)
        bad?.let { b ->
            throw BioException("Not ${article(what)} $what sequence: unexpected ${b.take(8).joinToString(", ") { "'$it'" }}${if (b.size > 8) ", ..." else ""}")
        }
        return letters
    }

    /**
     * `a` or `an` before a kind's name, read as it is said: an initialism goes by
     * its first letter's sound, so it is "a DNA" and "a TNA" but "an RNA" and
     * "an LNA".
     */
    /** A kind's name in prose: an initialism keeps its capitals (`DNA`), a word goes lowercase (`protein`). */
    fun prose(kind: String): String = if (kind == kind.uppercase()) kind else kind.lowercase()

    fun article(what: String): String {
        val first = what.first()
        val vowelSound = if (what == what.uppercase()) first in "AEFHILMNORSX" else first.lowercaseChar() in "aeiou"
        return if (vowelSound) "an" else "a"
    }

    /**
     * G + C as a percentage of the bases known to be one or the other (G, C and S
     * against A, T, U and W; N, R, Y and the other codes are left out, as
     * Biopython's `gc_fraction` does by default), unrounded; null when there are
     * none.
     */
    fun gc(letters: String): Double? {
        val strong = letters.count { it == 'G' || it == 'C' || it == 'S' }
        val weak = letters.count { it == 'A' || it == 'T' || it == 'U' || it == 'W' }
        return if (strong + weak == 0) null else strong * 100.0 / (strong + weak)
    }

    /**
     * Text as pasted from anywhere: a FASTA header line (`>` or `;`), whitespace
     * and digits (GenBank numbering) dropped, uppercased.
     */
    fun clean(text: String): String {
        val sb = StringBuilder(text.length)
        for (line in text.lineSequence()) {
            val t = line.trim()
            if (t.startsWith(">") || t.startsWith(";")) continue
            for (ch in t) if (!ch.isWhitespace() && !ch.isDigit()) sb.append(ch.uppercaseChar())
        }
        return sb.toString()
    }

    /**
     * The kind of cleaned letters when none is given: RNA when only A, C, G, U and N
     * with a U; DNA when only A, C, G, T and N; protein otherwise. Every IUPAC code
     * is an amino acid letter too, so a sequence with them is DNA or RNA only when
     * declared.
     */
    fun detect(cleaned: String): Kind = when {
        cleaned.all { it in "ACGUN" } && 'U' in cleaned -> Kind.RNA
        cleaned.all { it in "ACGTN" } -> Kind.DNA
        else -> Kind.PROTEIN
    }

    enum class Kind { DNA, RNA, PROTEIN }
}
