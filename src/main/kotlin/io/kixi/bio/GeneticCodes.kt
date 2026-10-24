package io.kixi.bio

/**
 * One of NCBI's genetic codes: how codons read as amino acids, which start and
 * which stop. [id] is NCBI's number (1 the standard code, 2 vertebrate
 * mitochondrial, 11 bacterial and plastid, ...).
 */
public class GeneticCode internal constructor(
    public val id: Int,
    public val name: String,
    private val aminoAcids: String,
    private val starts: String,
    private val stops: String,
) {
    /**
     * The amino acid a codon reads as (`*` for a stop). An ambiguous codon (`GGN`)
     * reads as the amino acid every codon it stands for agrees on, else `X`.
     * T and U are the same.
     */
    public fun translate(codon: CharSequence): Char {
        require(codon.length == 3) { "A codon is three letters, not '$codon'" }
        val i = index(codon)
        if (i >= 0) return aminoAcids[i]
        val m0 = Letters.mask(codon[0]); val m1 = Letters.mask(codon[1]); val m2 = Letters.mask(codon[2])
        if (m0 == 0 || m1 == 0 || m2 == 0) return 'X'
        var found = 0.toChar()
        for (a in 0..3) if ((m0 shr BASE_BITS[a]) and 1 == 1)
            for (b in 0..3) if ((m1 shr BASE_BITS[b]) and 1 == 1)
                for (c in 0..3) if ((m2 shr BASE_BITS[c]) and 1 == 1) {
                    val aa = aminoAcids[a * 16 + b * 4 + c]
                    if (found == 0.toChar()) found = aa else if (found != aa) return 'X'
                }
        return found
    }

    /** True when [codon] (plain bases) is a stop codon of this code. */
    public fun isStop(codon: CharSequence): Boolean = index(codon).let { it >= 0 && stops[it] == '*' }

    /** True when [codon] (plain bases) is one of this code's start codons. */
    public fun isStart(codon: CharSequence): Boolean = index(codon).let { it >= 0 && starts[it] == 'M' }

    override fun toString(): String = "$id $name"

    private fun index(codon: CharSequence): Int {
        if (codon.length != 3) return -1
        var i = 0
        for (k in 0..2) {
            val b = when (codon[k]) { 'T', 'U' -> 0; 'C' -> 1; 'A' -> 2; 'G' -> 3; else -> return -1 }
            i = i * 4 + b
        }
        return i
    }

    private companion object {
        /** TCAG, the order of NCBI's tables, as bits of a base mask (A 1, C 2, G 4, T 8). */
        val BASE_BITS = intArrayOf(3, 1, 0, 2)
    }
}

/**
 * NCBI's 27 genetic codes (1 to 6, 9 to 16, 21 to 33), from `genetic-codes.txt`,
 * which is generated from NCBI's `gc.prt` (by way of Biopython's copy).
 */
public object GeneticCodes {

    public val all: List<GeneticCode> by lazy {
        val text = GeneticCodes::class.java.getResourceAsStream("genetic-codes.txt")?.bufferedReader()?.readText()
            ?: error("genetic-codes.txt is missing from Ki.Bio's resources")
        text.lineSequence().filter { it.isNotBlank() && !it.startsWith("#") }.map { line ->
            val f = line.split('\t')
            check(f.size == 5 && f[2].length == 64 && f[3].length == 64 && f[4].length == 64) { "Bad genetic code line: $line" }
            GeneticCode(f[0].toInt(), f[1], f[2], f[3], f[4])
        }.toList()
    }

    private val byId: Map<Int, GeneticCode> by lazy { all.associateBy { it.id } }

    /** The standard code, NCBI 1. */
    public val standard: GeneticCode get() = byId.getValue(1)

    /** The code numbered [id], or null. */
    @JvmStatic
    public fun byId(id: Int): GeneticCode? = byId[id]

    /** The code numbered [id] (an Int, or a whole number), or a [BioException] naming the ids there are. */
    @JvmStatic
    public fun of(id: Any?): GeneticCode {
        val n = when (id) {
            null -> return standard
            is GeneticCode -> return id
            is Number -> id.toDouble().takeIf { it == Math.rint(it) }?.toInt()
            else -> null
        } ?: throw BioException("A genetic code is an NCBI table number, not $id")
        return byId(n) ?: throw BioException("There is no NCBI genetic code $n (there are ${all.joinToString(", ") { it.id.toString() }})")
    }
}
