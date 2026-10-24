package io.kixi.bio

import io.kixi.sci.Formula
import io.kixi.uom.Concentration
import io.kixi.uom.MolarMass
import io.kixi.uom.Quantity
import io.kixi.uom.Temperature
import java.math.BigDecimal
import java.security.MessageDigest
import java.util.Base64
import io.kixi.uom.Unit as KiUnit

/**
 * A DNA, RNA or protein sequence (Ki type bags design, 5.5): its letters, kept as
 * written and uppercased. The schema and KS type `Sequence` is any of the three.
 *
 * Positions are 1-based, as in GenBank and every lab tool.
 */
public sealed class BioSequence(letters: String) {

    /** The letters, uppercase. */
    public val letters: String = letters

    /** `DNA`, `RNA` or `Protein`. */
    public abstract val kind: String

    /** Residues: nucleotides, base pairs or amino acids. */
    public val length: Int get() = letters.length

    /** True for a circular molecule (a plasmid). Proteins are linear. */
    public open val circular: Boolean get() = false

    /** True when [s] (a letter of this sequence) matches the pattern letter [p]. */
    internal abstract fun matches(s: Char, p: Char): Boolean

    /** The letters a pattern may use. */
    internal abstract val patternAlphabet: String

    /** A sequence like this one with other letters (and the same flags unless given). */
    internal abstract fun withLetters(letters: String, circular: Boolean = this.circular): BioSequence

    /**
     * The 1-based start of every match of [pattern] (text or a sequence),
     * overlapping matches included; across the origin when circular. The pattern
     * may use IUPAC codes: `GAANNNNTTC` is a site.
     */
    public fun find(pattern: Any?): List<Int> = search(patternLetters(pattern), firstOnly = false)

    /** How many times [pattern] occurs, overlaps included. */
    public fun count(pattern: Any?): Int = find(pattern).size

    /** True when [pattern] occurs: what KS and KQ `pattern in seq` ask. */
    public fun contains(pattern: Any?): Boolean = search(patternLetters(pattern), firstOnly = true).isNotEmpty()

    private fun search(p: String, firstOnly: Boolean): List<Int> {
        val n = letters.length
        val m = p.length
        if (m > n) return emptyList()
        val out = ArrayList<Int>()
        val last = if (circular) n - 1 else n - m
        for (i in 0..last) {
            var ok = true
            for (j in 0 until m) {
                val k = i + j
                if (!matches(letters[if (k < n) k else k - n], p[j])) { ok = false; break }
            }
            if (ok) { out.add(i + 1); if (firstOnly) break }
        }
        return out
    }

    /**
     * The letters [from] to [to], 1-based and inclusive. On a circular sequence
     * [from] after [to] runs across the origin. The result is linear.
     */
    public fun slice(from: Int, to: Int): BioSequence {
        val n = letters.length
        if (from !in 1..n || to !in 1..n) throw BioException("A slice of a $n-residue sequence runs from 1 to $n, not $from to $to")
        val s = when {
            from <= to -> letters.substring(from - 1, to)
            circular -> letters.substring(from - 1) + letters.substring(0, to)
            else -> throw BioException("A slice from $from to $to runs backwards; only a circular sequence wraps")
        }
        return withLetters(s, circular = false)
    }

    /** The letters reversed (not complemented). */
    public fun reverse(): BioSequence = withLetters(letters.reversed())

    /** The chemical formula, or null when a letter is ambiguous (or a stop is inside a protein). */
    public abstract val formula: Formula?

    /** The molar mass in Da from [formula] and the standard atomic weights, or null when there is none. */
    public val mass: Quantity<MolarMass>? get() = formula?.mass

    /** The monoisotopic mass in Da, or null when there is no [formula]. */
    public val monoisotopicMass: Quantity<MolarMass>? get() = formula?.monoisotopicMass

    /**
     * The GA4GH refget identifier: `SQ.` and the base64url of the first 24 bytes
     * of the letters' SHA-512. The same letters give the same digest whatever the
     * kind; index it for long sequences, whose own keys are as long as they are.
     */
    public val digest: String by lazy {
        val h = MessageDigest.getInstance("SHA-512").digest(letters.toByteArray(Charsets.US_ASCII))
        "SQ." + Base64.getUrlEncoder().withoutPadding().encodeToString(h.copyOf(24))
    }

    /**
     * True when [other] is the same molecule: the same kind and letters, or, for a
     * circular sequence, the letters from another origin, or, for a double-stranded
     * one, the other strand.
     */
    public abstract fun sameMolecule(other: Any?): Boolean

    private fun patternLetters(pattern: Any?): String {
        val p = when (pattern) {
            is BioSequence -> {
                if ((pattern is Protein) != (this is Protein)) throw BioException("Cannot look for a ${pattern.kind} in a $kind")
                pattern.letters
            }
            is CharSequence -> pattern.filterNot { it.isWhitespace() }.toString().uppercase()
            null -> throw BioException("A pattern to find, not nil")
            else -> throw BioException("A pattern is text or a sequence, not ${pattern::class.simpleName}")
        }
        if (p.isEmpty()) throw BioException("An empty pattern")
        for (ch in p) if (patternAlphabet.indexOf(ch) < 0) throw BioException("'$ch' cannot be in a $kind pattern")
        return p
    }

    override fun equals(other: Any?): Boolean =
        other is BioSequence && other.javaClass == javaClass && other.letters == letters && other.flags == flags

    override fun hashCode(): Int = (javaClass.hashCode() * 31 + letters.hashCode()) * 31 + flags

    /** The letters. */
    override fun toString(): String = letters

    /** Bit 0 double-stranded, bit 1 circular. */
    internal open val flags: Int get() = 0
}

/** DNA or RNA: a strand, or a duplex given by its top strand ([ds]), linear or [circular]. */
public sealed class NucleicAcid(letters: String, public val ds: Boolean, final override val circular: Boolean) : BioSequence(letters) {

    override val flags: Int get() = (if (ds) 1 else 0) or (if (circular) 2 else 0)

    override fun matches(s: Char, p: Char): Boolean = Letters.nucleotideMatches(s, p)

    override val patternAlphabet: String get() = "ACGTURYSWKMBDHVN"

    internal abstract fun with(letters: String, ds: Boolean = this.ds, circular: Boolean = this.circular): NucleicAcid

    override fun withLetters(letters: String, circular: Boolean): BioSequence = with(letters, circular = circular)

    internal abstract val isRna: Boolean

    /**
     * G + C as a percentage of the bases that are known to be one or the other
     * (G, C and S against A, T, U and W; N, R, Y and the other codes are left out,
     * as Biopython's `gc_fraction` does by default), unrounded; null when there are
     * none.
     */
    public val gc: Double? get() {
        val strong = letters.count { it == 'G' || it == 'C' || it == 'S' }
        val weak = letters.count { it == 'A' || it == 'T' || it == 'U' || it == 'W' }
        return if (strong + weak == 0) null else strong * 100.0 / (strong + weak)
    }

    /** The complementary strand, same direction; IUPAC codes complement too (R and Y, ...). */
    public fun complement(): NucleicAcid = with(complementLetters(letters))

    /** The other strand read 5' to 3'. */
    public fun reverseComplement(): NucleicAcid = with(complementLetters(letters).reversed())

    private fun complementLetters(s: String): String {
        val sb = StringBuilder(s.length)
        for (ch in s) {
            val m = Letters.complementMask(Letters.mask(ch))
            sb.append(if (isRna) Letters.rnaLetter(m) else Letters.dnaLetter(m))
        }
        return sb.toString()
    }

    /**
     * The codons of reading [frame]: 1, 2 or 3 on this strand, -1, -2 or -3 on the
     * other. A partial codon at the end is dropped.
     */
    public fun codons(frame: Int = 1): List<String> {
        val strand = strand(frame)
        val out = ArrayList<String>(strand.length / 3)
        var i = Math.abs(frame) - 1
        while (i + 3 <= strand.length) { out.add(strand.substring(i, i + 3)); i += 3 }
        return if (isRna) out.map { it.replace('T', 'U') } else out
    }

    /** The protein reading [frame] gives with the genetic [code] (an NCBI table number, 1 by default); stops are `*`. */
    public fun translate(frame: Int = 1, code: Any? = null): Protein {
        val gc = GeneticCodes.of(code)
        val strand = strand(frame)
        val sb = StringBuilder(strand.length / 3)
        var i = Math.abs(frame) - 1
        while (i + 3 <= strand.length) { sb.append(gc.translate(strand.subSequence(i, i + 3))); i += 3 }
        if (sb.isEmpty()) throw BioException("Frame $frame of a ${letters.length}-nt sequence has no whole codon")
        return Protein(sb.toString())
    }

    /**
     * Every open reading frame of at least [minLength] nucleotides (stop codon
     * included) on either strand, longest first: from an ATG to the next stop in
     * frame (stops from the genetic [code]), the longest one for each stop. On a
     * circular sequence frames run across the origin and an ORF is at most the
     * circle's length. Each is a linear, single strand.
     */
    public fun orfs(minLength: Int = 100, code: Any? = null): List<NucleicAcid> {
        if (minLength < 3) throw BioException("An ORF is at least 3 nt, not $minLength")
        val gc = GeneticCodes.of(code)
        val n = letters.length
        if (n < 6) return emptyList()
        class Orf(val strand: Int, val start: Int, val text: String)
        val best = HashMap<Long, Orf>()   // the longest ORF for each strand and stop codon
        for ((si, s0) in listOf(letters, reverseComplement().letters).withIndex()) {
            val s = if (isRna) s0.replace('U', 'T') else s0
            val ring = if (circular) s + s.substring(0, 2) else s
            fun codon(i: Int): CharSequence? = if (i + 3 <= ring.length) ring.subSequence(i, i + 3) else null
            val stop = BooleanArray(n) { i -> codon(i)?.let { gc.isStop(it) } ?: false }
            // toStop[i]: nucleotides from i to the first in-frame stop codon at or after i, or -1.
            val toStop = IntArray(n) { -1 }
            if (circular) {
                // i -> i + 3 (mod n) runs around one cycle, or three when 3 divides n.
                val seen = BooleanArray(n)
                for (first in 0 until n) {
                    if (seen[first]) continue
                    val cycle = ArrayList<Int>()
                    var i = first
                    while (!seen[i]) { seen[i] = true; cycle.add(i); i = (i + 3) % n }
                    val anchor = cycle.indexOfFirst { stop[it] }
                    if (anchor < 0) continue
                    for (t in 0 until cycle.size) {
                        val j = Math.floorMod(anchor - t, cycle.size)
                        val pos = cycle[j]
                        toStop[pos] = if (stop[pos]) 0 else toStop[cycle[(j + 1) % cycle.size]] + 3
                    }
                }
            } else {
                for (i in n - 1 downTo 0) toStop[i] = if (stop[i]) 0 else if (i + 3 < n && toStop[i + 3] >= 0) toStop[i + 3] + 3 else -1
            }
            val doubled = if (circular) s + s else s
            for (i in 0 until n) {
                if (toStop[i] < 0 || codon(i)?.toString() != "ATG") continue
                val len = toStop[i] + 3
                if (circular && len > n) continue
                val key = si.toLong() * (n.toLong() + 1) + (i + toStop[i]) % n
                val prev = best[key]
                if (prev == null || prev.text.length < len) best[key] = Orf(si, i, doubled.substring(i, i + len))
            }
        }
        return best.values.filter { it.text.length >= minLength }
            .sortedWith(compareByDescending<Orf> { it.text.length }.thenBy { it.strand }.thenBy { it.start })
            .map { with(if (isRna) it.text.replace('T', 'U') else it.text, ds = false, circular = false) }
    }

    /** The longest ORF on either strand (at least one codon and a stop), or a [BioException]. */
    public fun orf(code: Any? = null): NucleicAcid =
        orfs(3, code).firstOrNull() ?: throw BioException("No open reading frame (ATG to a stop) on either strand")

    private fun strand(frame: Int): String {
        if (frame == 0 || frame > 3 || frame < -3) throw BioException("A frame is 1, 2, 3, -1, -2 or -3, not $frame")
        val s = if (frame > 0) letters else reverseComplement().letters
        return if (isRna) s.replace('U', 'T') else s
    }

    override val formula: Formula? by lazy {
        val top = Masses.strand(letters, isRna, circular) ?: return@lazy null
        val counts = if (ds) Masses.plus(top, Masses.strand(reverseComplement().letters, isRna, circular) ?: return@lazy null) else top
        Masses.formula(counts)
    }

    override fun sameMolecule(other: Any?): Boolean {
        if (other !is NucleicAcid || other.javaClass != javaClass || other.ds != ds || other.circular != circular) return false
        if (other.letters.length != letters.length) return false
        val candidates = if (ds) listOf(other.letters, other.reverseComplement().letters) else listOf(other.letters)
        return candidates.any { c -> if (circular) (letters + letters).contains(c) else c == letters }
    }
}

/**
 * DNA: `.dna("GAATTC")`, `.dna(ds "…")` for a duplex, `.dna(ds circular "…")` for a
 * plasmid. Letters: A C G T and the IUPAC codes.
 */
public class DNA(letters: String, ds: Boolean = false, circular: Boolean = false) :
    NucleicAcid(Letters.letters(letters, Letters.DNA, "DNA"), ds, circular) {

    override val kind: String get() = "DNA"
    override val isRna: Boolean get() = false
    override fun with(letters: String, ds: Boolean, circular: Boolean): DNA = DNA(letters, ds, circular)

    /** The RNA this DNA is transcribed into: T becomes U. */
    public fun transcribe(): RNA = RNA(letters.replace('T', 'U'), ds, circular)

    /**
     * The melting temperature in °C with the complement, as IDT's OligoAnalyzer
     * computes it (see [Melting]): [na] monovalent cations (default 50mM), [oligo]
     * the oligo (250nM), [mg] Mg2+ (0), [dntp] dNTPs (0); the order RevGrid's
     * `tm(dna, salt, oligo)` always had, with the buffer's Mg2+ and dNTPs after.
     * Each is a concentration quantity or a number in mol/L. Plain A, C, G and T only.
     */
    public fun tm(na: Any? = null, oligo: Any? = null, mg: Any? = null, dntp: Any? = null): Quantity<Temperature> {
        if (circular) throw BioException("A melting temperature is for a linear oligo, not a circular sequence")
        if (letters.length < 2) throw BioException("A melting temperature needs at least 2 nt")
        letters.firstOrNull { it !in "ACGT" }?.let { throw BioException("A melting temperature needs plain A, C, G and T; this has '$it'") }
        val naM = molar(na, 0.05, "Na+")
        val oligoM = molar(oligo, 250e-9, "The oligo concentration")
        val mgM = molar(mg, 0.0, "Mg2+")
        val dntpM = molar(dntp, 0.0, "dNTP")
        if (oligoM == 0.0) throw BioException("The oligo concentration must be more than 0")
        if (naM == 0.0 && mgM == 0.0) throw BioException("A melting temperature needs Na+ or Mg2+ above 0")
        val celsius = Melting.tm(letters, naM, mgM, dntpM, oligoM)
        if (!celsius.isFinite()) throw BioException("No melting temperature for these conditions (Na+ $na, oligo $oligo, Mg2+ $mg, dNTPs $dntp)")
        return Quantity(BigDecimal(celsius.toString()), KiUnit.dC)
    }

    public companion object {
        /** A concentration in mol/L from a quantity or a number, or [default] for nil. */
        internal fun molar(x: Any?, default: Double, what: String): Double {
            val v = when (x) {
                null -> default
                is Quantity<*> -> {
                    if (x.unit !is Concentration) throw BioException("$what is a concentration such as 50mM, not ${x}")
                    @Suppress("UNCHECKED_CAST")
                    ((x as Quantity<Concentration>) convertTo KiUnit.M).value.toDouble()
                }
                is Number -> x.toDouble()
                else -> throw BioException("$what is a concentration such as 50mM, not ${x::class.simpleName}")
            }
            if (!v.isFinite() || v < 0) throw BioException("$what cannot be $x")
            return v
        }
    }
}

/** RNA: `.rna("AUGGCC")` (also `ds`, `circular`). Letters: A C G U and the IUPAC codes. */
public class RNA(letters: String, ds: Boolean = false, circular: Boolean = false) :
    NucleicAcid(Letters.letters(letters, Letters.RNA, "RNA"), ds, circular) {

    override val kind: String get() = "RNA"
    override val isRna: Boolean get() = true
    override fun with(letters: String, ds: Boolean, circular: Boolean): RNA = RNA(letters, ds, circular)

    /** The DNA this RNA is reverse transcribed into (the same strand, U becomes T). */
    public fun reverseTranscribe(): DNA = DNA(letters.replace('U', 'T'), ds, circular)
}

/**
 * A protein: `.protein("MKTAYIAKQR")`. Letters: the 20, U (selenocysteine),
 * O (pyrrolysine), the ambiguity codes B Z J X, and `*` for a stop.
 */
public class Protein(letters: String) : BioSequence(Letters.letters(letters, Letters.PROTEIN, "protein")) {

    override val kind: String get() = "Protein"
    override fun matches(s: Char, p: Char): Boolean = Letters.residueMatches(s, p)
    override val patternAlphabet: String get() = Letters.PROTEIN
    override fun withLetters(letters: String, circular: Boolean): Protein = Protein(letters)

    override val formula: Formula? by lazy { Masses.protein(letters)?.let { Masses.formula(it) } }

    /** The residues that count for charge: the letters less a final stop; null when a stop is inside. */
    private val chain: String? get() = letters.removeSuffix("*").takeIf { it.isNotEmpty() && '*' !in it }

    /**
     * The isoelectric point, as ExPASy's Compute pI/Mw finds it (Bjellqvist et al.
     * 1993, 1994): side-chain pKs, the N-terminal pK of the first residue and the
     * C-terminal pK of the last, bisected between pH 0 and 14 to 0.0001. Null when
     * a stop is inside the protein.
     */
    public val pI: Double? by lazy {
        val c = chain ?: return@lazy null
        var lo = 0.0
        var hi = 14.0
        var loops = 0
        while (hi - lo >= 0.0001 && loops++ < 2000) {
            val mid = (lo + hi) / 2
            if (charge(c, mid) > 0) lo = mid else hi = mid
        }
        (lo + hi) / 2
    }

    /** The net charge at [pH] (Bjellqvist pKs, as [pI]); null when a stop is inside. */
    public fun charge(pH: Double): Double? = chain?.let { charge(it, pH) }

    /**
     * The molar extinction coefficient at 280 nm in water, M⁻¹cm⁻¹ (Pace et al.
     * 1995): 5500 per Trp, 1490 per Tyr, 125 per cystine, every pair of Cys taken
     * as a cystine. ExPASy ProtParam's first value.
     */
    public val extinction280: Int get() = 5500 * letters.count { it == 'W' } + 1490 * letters.count { it == 'Y' } + 125 * (letters.count { it == 'C' } / 2)

    /** [extinction280] with every Cys reduced (no cystines): ProtParam's second value. */
    public val extinction280Reduced: Int get() = 5500 * letters.count { it == 'W' } + 1490 * letters.count { it == 'Y' }

    override fun sameMolecule(other: Any?): Boolean = other is Protein && other.letters == letters

    private companion object {
        fun charge(c: String, ph: Double): Double {
            fun pos(pk: Double) = 1.0 / (1.0 + Math.pow(10.0, ph - pk))
            fun neg(pk: Double) = 1.0 / (1.0 + Math.pow(10.0, pk - ph))
            val nTerm = when (c.first()) { 'A' -> 7.59; 'M' -> 7.0; 'S' -> 6.93; 'P' -> 8.36; 'T' -> 6.82; 'V' -> 7.44; 'E' -> 7.7; else -> 7.5 }
            val cTerm = when (c.last()) { 'D' -> 4.55; 'E' -> 4.75; else -> 3.55 }
            var k = 0; var r = 0; var h = 0; var d = 0; var e = 0; var cys = 0; var y = 0
            for (ch in c) when (ch) { 'K' -> k++; 'R' -> r++; 'H' -> h++; 'D' -> d++; 'E' -> e++; 'C' -> cys++; 'Y' -> y++ }
            return pos(nTerm) + k * pos(10.0) + r * pos(12.0) + h * pos(5.98) -
                    neg(cTerm) - d * neg(4.05) - e * neg(4.45) - cys * neg(9.0) - y * neg(10.0)
        }
    }
}
