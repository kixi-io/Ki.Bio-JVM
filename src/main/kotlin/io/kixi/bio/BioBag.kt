package io.kixi.bio

import io.kixi.text.ParseException
import io.kixi.types.BagFunction
import io.kixi.types.Member
import io.kixi.types.TypeBag
import io.kixi.types.TypeBags
import io.kixi.types.ValueType
import kotlin.reflect.KClass

/**
 * The `bio` type bag: DNA, RNA and protein sequences for KD, KS, KQ and RevDB (Ki
 * type bags design, 5.5 and 11).
 *
 * - `DNA` `.dna("GAATTC")`, `.dna(ds circular "…")`; `RNA` `.rna("AUG")`; `Protein`
 *   `.protein("MKT")`; `Sequence`, any of the three, for schemas and `is`.
 * - Functions `dna`, `rna`, `protein`, `seq`, `revcomp`, `transcribe`, `translate`,
 *   `gc`, `tm`, `pI`, `a260`, `hamming`, `levenshtein`.
 *
 * The bag is found on the classpath through `META-INF/services`. It depends on
 * Ki.Sci for formulas and masses. Its Kotlin API is [DNA], [RNA], [Protein],
 * [BioSequence] and [GeneticCodes].
 */
public class BioBag : TypeBag {
    override val name: String get() = NAME
    override val version: Int get() = VERSION
    override val types: List<ValueType<*>> get() = TYPES
    override val functions: List<BagFunction> get() = FUNCTIONS

    public companion object {
        public const val NAME: String = "bio"
        public const val VERSION: Int = 1

        private val TYPES: List<ValueType<*>> = listOf(SequenceType, DNAType, RNAType, ProteinType)

        private val FUNCTIONS: List<BagFunction> = listOf(
            BagFunction("dna", 1..1, "DNA") { a -> dnaOf(a[0]) },
            BagFunction("rna", 1..1, "RNA") { a -> rnaOf(a[0]) },
            BagFunction("protein", 1..1, "Protein") { a -> proteinOf(a[0]) },
            BagFunction("seq", 1..2, "Sequence") { a -> seqOf(a[0], a.getOrNull(1)) },
            BagFunction("revcomp", 1..1, "Sequence") { a -> nucleicOf(a[0], "revcomp").reverseComplement() },
            BagFunction("transcribe", 1..1, "RNA") { a -> dnaOf(a[0]).transcribe() },
            BagFunction("translate", 1..3, "Protein") { a -> nucleicOf(a[0], "translate").translate(frameArg(a.getOrNull(1)), a.getOrNull(2)) },
            BagFunction("gc", 1..1, "Double") { a -> nucleicOf(a[0], "gc").gc },
            BagFunction("tm", 1..5, "Quantity") { a -> dnaOf(a[0]).tm(a.getOrNull(1), a.getOrNull(2), a.getOrNull(3), a.getOrNull(4)) },
            BagFunction("pI", 1..1, "Double") { a -> proteinOf(a[0]).pI },
            BagFunction("a260", 2..2, "Double") { a -> a260(a[0], a[1]) },
            BagFunction("hamming", 2..2, "Int") { a -> hamming(lettersOf(a[0]), lettersOf(a[1])) },
            BagFunction("levenshtein", 2..2, "Int") { a -> levenshtein(lettersOf(a[0]), lettersOf(a[1])) },
        )

        private fun name(x: Any?): String = x?.let { it::class.simpleName } ?: "nil"

        /** DNA from text (a FASTA header, whitespace and numbering dropped) or DNA. */
        @JvmStatic
        public fun dnaOf(x: Any?): DNA = when (x) {
            is DNA -> x
            is CharSequence -> DNA(Letters.clean(x.toString()))
            else -> throw BioException("Expected DNA, got ${name(x)}")
        }

        /** RNA from text or RNA. */
        @JvmStatic
        public fun rnaOf(x: Any?): RNA = when (x) {
            is RNA -> x
            is CharSequence -> RNA(Letters.clean(x.toString()))
            else -> throw BioException("Expected RNA, got ${name(x)}")
        }

        /** A protein from text or a protein. */
        @JvmStatic
        public fun proteinOf(x: Any?): Protein = when (x) {
            is Protein -> x
            is CharSequence -> Protein(Letters.clean(x.toString()))
            else -> throw BioException("Expected a protein, got ${name(x)}")
        }

        /**
         * A sequence from text, its kind detected unless [kind] names one (`dna`,
         * `dsdna`, `rna`, `dsrna`, `protein`; also `ssdna`, `peptide`, `aa`), or a
         * sequence (read as [kind] when one is given).
         */
        @JvmStatic
        public fun seqOf(x: Any?, kind: Any? = null): BioSequence {
            val letters = when (x) {
                is BioSequence -> if (kind == null) return x else x.letters
                is CharSequence -> Letters.clean(x.toString())
                else -> throw BioException("Expected a sequence, got ${name(x)}")
            }
            if (letters.isEmpty()) throw BioException("An empty sequence")
            return when (kind?.toString()?.trim()?.lowercase()?.replace("-", "")) {
                null -> when (Letters.detect(letters)) {
                    Letters.Kind.DNA -> DNA(letters)
                    Letters.Kind.RNA -> RNA(letters)
                    Letters.Kind.PROTEIN -> Protein(letters)
                }
                "dna", "ssdna" -> DNA(letters)
                "dsdna" -> DNA(letters, ds = true)
                "rna", "ssrna" -> RNA(letters)
                "dsrna" -> RNA(letters, ds = true)
                "protein", "peptide", "aa", "amino" -> Protein(letters)
                else -> throw BioException("Unknown sequence kind '$kind' (dna, dsdna, rna, dsrna, protein)")
            }
        }

        /** DNA or RNA from a sequence or text; [what] names the operation in errors. */
        @JvmStatic
        public fun nucleicOf(x: Any?, what: String): NucleicAcid = when (val s = seqOf(x)) {
            is NucleicAcid -> s
            else -> throw BioException("$what is for DNA and RNA, not a protein")
        }

        /**
         * Nucleic acid concentration in ng/µL from an A260 reading: 50 per unit for
         * double-stranded DNA, 33 for single-stranded DNA, 40 for RNA. [kind] is a
         * sequence or `dsdna`, `ssdna` (or `dna`), `rna`.
         */
        @JvmStatic
        public fun a260(absorbance: Any?, kind: Any?): Double {
            val a = (absorbance as? Number)?.toDouble() ?: throw BioException("An absorbance is a number, not ${name(absorbance)}")
            val factor = when (kind) {
                is DNA -> if (kind.ds) 50.0 else 33.0
                is RNA -> 40.0
                is Protein -> throw BioException("A260 conversion is for nucleic acids; use extinction280 for a protein")
                is CharSequence -> when (kind.toString().trim().lowercase().replace("-", "")) {
                    "dsdna" -> 50.0
                    "dna", "ssdna" -> 33.0
                    "rna", "ssrna" -> 40.0
                    else -> throw BioException("A260 conversion needs dsdna, ssdna or rna, not '$kind'")
                }
                else -> throw BioException("A260 conversion needs dsdna, ssdna or rna, not ${name(kind)}")
            }
            return a * factor
        }

        /** Positions at which two equally long sequences differ. */
        @JvmStatic
        public fun hamming(a: String, b: String): Int {
            if (a.length != b.length) throw BioException("Hamming distance needs equal lengths, not ${a.length} and ${b.length}")
            var d = 0
            for (i in a.indices) if (a[i] != b[i]) d++
            return d
        }

        /** Insertions, deletions and substitutions that turn [a] into [b]. Up to 10^8 cells (a 10 kb pair). */
        @JvmStatic
        public fun levenshtein(a: String, b: String): Int {
            if (a.length.toLong() * b.length > 100_000_000L) throw BioException("Levenshtein distance of ${a.length} by ${b.length} letters is too large (10^8 cells at most)")
            var prev = IntArray(b.length + 1) { it }
            var cur = IntArray(b.length + 1)
            for (i in 1..a.length) {
                cur[0] = i
                for (j in 1..b.length) {
                    val cost = if (a[i - 1] == b[j - 1]) 0 else 1
                    cur[j] = minOf(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost)
                }
                val t = prev; prev = cur; cur = t
            }
            return prev[b.length]
        }

        private fun lettersOf(x: Any?): String = when (x) {
            is BioSequence -> x.letters
            is CharSequence -> Letters.clean(x.toString())
            else -> throw BioException("Expected a sequence or text, got ${name(x)}")
        }

        internal fun intArg(x: Any?, what: String): Int = when (x) {
            is Number -> x.toDouble().takeIf { it == Math.rint(it) && it >= Int.MIN_VALUE && it <= Int.MAX_VALUE }?.toInt()
            else -> null
        } ?: throw BioException("$what is a whole number, not ${x ?: "nil"}")

        internal fun frameArg(x: Any?): Int = if (x == null) 1 else intArg(x, "A frame")
    }
}

/** Reads the inside of `.dna(…)`: `ds` and `circular` before a quoted string, or a bare token. */
private class LiteralParts(literal: String, allowFlags: Boolean, what: String) {
    val ds: Boolean
    val circular: Boolean
    val text: String

    init {
        var rest = TypeBags.inside(literal)
        var ds = false
        var circular = false
        while (true) {
            val m = FLAG.find(rest) ?: break
            val word = m.groupValues[1]
            if (!allowFlags) throw ParseException("A $what has no '$word'")
            when (word) {
                "ds" -> { if (ds) throw ParseException("'ds' twice in $literal"); ds = true }
                "circular" -> { if (circular) throw ParseException("'circular' twice in $literal"); circular = true }
            }
            rest = rest.substring(m.range.last + 1)
        }
        this.ds = ds
        this.circular = circular
        text = TypeBags.stringOrToken(rest)
    }

    private companion object { val FLAG = Regex("""^\s*(ds|circular)\s+""") }
}

/** Reads a literal, turning Ki.Bio's exceptions into Ki.Core's ParseException. */
private inline fun <T> readLiteral(read: () -> T): T = try { read() } catch (e: BioException) { throw ParseException(e.message) }

private fun nucleicLiteral(name: String, s: NucleicAcid): String {
    val flags = (if (s.ds) "ds " else "") + (if (s.circular) "circular " else "")
    return ".$name($flags\"${s.letters}\")"
}

// Members shared by the four types, typed for each.
private fun <T : BioSequence> common(): List<Member<T>> = listOf(
    Member("text", null, "String") { s, _ -> s.letters },
    Member("kind", null, "String") { s, _ -> s.kind },
    Member("length", null, "Int") { s, _ -> s.length },
    Member("formula", null, "Formula") { s, _ -> s.formula },
    Member("mass", null, "Quantity") { s, _ -> s.mass },
    Member("monoisotopicMass", null, "Quantity") { s, _ -> s.monoisotopicMass },
    Member("digest", null, "String") { s, _ -> s.digest },
    Member("find", 1..1, "List") { s, a -> s.find(a[0]) },
    Member("count", 1..1, "Int") { s, a -> s.count(a[0]) },
    Member("contains", 1..1, "Bool") { s, a -> s.contains(a[0]) },
    Member("slice", 2..2, "Sequence") { s, a -> s.slice(BioBag.intArg(a[0], "A slice's start"), BioBag.intArg(a[1], "A slice's end")) },
    Member("reverse", 0..0, "Sequence") { s, _ -> s.reverse() },
    Member("sameMolecule", 1..1, "Bool") { s, a -> s.sameMolecule(a[0]) },
)

private fun <T : NucleicAcid> nucleic(): List<Member<T>> = common<T>() + listOf(
    Member("ds", null, "Bool") { s, _ -> s.ds },
    Member("circular", null, "Bool") { s, _ -> s.circular },
    Member("gc", null, "Double") { s, _ -> s.gc },
    Member("complement", 0..0, "Sequence") { s, _ -> s.complement() },
    Member("reverseComplement", 0..0, "Sequence") { s, _ -> s.reverseComplement() },
    Member("codons", 0..1, "List") { s, a -> s.codons(BioBag.frameArg(a.getOrNull(0))) },
    Member("translate", 0..2, "Protein") { s, a -> s.translate(BioBag.frameArg(a.getOrNull(0)), a.getOrNull(1)) },
    Member("orfs", 0..2, "List") { s, a -> s.orfs(a.getOrNull(0)?.let { BioBag.intArg(it, "A minimum length") } ?: 100, a.getOrNull(1)) },
    Member("orf", 0..1, "Sequence") { s, a -> s.orf(a.getOrNull(0)) },
)

/** `Sequence`: any DNA, RNA or protein, for schemas (`type=Sequence`) and `is`. It has no literal of its own. */
public object SequenceType : ValueType<BioSequence> {
    override val name: String = "Sequence"
    override val kotlinClass: KClass<BioSequence> = BioSequence::class
    override val literals: List<String> = emptyList()
    override fun parse(literal: String): BioSequence = throw ParseException("Sequence has no literal; write .dna(…), .rna(…) or .protein(…)")
    override fun format(value: BioSequence): String = when (value) {
        is DNA -> DNAType.format(value)
        is RNA -> RNAType.format(value)
        is Protein -> ProteinType.format(value)
    }
    override fun encode(value: BioSequence): ByteArray = Packing.encode(value)
    override fun decode(bytes: ByteArray, version: Int): BioSequence = Packing.decode(bytes, bytes.firstOrNull()?.toInt()?.and(3) ?: -1)
    /**
     * Every member of DNA, RNA and protein, so a query over an attribute typed
     * `Sequence` checks against what any of them has (`c.insert.translate(1)`
     * where `c.insert.kind == "DNA"`). A value answers with its own type's members
     * when the query runs; this list is for the checker, and when used it hands
     * the call to the value's type.
     */
    override val members: List<Member<BioSequence>> by lazy {
        val seen = LinkedHashMap<String, Member<*>>()
        for (t in listOf(DNAType, RNAType, ProteinType)) for (m in t.members) seen.putIfAbsent(m.name, m)
        seen.values.map { m ->
            Member<BioSequence>(m.name, m.arity, m.resultType) { s, a ->
                val own = TypeBags.member(s, m.name) ?: throw BioException("A ${s.kind} has no ${m.name}")
                own.invoke(s, a)
            }
        }
    }
    override val samples: List<BioSequence> = emptyList()
}

/** `DNA`, `.dna("GAATTC")`, `.dna(ds circular "…")`. Orders by its letters. */
public object DNAType : ValueType<DNA> {
    override val name: String = "DNA"
    override val kotlinClass: KClass<DNA> = DNA::class
    override val literals: List<String> = listOf("dna")
    override fun parse(literal: String): DNA = readLiteral { LiteralParts(literal, true, "DNA").let { DNA(it.text, it.ds, it.circular) } }
    override fun format(value: DNA): String = nucleicLiteral("dna", value)
    override fun encode(value: DNA): ByteArray = Packing.encode(value)
    override fun decode(bytes: ByteArray, version: Int): DNA = Packing.decode(bytes, 0) as DNA
    override val ordered: Boolean = true
    override fun sortKey(value: DNA): ByteArray = Packing.sortKey(value)

    override val members: List<Member<DNA>> = nucleic<DNA>() + listOf(
        Member("transcribe", 0..0, "RNA") { s, _ -> s.transcribe() },
        Member("tm", 0..4, "Quantity") { s, a -> s.tm(a.getOrNull(0), a.getOrNull(1), a.getOrNull(2), a.getOrNull(3)) },
    )

    override val samples: List<DNA> by lazy {
        listOf(
            DNA("A"), DNA("ACGT"), DNA("GAATTC"), DNA("GAATTC", ds = true), DNA("ACGTRYSWKMBDHVN"),
            DNA("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"), DNA("TTTTTTTTTTTTTTTTTTTT"), DNA("NNNNNACGT"),
            DNA("GGATCCATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTC", ds = true, circular = true),
            DNA("ATGAAACCCGGGTTTTAA", circular = true), DNA(SAMPLE_313),
        )
    }
}

/** `RNA`, `.rna("AUGGCC")`. Orders by its letters. */
public object RNAType : ValueType<RNA> {
    override val name: String = "RNA"
    override val kotlinClass: KClass<RNA> = RNA::class
    override val literals: List<String> = listOf("rna")
    override fun parse(literal: String): RNA = readLiteral { LiteralParts(literal, true, "RNA").let { RNA(it.text, it.ds, it.circular) } }
    override fun format(value: RNA): String = nucleicLiteral("rna", value)
    override fun encode(value: RNA): ByteArray = Packing.encode(value)
    override fun decode(bytes: ByteArray, version: Int): RNA = Packing.decode(bytes, 1) as RNA
    override val ordered: Boolean = true
    override fun sortKey(value: RNA): ByteArray = Packing.sortKey(value)

    override val members: List<Member<RNA>> = nucleic<RNA>() + listOf(
        Member("reverseTranscribe", 0..0, "DNA") { s, _ -> s.reverseTranscribe() },
    )

    override val samples: List<RNA> by lazy {
        listOf(RNA("U"), RNA("AUG"), RNA("AUGGCCAUUGUAAUGGGCCGCUGA"), RNA("ACGURYSN"), RNA("GGGAAACCC", ds = true), RNA("AUGUUUUAA", circular = true))
    }
}

/** `Protein`, `.protein("MKTAYIAKQR")`. Orders by its letters. */
public object ProteinType : ValueType<Protein> {
    override val name: String = "Protein"
    override val kotlinClass: KClass<Protein> = Protein::class
    override val literals: List<String> = listOf("protein")
    override fun parse(literal: String): Protein = readLiteral { Protein(LiteralParts(literal, false, "protein").text) }
    override fun format(value: Protein): String = ".protein(\"${value.letters}\")"
    override fun encode(value: Protein): ByteArray = Packing.encode(value)
    override fun decode(bytes: ByteArray, version: Int): Protein = Packing.decode(bytes, 2) as Protein
    override val ordered: Boolean = true
    override fun sortKey(value: Protein): ByteArray = Packing.sortKey(value)

    override val members: List<Member<Protein>> = common<Protein>() + listOf(
        Member("pI", null, "Double") { s, _ -> s.pI },
        Member("charge", 1..1, "Double") { s, a -> s.charge((a[0] as? Number)?.toDouble() ?: throw BioException("A pH is a number, not ${a[0] ?: "nil"}")) },
        Member("extinction280", null, "Int") { s, _ -> s.extinction280 },
        Member("extinction280Reduced", null, "Int") { s, _ -> s.extinction280Reduced },
    )

    override val samples: List<Protein> by lazy {
        listOf(Protein("M"), Protein("MKTAYIAKQRQISFVKSHFSRQ"), Protein(Letters.PROTEIN), Protein("MKV*"), Protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"),
            Protein("PEPTIDE"), Protein("CCWWYY"))
    }
}

/** A fixed 313-nt sample (the start of lacZ alpha, as in pUC19), for the packing of a longer value. */
private const val SAMPLE_313 =
    "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGTACAATCTGCTCTG"
