package io.kixi.bio

/**
 * RevDB's stored form of a sequence (bag version 2), and its sort key.
 *
 * Stored: a head byte (the kind in bits 0-1 with bit 6 as its high bit: 0 DNA,
 * 1 RNA, 2 protein, 3 TNA, 4 PNA, 5 LNA; bit 2 double-stranded; bit 3 circular;
 * bits 4-5 the packing: 0 two bits a base, for plain A C G T/U; 1 four bits, the
 * IUPAC base set; 2 five bits a residue, protein), the length as an unsigned
 * LEB128 varint, then the letters packed from the high bits down, the last byte's
 * spare bits zero. A 10 kb plasmid takes 2.5 KB.
 *
 * Version 2 widened the kind into bit 6, which version 1 reserved and wrote as
 * zero, so every value version 1 stored encodes to the same bytes here and decodes
 * unchanged: nothing needs rewriting. Bit 7 stays reserved.
 *
 * Sort key: the letters' ASCII bytes, a 0 byte, then the flags byte, so a type
 * orders by its letters (A < C < G < T, lexical) and two keys are equal exactly
 * when the values are. The kind is not in the key and need not be: RevDB writes
 * the bag and type names ahead of it, taken from the value's own type, so a TNA
 * and a DNA with the same letters are different keys even in one attribute
 * declared `type=Sequence`.
 */
internal object Packing {

    private const val PLAIN = 0
    private const val IUPAC = 1
    private const val RESIDUES = 2

    /** Longest sequence decoded: 2^31 - 1 letters is the String limit anyway. */
    private const val MAX_LENGTH = Int.MAX_VALUE - 8

    /** The stored kind of a value: 0 DNA, 1 RNA, 2 protein, 3 TNA, 4 PNA, 5 LNA. */
    fun kindOf(s: BioSequence): Int =
        when (s) { is DNA -> 0; is RNA -> 1; is Protein -> 2; is TNA -> 3; is PNA -> 4; is LNA -> 5 }

    /** The kind [bytes] holds, or -1 when there is no head byte; -1 never matches an expected kind. */
    fun kindIn(bytes: ByteArray): Int {
        val head = bytes.firstOrNull()?.toInt()?.and(0xFF) ?: return -1
        return (head and 3) or (((head ushr 6) and 1) shl 2)
    }

    fun encode(s: BioSequence): ByteArray {
        val kind = kindOf(s)
        val letters = s.letters
        val packing = when {
            s is Protein -> RESIDUES
            letters.all { Letters.isPlainBase(it) } -> PLAIN
            else -> IUPAC
        }
        val bits = when (packing) { PLAIN -> 2; IUPAC -> 4; else -> 5 }
        val head = (kind and 3) or ((s.flags and 3) shl 2) or (packing shl 4) or (((kind ushr 2) and 1) shl 6)
        val out = java.io.ByteArrayOutputStream(2 + 5 + (letters.length.toLong() * bits / 8 + 1).toInt())
        out.write(head)
        var n = letters.length
        while (n >= 0x80) { out.write((n and 0x7F) or 0x80); n = n ushr 7 }
        out.write(n)
        var acc = 0L
        var have = 0
        for (ch in letters) {
            val code = when (packing) {
                PLAIN -> when (ch) { 'A' -> 0; 'C' -> 1; 'G' -> 2; else -> 3 }
                IUPAC -> Letters.mask(ch)
                else -> Letters.PROTEIN.indexOf(ch)
            }
            acc = (acc shl bits) or code.toLong()
            have += bits
            while (have >= 8) { have -= 8; out.write(((acc ushr have) and 0xFF).toInt()) }
        }
        if (have > 0) out.write(((acc shl (8 - have)) and 0xFF).toInt())
        return out.toByteArray()
    }

    /** Decodes [bytes] as a value of [expectedKind] ([kindOf]); any damage is an IllegalArgumentException. */
    fun decode(bytes: ByteArray, expectedKind: Int): BioSequence {
        require(bytes.isNotEmpty()) { "A stored sequence is not empty" }
        val head = bytes[0].toInt() and 0xFF
        val kind = (head and 3) or (((head ushr 6) and 1) shl 2)
        val flags = (head ushr 2) and 3
        val packing = (head ushr 4) and 3
        require(head ushr 7 == 0) { "Unknown bits in a stored sequence's head" }
        require(kind <= 5) { "A stored sequence of unknown kind $kind" }
        require(kind == expectedKind) { "A stored sequence of kind $kind where kind $expectedKind was expected" }
        require(
            when (kind) {
                2 -> packing == RESIDUES && flags == 0
                4 -> (packing == PLAIN || packing == IUPAC) && flags == 0   // a PNA is one strand, never circular
                else -> packing == PLAIN || packing == IUPAC
            }
        ) { "A stored sequence's packing does not suit its kind" }
        var pos = 1
        var len = 0L
        var shift = 0
        while (true) {
            require(pos < bytes.size && shift <= 28) { "A stored sequence's length is cut short" }
            val b = bytes[pos++].toInt() and 0xFF
            len = len or ((b and 0x7F).toLong() shl shift)
            if (b and 0x80 == 0) break
            shift += 7
        }
        require(len in 1..MAX_LENGTH.toLong()) { "A stored sequence of length $len" }
        val n = len.toInt()
        val bits = when (packing) { PLAIN -> 2; IUPAC -> 4; else -> 5 }
        val need = ((n.toLong() * bits + 7) / 8)
        require(bytes.size - pos.toLong() == need) { "A stored sequence of $n letters needs $need bytes, not ${bytes.size - pos}" }
        val sb = StringBuilder(n)
        var acc = 0L
        var have = 0
        val mask = (1 shl bits) - 1
        for (i in 0 until n) {
            while (have < bits) { acc = (acc shl 8) or (bytes[pos++].toLong() and 0xFF); have += 8 }
            have -= bits
            val code = ((acc ushr have).toInt()) and mask
            sb.append(when (packing) {
                PLAIN -> if (kind == 1) "ACGU"[code] else "ACGT"[code]
                IUPAC -> {
                    require(code != 0) { "A stored base set is empty" }
                    if (kind == 1) Letters.rnaLetter(code) else Letters.dnaLetter(code)
                }
                else -> {
                    require(code < Letters.PROTEIN.length) { "A stored residue code $code" }
                    Letters.PROTEIN[code]
                }
            })
            acc = acc and ((1L shl have) - 1)
        }
        require(acc == 0L) { "A stored sequence's spare bits are not zero" }
        val letters = sb.toString()
        val ds = flags and 1 != 0
        val circular = flags and 2 != 0
        return when (kind) {
            0 -> DNA(letters, ds, circular)
            1 -> RNA(letters, ds, circular)
            2 -> Protein(letters)
            3 -> TNA(letters, ds, circular)
            4 -> PNA(letters)
            else -> LNA(letters, ds, circular)
        }
    }

    fun sortKey(s: BioSequence): ByteArray {
        val l = s.letters
        val out = ByteArray(l.length + 2)
        for (i in l.indices) out[i] = l[i].code.toByte()
        out[l.length] = 0
        out[l.length + 1] = s.flags.toByte()
        return out
    }
}
