package io.kixi.bio

/**
 * RevDB's stored form of a sequence (bag version 1), and its sort key.
 *
 * Stored: a head byte (bits 0-1 the kind: 0 DNA, 1 RNA, 2 protein; bit 2
 * double-stranded; bit 3 circular; bits 4-5 the packing: 0 two bits a base,
 * for plain A C G T/U; 1 four bits, the IUPAC base set; 2 five bits a residue,
 * protein), the length as an unsigned LEB128 varint, then the letters packed
 * from the high bits down, the last byte's spare bits zero. A 10 kb plasmid
 * takes 2.5 KB.
 *
 * Sort key: the letters' ASCII bytes, a 0 byte, then the flags byte, so a type
 * orders by its letters (A < C < G < T, lexical) and two keys are equal exactly
 * when the values are.
 */
internal object Packing {

    private const val PLAIN = 0
    private const val IUPAC = 1
    private const val RESIDUES = 2

    /** Longest sequence decoded: 2^31 - 1 letters is the String limit anyway. */
    private const val MAX_LENGTH = Int.MAX_VALUE - 8

    fun encode(s: BioSequence): ByteArray {
        val kind = when (s) { is DNA -> 0; is RNA -> 1; is Protein -> 2 }
        val letters = s.letters
        val packing = when {
            s is Protein -> RESIDUES
            letters.all { Letters.isPlainBase(it) } -> PLAIN
            else -> IUPAC
        }
        val bits = when (packing) { PLAIN -> 2; IUPAC -> 4; else -> 5 }
        val head = kind or ((s.flags and 3) shl 2) or (packing shl 4)
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

    /** Decodes [bytes] as a value of [expectedKind] (0 DNA, 1 RNA, 2 protein); any damage is an IllegalArgumentException. */
    fun decode(bytes: ByteArray, expectedKind: Int): BioSequence {
        require(bytes.isNotEmpty()) { "A stored sequence is not empty" }
        val head = bytes[0].toInt() and 0xFF
        val kind = head and 3
        val flags = (head ushr 2) and 3
        val packing = (head ushr 4) and 3
        require(head ushr 6 == 0) { "Unknown bits in a stored sequence's head" }
        require(kind == expectedKind) { "A stored sequence of kind $kind where kind $expectedKind was expected" }
        require(if (kind == 2) packing == RESIDUES && flags == 0 else packing == PLAIN || packing == IUPAC) { "A stored sequence's packing does not suit its kind" }
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
            else -> Protein(letters)
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
