package io.kixi.bio

import io.kixi.kd.KD
import io.kixi.revdb.RevDB
import io.kixi.types.TypeBags
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldBeEmpty
import io.kotest.matchers.doubles.plusOrMinus
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import io.kotest.matchers.string.shouldContain

/**
 * The synthetic backbones TNA, LNA and PNA (2026-10-29): each its own type rather
 * than a flag on DNA, so a melting temperature that does not hold is refused
 * rather than guessed, and a derived mass warns.
 */
class SyntheticBackbonesTest : FunSpec({

    test("the three are registered with their literals") {
        TypeBags.literal("tna") shouldBe TNAType
        TypeBags.literal("lna") shouldBe LNAType
        TypeBags.literal("pna") shouldBe PNAType
        TypeBags.typeFor(TNA("ACGT")) shouldBe TNAType
        TypeBags.typeFor(PNA("ACGT")) shouldBe PNAType
        BioBag.VERSION shouldBe 2
    }

    test("KD reads the letters bare or quoted, and writes them bare") {
        val t = KD.read("""oligo tna=.tna(ATGCGTACGTTA) lna=.lna("gcatgc") pna=.pna(CATTGC) ds=.tna(ds circular ACGT)""")
        t["tna"] shouldBe TNA("ATGCGTACGTTA")
        t["lna"] shouldBe LNA("GCATGC")
        t["pna"] shouldBe PNA("CATTGC")
        t["ds"] shouldBe TNA("ACGT", ds = true, circular = true)
        t.toString() shouldContain "tna=.tna(ATGCGTACGTTA)"
        t.toString() shouldContain "pna=.pna(CATTGC)"
        t.toString() shouldContain "ds=.tna(ds circular ACGT)"
    }

    test("the letters are checked as DNA's, and a PNA takes no strand flags") {
        shouldThrow<Exception> { KD.read("t s=.tna(AUGC)") }.message!! shouldContain "'U'"
        shouldThrow<Exception> { KD.read("t s=.lna(ATGX)") }.message!! shouldContain "'X'"
        shouldThrow<Exception> { KD.read("t s=.pna(ATGU)") }.message!! shouldContain "'U'"
        // An initialism takes its article from how it is said: an RNA, an LNA, a TNA.
        shouldThrow<Exception> { KD.read("t s=.lna(AU)") }.message!! shouldContain "Not an LNA sequence"
        shouldThrow<Exception> { KD.read("t s=.tna(AU)") }.message!! shouldContain "Not a TNA sequence"
        shouldThrow<Exception> { KD.read("t s=.pna(ds ACGT)") }.message!! shouldContain "A PNA has no 'ds'"
        // IUPAC codes are DNA's.
        TNA("ACGTRYSWKMBDHVN").length shouldBe 15
        PNA("ACGTRYSWKMBDHVN").length shouldBe 15
    }

    test("a TNA is not a DNA with the same letters") {
        TNA("ACGT") shouldNotBe DNA("ACGT")
        DNA("ACGT") shouldNotBe TNA("ACGT")
        TNA("ACGT").sameMolecule(DNA("ACGT")) shouldBe false
        LNA("ACGT") shouldNotBe TNA("ACGT")
        PNA("ACGT").sameMolecule(PNA("ACGT")) shouldBe true
        TNA("ACGT").kind shouldBe "TNA"
        LNA("ACGT").kind shouldBe "LNA"
        PNA("ACGT").kind shouldBe "PNA"
    }

    test("letter operations carry over; a melting temperature does not") {
        TNA("GAATTC").reverseComplement() shouldBe TNA("GAATTC")
        LNA("ATGC").complement() shouldBe LNA("TACG")
        LNA("ATGC").reverseComplement() shouldBe LNA("GCAT")
        TNA("GGCC").gc shouldBe 100.0
        PNA("ATGC").gc shouldBe 50.0
        PNA("ATGC").reverseComplement() shouldBe PNA("GCAT")
        PNA("GAATTCNN").find("GAATTC") shouldBe listOf(1)
        TNA("ATGCGT").slice(2, 4) shouldBe TNA("TGC")
        // No tm on any of them: DNA alone has one, since no model here fits these backbones.
        TypeBags.member(TNA("ACGTACGT"), "tm") shouldBe null
        TypeBags.member(LNA("ACGTACGT"), "tm") shouldBe null
        TypeBags.member(PNA("ACGTACGT"), "tm") shouldBe null
        TypeBags.member(DNA("ACGTACGT"), "tm") shouldNotBe null
        // Nor transcription, which is DNA's, or reverse transcription, which is RNA's.
        TypeBags.member(TNA("ACGT"), "transcribe") shouldBe null
        TypeBags.member(PNA("ACGT"), "circular") shouldBe null
        TypeBags.member(PNA("ACGT"), "ds") shouldBe null
    }

    test("derived masses, and the warning that says so once") {
        // TNA is each DNA residue less CH2 (threose for deoxyribose).
        val dna = DNA("ACGT").mass!!.value.toDouble()
        val tna = TNA("ACGT").mass!!.value.toDouble()
        (dna - tna) shouldBe (4 * 14.0266 plusOrMinus 0.01)
        // LNA is each RNA residue plus the bridge's carbon, and uses T rather than U.
        val lna = LNA("ACG").mass!!.value.toDouble()
        val rna = RNA("ACG").mass!!.value.toDouble()
        (lna - rna) shouldBe (3 * 12.011 plusOrMinus 0.01)
        // The PNA residues reproduce the values PNA mass calculators use.
        PNA("A").mass!!.value.toDouble() shouldBe (275.27 + 18.015 plusOrMinus 0.02)
        PNA("ACGT").mass!!.value.toDouble() shouldBe (275.27 + 251.25 + 291.27 + 266.26 + 18.015 plusOrMinus 0.05)
        // An ambiguous letter has no formula, as for DNA.
        TNA("ACGN").formula shouldBe null
        PNA("ACGN").mass shouldBe null

        val lines = ArrayList<String>()
        val previous = io.kixi.logHandler
        io.kixi.logHandler = { level, message -> if (level == io.kixi.LogLevel.WARN) lines += message }
        try {
            // Already warned above, once per kind per JVM, so these add nothing.
            TNA("ACGTACGT").mass
            PNA("ACGTACGT").mass
            Masses.warnUnverified("TNA")
            lines.shouldBeEmpty()
            Masses.warnUnverified("XNA-for-this-test")
            lines.size shouldBe 1
            lines[0] shouldContain "are derived from the backbone's chemistry and are not verified"
            lines[0] shouldContain "provisional"
        } finally {
            io.kixi.logHandler = previous
        }
        // The article follows how an initialism is said, so the warning reads
        // "an LNA mass" and "a TNA mass". Two kinds not warned about yet:
        val said = ArrayList<String>()
        io.kixi.logHandler = { _, message -> said += message }
        try {
            Masses.warnUnverified("ZNA")
            said.last() shouldContain "Treat a ZNA mass"
            Masses.warnUnverified("FNA")
            said.last() shouldContain "Treat an FNA mass"
        } finally {
            io.kixi.logHandler = previous
        }
        Letters.article("LNA") shouldBe "an"
        Letters.article("TNA") shouldBe "a"
        Letters.article("RNA") shouldBe "an"
        Letters.article("DNA") shouldBe "a"
        Letters.article("protein") shouldBe "a"
    }

    test("the stored form is a byte wider in its kind alone, and version 1's bytes are unchanged") {
        // The head byte still holds everything: one byte, whatever the kind.
        Packing.encode(TNA("ACGT")).size shouldBe 1 + 1 + 1
        Packing.encode(DNA("ACGT")).size shouldBe 1 + 1 + 1
        Packing.encode(PNA("ACGTRY")).size shouldBe 1 + 1 + 3
        Packing.encode(TNA("ACGT".repeat(2500), ds = true, circular = true)).size shouldBe 1 + 2 + 2500
        // Version 1 wrote the kind in bits 0-1 with bit 6 zero, and still reads back.
        Packing.encode(DNA("ACGT"))[0].toInt() and 0x40 shouldBe 0
        Packing.encode(RNA("ACGU"))[0].toInt() and 0x40 shouldBe 0
        Packing.encode(Protein("MKT"))[0].toInt() and 0x40 shouldBe 0
        // TNA fits the value version 1 left spare; PNA and LNA use bit 6.
        Packing.kindOf(TNA("ACGT")) shouldBe 3
        Packing.kindOf(PNA("ACGT")) shouldBe 4
        Packing.kindOf(LNA("ACGT")) shouldBe 5
        Packing.encode(PNA("ACGT"))[0].toInt() and 0x40 shouldNotBe 0

        for (s in listOf(
            TNA("A"), TNA("ACGTRYSWKMBDHVN"), TNA("GAATTC", ds = true), TNA("ACGT", ds = true, circular = true),
            LNA("GCATGC"), LNA("CTGATC", ds = true), PNA("A"), PNA("CATTGC"), PNA("ACGTRY"),
        )) {
            val bytes = Packing.encode(s)
            Packing.decode(bytes, Packing.kindIn(bytes)) shouldBe s
            SequenceType.decode(bytes, BioBag.VERSION) shouldBe s
        }
        // A kind's own decode refuses another kind's bytes.
        shouldThrow<IllegalArgumentException> { Packing.decode(Packing.encode(TNA("ACGT")), 0) }
        shouldThrow<IllegalArgumentException> { Packing.decode(Packing.encode(DNA("ACGT")), 3) }
        // A stored PNA with a strand flag is damage: a PNA is one strand.
        val flagged = Packing.encode(PNA("ACGT")).also { it[0] = (it[0].toInt() or 0x04).toByte() }
        shouldThrow<IllegalArgumentException> { Packing.decode(flagged, 4) }
    }

    test("KS: literals, members and the bag functions") {
        ks(
            """
            use ki.bio.*
            say .tna(ATGCGTACGTTA).length
            say .lna(ATGC).reverseComplement()
            say .pna(CATTGC) is PNA
            say .tna(ACGT) == .dna(ACGT)
            say tna("acgt atgc")
            say pna("CATTGC").gc
            say revcomp(.pna(ATGC))
            say seq("ACGT", "tna").kind
            """.trimIndent()
        ) shouldBe "12\nGCAT\ntrue\nfalse\nACGTATGC\n50\nGCAT\nTNA"
        // A sequence of another kind is not quietly re-backboned.
        shouldThrow<Exception> { ks("use ki.bio.*\nsay tna(.dna(ACGT))") }.message!! shouldContain "Expected TNA"
        shouldThrow<Exception> { ks("use ki.bio.*\nsay revcomp(.protein(\"MKT\"))") }.message!! shouldContain "not a protein"
    }

    test("KQ and RevDB: a TNA index, and a Sequence attribute that keeps the kinds apart") {
        RevDB.inMemory().use { db ->
            db.writeBlocking {
                createCollection(
                    """
                    rev:collection oligos {
                        schema additionalAttributes=true { attribute seq type=Sequence required=true }
                        index seq
                    }
                    """.trimIndent()
                )
                createCollection(
                    """
                    rev:collection aptamers {
                        schema additionalAttributes=true { attribute seq type=TNA required=true }
                        index seq unique=true
                    }
                    """.trimIndent()
                )
                // The same letters in four backbones: four different molecules.
                insert("oligos", "oligo name=\"d\" seq=.dna(ATGCGT)")
                insert("oligos", "oligo name=\"t\" seq=.tna(ATGCGT)")
                insert("oligos", "oligo name=\"l\" seq=.lna(ATGCGT)")
                insert("oligos", "oligo name=\"p\" seq=.pna(ATGCGT)")
                insert("aptamers", "aptamer name=\"A1\" seq=.tna(GGTTGGTGTGGTTGG)")
                shouldThrow<Exception> { insert("aptamers", "aptamer name=\"dup\" seq=.tna(ggttggtgtggttgg)") }
                    .message!! shouldContain "Unique"
                shouldThrow<Exception> { insert("aptamers", "aptamer name=\"wrong\" seq=.dna(ACGT)") }
                    .message!! shouldContain "TNA"
            }
            // An index lookup finds the one backbone asked for, not the letters.
            db.kq("from o in oligos where o.seq == .tna(ATGCGT) select o.name").map { it["name"] } shouldBe listOf("t")
            db.kq("from o in oligos where o.seq == .pna(ATGCGT) select o.name").map { it["name"] } shouldBe listOf("p")
            db.readBlocking {
                kq("explain from o in oligos where o.seq == .tna(ATGCGT) select o.name").toList().toString()
            } shouldContain "lookup"
            db.kq("from o in oligos orderby o.kind, o.name select o.seq.kind as k").map { it["k"] }.toSet() shouldBe
                setOf("DNA", "TNA", "LNA", "PNA")
            db.kq("from a in aptamers where a.seq.gc > 50 select a.name").map { it["name"] } shouldBe listOf("A1")
            db.kq("from o in oligos where o.seq.kind == \"TNA\" select o.seq.reverseComplement() as rc")
                .map { it["rc"] } shouldBe listOf(TNA("ACGCAT"))
            db.kq("from a in aptamers select a.seq.length as n").map { it["n"] } shouldBe listOf(15)
        }
    }
})
