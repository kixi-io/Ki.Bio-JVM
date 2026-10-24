package io.kixi.bio

import io.kixi.kd.KD
import io.kixi.revdb.RevDB
import io.kixi.sci.Formula
import io.kixi.types.TypeBags
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldBeEmpty
import io.kotest.matchers.doubles.plusOrMinus
import io.kotest.matchers.shouldBe
import io.kotest.matchers.string.shouldContain
import io.kotest.matchers.types.shouldBeInstanceOf

/** The `bio` bag as KD, KS, KQ and RevDB see it, and the sequence operations. */
class BioBagTest : FunSpec({

    val bag = TypeBags.byName("bio")!!

    test("the bag is discovered from the classpath, with no problems") {
        bag.shouldBeInstanceOf<BioBag>()
        TypeBags.problems.filter { "bio" in it }.shouldBeEmpty()
        TypeBags.literal("dna") shouldBe DNAType
        TypeBags.literal("RNA") shouldBe RNAType
        TypeBags.literal("protein") shouldBe ProteinType
        TypeBags.typeNamed("Sequence") shouldBe SequenceType
        TypeBags.typeFor(DNA("ACGT")) shouldBe DNAType
    }

    test("every sample passes the contract's self-check and the corpus round trip") {
        corpusProblems(bag).shouldBeEmpty()
    }

    test("KD reads the literals, with ds and circular before a quoted string") {
        val t = KD.read("""construct name="pUC-insert" insert=.dna(ds circular "ggatccATGAAAtaa") primer=.dna(GTAAAACGACGGCCAGT) mrna=.rna("AUGGCC") p=.protein("mktay")""")
        t["insert"] shouldBe DNA("GGATCCATGAAATAA", ds = true, circular = true)
        t["primer"] shouldBe DNA("GTAAAACGACGGCCAGT")
        t["mrna"] shouldBe RNA("AUGGCC")
        t["p"] shouldBe Protein("MKTAY")
        t.toString() shouldContain "insert=.dna(ds circular \"GGATCCATGAAATAA\")"
        KD.read("t s=.dna(ds AC)")["s"] shouldBe DNA("AC", ds = true)
        shouldThrow<Exception> { KD.read("t s=.dna(\"ACGU\")") }.message!! shouldContain "'U'"
        shouldThrow<Exception> { KD.read("t s=.protein(ds \"MKT\")") }.message!! shouldContain "no 'ds'"
        shouldThrow<Exception> { KD.read("t s=.dna(ds ds \"A\")") }.message!! shouldContain "twice"
        shouldThrow<Exception> { KD.read("t s=.dna(\"\")") }
        // The literal is letters only (design 5.5): no spaces, and a flag is not a sequence.
        shouldThrow<Exception> { KD.read("t s=.dna(\"ds AC\")") }.message!! shouldContain "without spaces"
        shouldThrow<Exception> { KD.read("t s=.dna(ds)") }.message!! shouldContain "needs its letters"
        KD.read("t s=.dna(\"DS\")")["s"] shouldBe DNA("DS")
        DNA("ggatcc ATG") shouldBe DNA("GGATCCATG")                      // the constructor and dna(text) read looser text
    }

    test("KS: literals, members, functions and `in`") {
        val out = ks(
            """
            use ki.bio.*
            let p = .dna("GTAAAACGACGGCCAGT")
            say p.length
            say p.reverseComplement()
            say p is DNA
            say p is Sequence
            say "GACGG" in p
            say "GAYGG" in p
            say p.find("AAA")
            say revcomp("ATGC")
            say gc("GGCC")
            say translate("ATGGCCTAA")
            say .rna("AUGGCC").reverseTranscribe()
            say .protein("MKT").kind
            say seq(">x\n1 atgc").kind
            say a260(0.5, "dsdna")
            """.trimIndent()
        )
        out.lines() shouldBe listOf("17", "ACTGGCCGTCGTTTTAC", "true", "true", "true", "true", "[3, 4]", "GCAT", "100", "MA*", "ATGGCC", "Protein", "DNA", "25")
    }

    test("KQ and RevDB: a typed collection, an index lookup by sequence, members in queries") {
        RevDB.inMemory().use { db ->
            db.writeBlocking {
                createCollection("""
                    rev:collection primers {
                        schema additionalAttributes=true { attribute seq type=DNA required=true }
                        index seq unique=true
                    }
                """.trimIndent())
                createCollection("rev:collection constructs { schema additionalAttributes=true { attribute insert type=Sequence } }")
                insert("primers", "primer name=\"M13F\" seq=.dna(GTAAAACGACGGCCAGT)")
                insert("primers", "primer name=\"T7\" seq=.dna(TAATACGACTCACTATAGGG)")
                insert("primers", "primer name=\"M13R\" seq=.dna(CAGGAAACAGCTATGAC)")
                insert("constructs", "construct name=\"a\" insert=.dna(ds circular \"GGATCCATGAAATAA\")")
                insert("constructs", "construct name=\"b\" insert=.protein(\"MKT\")")
                shouldThrow<Exception> { insert("primers", "primer name=\"dup\" seq=.dna(gtaaaacgacggccagt)") }.message!! shouldContain "Unique"
                shouldThrow<Exception> { insert("primers", "primer name=\"text\" seq=\"ACGT\"") }.message!! shouldContain "DNA"
            }
            db.kq("from p in primers where p.seq == .dna(TAATACGACTCACTATAGGG) select p.name").map { it["name"] } shouldBe listOf("T7")
            db.readBlocking { kq("explain from p in primers where p.seq == .dna(TAATACGACTCACTATAGGG) select p.name").toList().toString() } shouldContain "lookup"
            db.kq("from p in primers orderby p.seq select p.name").map { it["name"] } shouldBe listOf("M13R", "M13F", "T7")
            db.kq("from p in primers where p.seq.gc > 45 orderby p.seq select p.name, p.seq.length as n").map { it["name"] } shouldBe listOf("M13R", "M13F")
            db.kq("from p in primers where \"ACGAC\" in p.seq select p.name").map { it["name"] }.toSet() shouldBe setOf("M13F", "T7")
            db.kq("from c in constructs where c.insert.kind == \"DNA\" select c.insert.translate(1).text as p").map { it["p"] } shouldBe listOf("GSMK*")
            db.kq("from p in primers where p.name == \"M13F\" select revcomp(p.seq) as rc").map { it["rc"] } shouldBe listOf(DNA("ACTGGCCGTCGTTTTAC"))
            // A Sequence-typed attribute compares with a DNA literal (RevDB's checker relates a type to its subtypes).
            db.kq("from c in constructs where c.insert == .dna(ds circular GGATCCATGAAATAA) select c.name").map { it["name"] } shouldBe listOf("a")
            db.kq("from c in constructs where c.insert != .dna(ACGT) orderby c.name select c.name").map { it["name"] } shouldBe listOf("a", "b")
        }
    }

    test("gc: G, C and S over the bases known to be G/C or A/T, as Biopython's gc_fraction") {
        DNA("GCAT").gc shouldBe 50.0
        DNA("GCNN").gc shouldBe 100.0
        DNA("SWNN").gc shouldBe 50.0
        RNA("GGCU").gc shouldBe 75.0
        DNA("NNNN").gc shouldBe null
        ks("use ki.bio.*\nsay gc(.dna(\"NNRY\"))") shouldBe "nil"
    }

    test("strands: complement with IUPAC codes, transcription, reverse transcription") {
        DNA("ACGTRYSWKMBDHVN").complement().letters shouldBe "TGCAYRSWMKVHDBN"
        DNA("AACG").reverseComplement().letters shouldBe "CGTT"
        RNA("AUGC").complement().letters shouldBe "UACG"
        DNA("ATGC", ds = true, circular = true).transcribe() shouldBe RNA("AUGC", ds = true, circular = true)
        RNA("AUGC").reverseTranscribe() shouldBe DNA("ATGC")
        DNA("ACGT").reverse().letters shouldBe "TGCA"
    }

    test("find, count and contains: IUPAC patterns, overlaps, across the origin of a circular sequence") {
        val s = DNA("GAATTCAAAGGATCC")
        s.find("GAATTC") shouldBe listOf(1)
        s.find("GGATYC") shouldBe listOf(10)
        s.find("AA") shouldBe listOf(2, 7, 8)
        s.count("N") shouldBe 15
        DNA("NNN").find("A") shouldBe emptyList()
        DNA("ACGN").find("ACGN") shouldBe listOf(1)
        val plasmid = DNA("TTCAAAAAAAGAA", circular = true)
        plasmid.find("GAATTC") shouldBe listOf(11)      // G A A at 11 to 13, T T C at 1 to 3
        DNA("TTCAAAAAAAGAA").find("GAATTC") shouldBe emptyList()
        plasmid.contains("GAATTC") shouldBe true
        RNA("AUGGCC").find("ATG") shouldBe listOf(1)
        Protein("MKTAYIAKQR").find("AXI") shouldBe listOf(4)
        Protein("MKDNE").find("B") shouldBe listOf(3, 4)
        shouldThrow<BioException> { DNA("ACGT").find("MKX*") }
        shouldThrow<BioException> { DNA("ACGT").find(Protein("MK")) }
    }

    test("slice: 1-based and inclusive, across the origin when circular") {
        DNA("ACGTACGT").slice(2, 4).letters shouldBe "CGT"
        DNA("ACGTACGT", ds = true, circular = true).slice(7, 2) shouldBe DNA("GTAC", ds = true)
        shouldThrow<BioException> { DNA("ACGT").slice(3, 2) }
        shouldThrow<BioException> { DNA("ACGT").slice(0, 2) }
        Protein("MKTAY").slice(2, 3) shouldBe Protein("KT")
    }

    test("codons and translation: frames on both strands, ambiguous codons, other codes") {
        val s = DNA("ATGGCCATTGTAATGGGCCGCTGA")
        s.codons() shouldBe listOf("ATG", "GCC", "ATT", "GTA", "ATG", "GGC", "CGC", "TGA")
        RNA("AUGGCCUAA").codons() shouldBe listOf("AUG", "GCC", "UAA")
        RNA("AUGGCCUAA").codons(-1) shouldBe listOf("UUA", "GGC", "CAU")
        s.translate().letters shouldBe "MAIVMGR*"
        s.translate(-1).letters shouldBe "SAAHYNGH"
        DNA("GGNTAR").translate().letters shouldBe "G*"
        DNA("NNN").translate().letters shouldBe "X"
        DNA("RAYSARMTHYAA").translate().letters shouldBe "BZJX"           // D or N, E or Q, I or L, Q or a stop (as Biopython)
        DNA("TGA").translate(1, 2).letters shouldBe "W"
        RNA("AUGUGA").translate().letters shouldBe "M*"
        shouldThrow<BioException> { s.translate(4) }.message!! shouldContain "frame"
        shouldThrow<BioException> { s.translate(1, 7) }.message!! shouldContain "no NCBI genetic code 7"
        shouldThrow<BioException> { DNA("AT").translate() }
    }

    test("ORFs: both strands, the longest for each stop, across the origin of a plasmid") {
        val linear = DNA("CCATGAAATTTTAGCCCTTACATTTCATGG")  // ATGAAATTTTAG on top; ATGTAA on the other strand
        // Top strand: ATG AAA TTT TAG. The other strand, CCATGAAATGTAAGGGCTAAAATTTCATGG: ATG AAA TGT AAG GGC TAA, and ATG TAA in another frame.
        linear.orfs(3).map { it.letters } shouldBe listOf("ATGAAATGTAAGGGCTAA", "ATGAAATTTTAG", "ATGTAA")
        linear.orfs(12).size shouldBe 2
        linear.orf().letters shouldBe "ATGAAATGTAAGGGCTAA"
        DNA("ATGATGAAATAG").orfs(3).map { it.letters } shouldBe listOf("ATGATGAAATAG")
        // An ORF that runs across the origin: TAG CCC | ATG AAA (origin) GGG ... read from the ATG round to the TAG.
        val plasmid = DNA("GGGTAGCCCATGAAA", circular = true)
        plasmid.orfs(3).first().letters shouldBe "ATGAAAGGGTAG"
        DNA("GGGTAGCCCATGAAA").orfs(3).none { it.letters == "ATGAAAGGGTAG" } shouldBe true
        shouldThrow<BioException> { DNA("CCCCCC").orf() }
        // On a circle each ATG is read on its own: an ORF from an earlier ATG that would be longer
        // than the circle does not hide a shorter one from a later ATG to the same stop.
        DNA("TAAATGATGATGAT", circular = true).orfs(3).map { it.letters } shouldBe listOf("ATGATTAAATGA")
        DNA("CCATGAAATGATGATGG", circular = true).orfs(3, 27).first().letters shouldBe "ATGGCCATGAAATGA"
        RNA("GGGUAGCCCAUGAAA", circular = true).orfs(3).first().letters shouldBe "AUGAAAGGGUAG"
    }

    test("masses from Ki.Sci formulas: oligos with a 5' hydroxyl, duplexes, circles, proteins") {
        // IDT's formula: dA 313.209, dC 289.184, dG 329.208, dT 304.196, less 61.964 (PO2H less H2).
        fun idt(s: String) = s.sumOf { when (it) { 'A' -> 313.209; 'C' -> 289.184; 'G' -> 329.208; else -> 304.196 } } - 61.964
        for (s in listOf("GTAAAACGACGGCCAGT", "TAATACGACTCACTATAGGG", "A", "ACGTACGTACGT")) {
            DNA(s).mass!!.value.toDouble() shouldBe (idt(s) plusOrMinus 0.002 * s.length)
        }
        DNA("A").formula shouldBe Formula.parse("C10H13N5O3")            // deoxyadenosine
        DNA("ACGT", circular = true).formula!!.count(io.kixi.sci.Elements.bySymbol("P")!!) shouldBe 4
        DNA("ACGT").formula!!.count(io.kixi.sci.Elements.bySymbol("P")!!) shouldBe 3
        DNA("AACG", ds = true).formula shouldBe DNA("AACG").formula!! + DNA("CGTT").formula!!
        Protein("G").formula shouldBe Formula.parse("C2H5NO2")            // glycine
        Protein("GG*").formula shouldBe Formula.parse("C4H8N2O3")         // a final stop is ignored
        Protein("MK*V").formula shouldBe null
        Protein("*").mass shouldBe null                                   // translate("TAA"): no residues
        DNA("ACGN").mass shouldBe null
        Protein("MKX").mass shouldBe null
        Protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG").mass!!.value.toDouble() shouldBe (8564.84 plusOrMinus 0.05)   // ubiquitin, ExPASy 8564.84
    }

    test("protein numbers: pI as ExPASy, extinction at 280 nm as Pace (ProtParam's two values)") {
        Protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG").pI!! shouldBe (6.56 plusOrMinus 0.005)
        Protein("CCWWYY").extinction280 shouldBe 2 * 5500 + 2 * 1490 + 125
        Protein("CCWWYY").extinction280Reduced shouldBe 2 * 5500 + 2 * 1490
        Protein("MK*V").pI shouldBe null
    }

    test("Tm: IDT's defaults, conditions as quantities, refusals") {
        val t = DNA("GTAAAACGACGGCCAGT").tm().value.toDouble()
        // The order is (Na+, oligo, Mg2+, dNTPs), RevGrid's tm(dna, salt, oligo) with the buffer after.
        val salty = ks("say .dna(\"GTAAAACGACGGCCAGT\").tm(50mM, 250nM, 2mM, 0.2mM) > .dna(\"GTAAAACGACGGCCAGT\").tm()")
        salty shouldBe "true"
        t shouldBe (DNA("GTAAAACGACGGCCAGT").tm(0.05, 250e-9, 0, 0).value.toDouble() plusOrMinus 1e-9)
        ks("use ki.bio.*\nsay tm(\"GTAAAACGACGGCCAGT\", 50mM, 250nM, 2mM, 0.2mM) == .dna(\"GTAAAACGACGGCCAGT\").tm(50mM, 250nM, 2mM, 0.2mM)") shouldBe "true"
        shouldThrow<BioException> { DNA("ACGT").tm(0.05, 0) }.message!! shouldContain "oligo concentration must be more than 0"
        shouldThrow<BioException> { DNA("ACGT").tm(0.05, -1e-9) }
        DNA("GTAAAACGACGGCCAGT").tm(0, 250e-9, 0.002, 0).value.toDouble() shouldBe (61.13560696810322 plusOrMinus 1e-9)   // Mg2+ alone (Biopython 1.88)
        shouldThrow<BioException> { DNA("ACGN").tm() }.message!! shouldContain "'N'"
        shouldThrow<BioException> { DNA("ACGT", circular = true).tm() }
        shouldThrow<BioException> { DNA("ACGT").tm(0, 250e-9) }.message!! shouldContain "Na+ or Mg2+"
    }

    test("digest, sameMolecule, equality") {
        DNA("ACGT").digest shouldBe "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"       // refget's sha512t24u of ACGT
        DNA("ACGTA", circular = true).sameMolecule(DNA("TAACG", circular = true)) shouldBe true
        DNA("ACGTA", circular = true).sameMolecule(DNA("TAACG")) shouldBe false
        DNA("AACG", ds = true).sameMolecule(DNA("CGTT", ds = true)) shouldBe true
        DNA("AACG").sameMolecule(DNA("CGTT")) shouldBe false
        DNA("ACGT") shouldBe DNA("acg t")
        (DNA("ACGT") == DNA("ACGT", ds = true)) shouldBe false
        (DNA("ACGT").hashCode() == DNA("ACGT").hashCode()) shouldBe true
    }

    test("functions: seq detection and kinds, a260, hamming, levenshtein") {
        BioBag.seqOf("ACGU").shouldBeInstanceOf<RNA>()
        BioBag.seqOf("ACGT").shouldBeInstanceOf<DNA>()
        BioBag.seqOf("KKRR").shouldBeInstanceOf<Protein>()
        BioBag.seqOf("ACGT", "protein").shouldBeInstanceOf<Protein>()
        (BioBag.seqOf("ACGT", "dsdna") as DNA).ds shouldBe true
        BioBag.a260(1.0, DNA("ACGT")) shouldBe 33.0
        BioBag.a260(1.0, DNA("ACGT", ds = true)) shouldBe 50.0
        BioBag.a260(1.0, "rna") shouldBe 40.0
        BioBag.hamming("GATTACA", "GACTATA") shouldBe 2
        BioBag.levenshtein("KITTEN", "SITTING") shouldBe 3
        shouldThrow<BioException> { BioBag.hamming("A", "AA") }
        shouldThrow<BioException> { BioBag.nucleicOf("MKT", "revcomp") }.message!! shouldContain "not a protein"
    }

    test("the stored form packs 2, 4 or 5 bits a letter, and damage is refused") {
        val plasmid = DNA("ACGT".repeat(2500), ds = true, circular = true)
        Packing.encode(plasmid).size shouldBe 1 + 2 + 2500
        Packing.decode(Packing.encode(plasmid), 0) shouldBe plasmid
        Packing.encode(DNA("ACGN")).size shouldBe 1 + 1 + 2
        Packing.encode(Protein("MKTAY")).size shouldBe 1 + 1 + 4
        shouldThrow<IllegalArgumentException> { Packing.decode(Packing.encode(DNA("ACGT")).dropLast(1).toByteArray(), 0) }
        shouldThrow<IllegalArgumentException> { Packing.decode(Packing.encode(DNA("ACGT")), 1) }
        val kind3 = Packing.encode(DNA("ACGT")).also { it[0] = (it[0].toInt() or 3).toByte() }
        shouldThrow<IllegalArgumentException> { Packing.decode(kind3, 3) }               // what Sequence's decode passes
        val padded = Packing.encode(DNA("ACG")).also { it[it.size - 1] = (it[it.size - 1].toInt() or 1).toByte() }
        shouldThrow<IllegalArgumentException> { Packing.decode(padded, 0) }
    }
})
