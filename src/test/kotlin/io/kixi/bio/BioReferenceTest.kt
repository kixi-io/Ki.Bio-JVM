package io.kixi.bio

import io.kixi.uom.Unit as KiUnit
import io.kixi.uom.Quantity
import io.kotest.assertions.withClue
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldBeEmpty
import io.kotest.matchers.shouldBe
import java.math.BigDecimal
import kotlin.math.abs

/**
 * Ki.Bio against Biopython 1.88, which implements the methods Ki.Bio follows
 * (`src/test/resources/io/kixi/bio/biopython-reference.tsv`, written by
 * `tools/reference.py`): Tm as `MeltingTemp.Tm_NN` with the Allawi and SantaLucia
 * table and Owczarzy 2008 salt correction, pI as `IsoelectricPoint` (ExPASy's
 * method), translation under every NCBI code, and average masses.
 */
class BioReferenceTest : FunSpec({

    val rows: List<List<String>> = BioReferenceTest::class.java.getResourceAsStream("biopython-reference.tsv")!!
        .bufferedReader().readLines().filter { it.isNotBlank() && !it.startsWith("#") }.map { it.split('\t') }

    test("Tm agrees with Biopython to 0.001 °C for 131 oligos under 9 buffers") {
        val bad = ArrayList<String>()
        var n = 0
        for (r in rows.filter { it[0] == "tm" }) {
            val (seq, na, mg, dntp, oligo) = r.subList(1, 6)
            val want = r[6].toDouble()
            val got = DNA(seq).tm(Quantity(BigDecimal(na), KiUnit.mM), Quantity(BigDecimal(mg), KiUnit.mM),
                Quantity(BigDecimal(dntp), KiUnit.mM), Quantity(BigDecimal(oligo), KiUnit.nM)).value.toDouble()
            n++
            if (abs(got - want) > 0.001) bad.add("$seq Na $na Mg $mg dNTP $dntp oligo $oligo: $got, Biopython $want")
        }
        withClue(bad.take(10).joinToString("\n")) { bad.shouldBeEmpty() }
        n shouldBe 1179
    }

    test("pI and charge agree with Biopython (ExPASy's method) for 160 proteins") {
        val bad = ArrayList<String>()
        for (r in rows.filter { it[0] == "pI" }) {
            val p = Protein(r[1])
            val want = r[2].toDouble()
            // Biopython searches pH 4.05 to 12 only, so a very acidic or basic protein comes back pinned at a bound;
            // ExPASy's code, and Ki.Bio, search 0 to 14. There, Ki.Bio's pI must lie beyond the bound.
            if (abs(want - 4.05) < 0.001 || abs(want - 12.0) < 0.001) {
                if (if (want < 5) p.pI!! > 4.05 else p.pI!! < 12.0) bad.add("${r[1]}: pI ${p.pI} inside Biopython's bound $want")
                continue
            }
            // Both stop at 0.0001, halving differently: within the precision, and equal to 2 decimals as ExPASy prints.
            if (abs(p.pI!! - want) > 0.0002) bad.add("${r[1]}: pI ${p.pI}, Biopython $want")
            if (abs(p.charge(7.0)!! - r[3].toDouble()) > 1e-9) bad.add("${r[1]}: charge ${p.charge(7.0)}, Biopython ${r[3]}")
        }
        withClue(bad.take(10).joinToString("\n")) { bad.shouldBeEmpty() }
    }

    test("translation agrees with Biopython under all 27 NCBI codes, in every frame") {
        val bad = ArrayList<String>()
        val codes = HashSet<Int>()
        for (r in rows.filter { it[0] == "code" }) {
            val code = r[1].toInt(); codes += code
            val got = DNA(r[3]).translate(r[2].toInt(), code).letters
            if (got != r[4]) bad.add("code $code frame ${r[2]} ${r[3]}: $got, Biopython ${r[4]}")
        }
        withClue(bad.take(10).joinToString("\n")) { bad.shouldBeEmpty() }
        codes shouldBe GeneticCodes.all.map { it.id }.toSet()
        codes.size shouldBe 27
    }

    test("average masses agree with Biopython's within the difference of their atomic weights") {
        // Biopython's DNA and RNA have a 5' phosphate (monophosphates less water per bond); Ki.Bio's a 5' hydroxyl,
        // as an oligo is synthesized (IDT): HPO3 less. Biopython uses older atomic weights; the difference stays
        // under 0.002 % (about 0.0006 Da a nucleotide, 0.003 Da an amino acid).
        val hpo3 = 79.9799
        val bad = ArrayList<String>()
        for (r in rows.filter { it[0] == "mass" }) {
            val s = when (r[1]) { "DNA" -> DNA(r[2]); "RNA" -> RNA(r[2]); else -> Protein(r[2]) }
            val want = r[3].toDouble() - if (s is Protein) 0.0 else hpo3
            val got = s.mass!!.value.toDouble()
            if (abs(got - want) > 2e-5 * want + 0.01) bad.add("${r[1]} ${r[2].take(20)}: $got, Biopython less HPO3 $want")
        }
        withClue(bad.take(10).joinToString("\n")) { bad.shouldBeEmpty() }
    }
})
