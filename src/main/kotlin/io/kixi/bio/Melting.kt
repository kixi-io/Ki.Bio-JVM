package io.kixi.bio

import kotlin.math.ln
import kotlin.math.sqrt

/**
 * Melting temperature of a DNA oligo with its complement, as IDT's OligoAnalyzer
 * computes it (Ki type bags design, 5.5, decision 15):
 *
 * - Nearest-neighbor two-state thermodynamics with the unified DNA/DNA parameters
 *   (Allawi and SantaLucia 1997; SantaLucia 1998): the stacks' ΔH and ΔS, the
 *   terminal A·T or G·C initiation at each end, and the symmetry term for a
 *   self-complementary oligo. The Tm at 1 M Na+ is ΔH / (ΔS + R ln k).
 * - The salt correction of Owczarzy et al. 2008 (Biochemistry 47, 5336), applied
 *   as 1/Tm = 1/Tm(1 M) + correction: with no Mg2+, or little against the
 *   monovalent ions (√[Mg2+]/[Mon+] < 0.22), the monovalent correction of
 *   Owczarzy et al. 2004 (Biochemistry 43, 3537); otherwise the Mg2+ correction,
 *   its a, d and g adjusted for the monovalent ions below a ratio of 6. dNTPs bind
 *   Mg2+ (Ka 3 × 10^4), so the free Mg2+ is used.
 *
 * The constants and the order of the arithmetic follow Biopython's
 * `MeltingTemp.Tm_NN` (table `DNA_NN3`, salt correction 7), which the tests
 * compare against.
 */
internal object Melting {

    /** ΔH (kcal/mol) and ΔS (cal/(K·mol)) per stack, 5'→3' on the top strand. */
    private val NN: Map<String, DoubleArray> = buildMap {
        fun nn(a: String, b: String, h: Double, s: Double) { put(a, doubleArrayOf(h, s)); put(b, doubleArrayOf(h, s)) }
        nn("AA", "TT", -7.9, -22.2)
        nn("AT", "AT", -7.2, -20.4)
        nn("TA", "TA", -7.2, -21.3)
        nn("CA", "TG", -8.5, -22.7)
        nn("GT", "AC", -8.4, -22.4)
        nn("CT", "AG", -7.8, -21.0)
        nn("GA", "TC", -8.2, -22.2)
        nn("CG", "CG", -10.6, -27.2)
        nn("GC", "GC", -9.8, -24.4)
        nn("GG", "CC", -8.0, -19.9)
    }

    private const val R = 1.987            // cal/(K·mol), as Biopython and the 1998 tables use

    /**
     * Tm in °C of [letters] (plain A, C, G, T; at least two) with its complement.
     * Concentrations in mol/L: [mon] monovalent cations (Na+, K+), [mg] Mg2+,
     * [dntp] dNTPs, [oligo] the oligo. The complement is taken as much less
     * concentrated than the oligo, as for a primer on its template, so
     * k = [oligo]; a self-complementary oligo uses k = [oligo] too, with the
     * symmetry term.
     */
    fun tm(letters: String, mon: Double, mg: Double, dntp: Double, oligo: Double): Double {
        require(letters.length >= 2)
        var dH = 0.0
        var dS = 0.0
        // Terminal base pairs.
        for (end in charArrayOf(letters.first(), letters.last())) {
            if (end == 'G' || end == 'C') { dH += 0.1; dS += -2.8 } else { dH += 2.3; dS += 4.1 }
        }
        for (i in 0 until letters.length - 1) {
            val p = NN[letters.substring(i, i + 2)] ?: throw BioException("No nearest-neighbor parameters for '${letters.substring(i, i + 2)}'")
            dH += p[0]; dS += p[1]
        }
        val selfComplementary = letters == reverseComplement(letters)
        if (selfComplementary) dS += -1.4
        val k = oligo * OLIGO_FACTOR
        val tm1M = (1000.0 * dH) / (dS + R * ln(k))                       // kelvin, at 1 M Na+
        val corr = saltCorrection(letters, mon, mg, dntp)
        return 1.0 / (1.0 / tm1M + corr) - 273.15
    }

    /**
     * How the oligo concentration enters ln k for a non-self-complementary oligo:
     * 1 when the complement is much less concentrated (a primer on its target).
     * IDT does not state its convention; this is checked against values read from
     * OligoAnalyzer (see the tests).
     */
    internal const val OLIGO_FACTOR: Double = 1.0

    /** Owczarzy 2008 (with 2004 for monovalent ions), as 1/Tm(1 M) + correction. */
    fun saltCorrection(letters: String, mon: Double, mg0: Double, dntp: Double): Double {
        val fGC = letters.count { it == 'G' || it == 'C' }.toDouble() / letters.length
        var mg = mg0
        if (dntp > 0 && mg > 0) {
            val ka = 3e4
            mg = (-(ka * dntp - ka * mg + 1.0) + sqrt((ka * dntp - ka * mg + 1.0).let { it * it } + 4.0 * ka * mg)) / (2.0 * ka)
        }
        var a = 3.92; val b = -0.911; val c = 6.26; var d = 1.42
        val e = -48.2; val f = 52.5; var g = 8.31
        if (mon > 0) {
            val ratio = sqrt(mg) / mon
            if (ratio < 0.22) {
                val l = ln(mon)
                return (4.29 * fGC - 3.95) * 1e-5 * l + 9.40e-6 * l * l
            } else if (ratio < 6.0) {
                val l = ln(mon)
                a = 3.92 * (0.843 - 0.352 * sqrt(mon) * l)
                d = 1.42 * (1.279 - 4.03e-3 * l - 8.03e-3 * l * l)
                g = 8.31 * (0.486 - 0.258 * l + 5.25e-3 * l * l * l)
            }
        }
        if (mg <= 0) throw BioException("A melting temperature needs Na+ (or K+) or Mg2+")
        val lm = ln(mg)
        return (a + b * lm + fGC * (c + d * lm) + (1.0 / (2.0 * (letters.length - 1))) * (e + f * lm + g * lm * lm)) * 1e-5
    }

    private fun reverseComplement(s: String): String {
        val sb = StringBuilder(s.length)
        for (i in s.indices.reversed()) sb.append(when (s[i]) { 'A' -> 'T'; 'T' -> 'A'; 'C' -> 'G'; 'G' -> 'C'; else -> s[i] })
        return sb.toString()
    }
}
