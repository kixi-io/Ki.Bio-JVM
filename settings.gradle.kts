// Ki.Bio-JVM/settings.gradle.kts

rootProject.name = "Ki.Bio-JVM"

// Composite build, as in the other Ki repos: edits to the sibling checkouts are
// live when building from here. Ki.Core and Ki.Sci are the library's
// dependencies (Ki.Sci for formulas and masses); Ki.KD, Ki.KS and RevDB are here
// for the tests that put the bag's values through every layer (the corpus
// round trip). The versions in build.gradle.kts are the published fallback for
// builds without local checkouts.
includeBuild("../Ki.Core-JVM")
includeBuild("../Ki.Sci-JVM")
includeBuild("../Ki.KD-JVM")
includeBuild("../Ki.KS-JVM")
includeBuild("../RevDB-JVM")
