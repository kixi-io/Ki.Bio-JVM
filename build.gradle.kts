// Ki.Bio-JVM/build.gradle.kts

plugins {
    `java-library`
    `maven-publish`
    kotlin("jvm") version "2.3.0"
}

group = "io.kixi"
version = "0.1.0"
description = "Ki.Bio: the `bio` type bag (DNA, RNA and protein sequences) for KD, KS, KQ and RevDB"

repositories {
    mavenCentral()
}

dependencies {
    // Ki.Core, and Ki.Sci for formulas and masses (Ki type bags design, 5.5).
    api("io.kixi:Ki.Core-JVM:2.5.0")
    api("io.kixi:Ki.Sci-JVM:0.1.0")

    // The corpus round trip: KD text, KS literals and members, RevDB storage and indexes.
    testImplementation("io.kixi:Ki.KD-JVM:2.4.0")
    testImplementation("io.kixi:Ki.KS-JVM:2.4.0")
    testImplementation("io.kixi:RevDB-JVM:0.9.0")

    testImplementation("io.kotest:kotest-runner-junit5:6.0.7")
    testImplementation("io.kotest:kotest-assertions-core:6.0.7")
    // Kotest 6 asks for kotlin-reflect 2.2.x; pin it to the compiler's version, as Ki.KS does.
    testImplementation(kotlin("reflect"))
}

kotlin {
    jvmToolchain(21)
    explicitApi()
}

java {
    withSourcesJar()
}

tasks.test {
    useJUnitPlatform()
    maxHeapSize = "2g"
    testLogging {
        events("passed", "skipped", "failed")
        showStandardStreams = true
    }
}

// ============================================================================
// Publishing (GitHub Packages, as the other Ki libraries)
// ============================================================================

publishing {
    repositories {
        maven {
            name = "GitHubPackages"
            url = uri("https://maven.pkg.github.com/kixi-io/Ki.Bio-JVM")
            credentials {
                username = providers.gradleProperty("gpr.user").orNull
                    ?: System.getenv("GITHUB_ACTOR")
                password = providers.gradleProperty("gpr.key").orNull
                    ?: System.getenv("GITHUB_TOKEN")
            }
        }
    }
    publications {
        create<MavenPublication>("mavenJava") {
            from(components["java"])

            pom {
                groupId = "${project.group}"
                artifactId = "Ki.Bio-JVM"
                name.set("Ki.Bio")
                description.set(project.description)
                url.set("https://kixi.io")

                licenses {
                    license {
                        name.set("MIT License")
                        url.set("https://raw.githubusercontent.com/kixi-io/Ki.Bio-JVM/master/LICENSE")
                        distribution.set("repo")
                    }
                }
                developers {
                    developer {
                        id.set("dleuck")
                        name.set("Daniel Leuck")
                    }
                }
                scm {
                    connection.set("scm:git:git://github.com/kixi-io/Ki.Bio-JVM.git")
                    developerConnection.set("scm:git:ssh://github.com/kixi-io/Ki.Bio-JVM.git")
                    url.set("https://github.com/kixi-io/Ki.Bio-JVM/")
                }
            }
        }
    }
}
