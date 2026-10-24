package io.kixi.bio

import io.kixi.KiException

/** Text Ki.Bio cannot read as a sequence, or a sequence it cannot work with (Tm of an ambiguous primer, a frame of 4). */
public open class BioException(message: String) : KiException(message)
