#ifndef SIXT_H
#define SIXT_H 1

#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <expat.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "headas_rand.h"

#include "sixt_random.h"
#include "sixt_string.h"


#define FILENAME_LENGTH 512 /**< Maximum length of a filename. */
#define MAXMSG 512 /**< Maximum length of a message string. */

#define SIXT_NAN (0./0.) /**< Not a Number. */

/** Seed for the HEAdas random number generator. */
#define SIXT_HD_RANDOM_SEED 59843


/** Macro returning the maximum of 2 values. */
#define MAX(a, b) ( (a)>(b) ? (a) : (b) )
/** Macro returning the minimum of 2 values. */
#define MIN(a, b) ( (a)<(b) ? (a) : (b) )


#endif /* SIXT_H */

