#ifndef SIXT_H
#define SIXT_H 1


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "headas_rand.h"


#define FILENAME_LENGTH 512 /**< Maximum length of a filename. */
#define MAXMSG 512 /**< Maximum length of a message string. */
#define HEADER_LENGTH 288001 /**< Maximum length of a FITS header string. */

#define SIXT_NAN (0./0.) /**< Not a Number. */


#endif /* SIXT_H */

