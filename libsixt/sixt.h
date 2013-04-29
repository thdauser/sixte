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
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <expat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ape/ape_trad.h"
#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "wcshdr.h"



/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/** Maximum length of a filename. */
#define MAXFILENAME (512)
/** Maximum length of a message string. */
#define MAXMSG (512)

/** Not a Number. */
#define SIXT_NAN (0./0.)


/////////////////////////////////////////////////////////////////
// Macro definitions.
/////////////////////////////////////////////////////////////////


/** Returns the maximum of 2 values. */
#define MAX(a, b) ( (a)>(b) ? (a) : (b) )
/** Returns the minimum of 2 values. */
#define MIN(a, b) ( (a)<(b) ? (a) : (b) )


// Error handling macros.
#define SIXT_ERROR(msg) (sixt_error(__func__, msg))

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_STATUS_VOID(status) \
  if (EXIT_SUCCESS!=status) return;


#define CHECK_NULL_VOID(a,status,msg) \
  if (NULL==a) { \
    SIXT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return;\
  }

#define CHECK_NULL_BREAK(a,status,msg) \
  if (NULL==a) { \
    SIXT_ERROR(msg); \
    status=EXIT_FAILURE; \
    break;\
  }

#define CHECK_NULL_RET(a,status,msg,ret) \
  if (NULL==a) { \
    SIXT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return(ret);	 \
  }

#define CHECK_NULL(a,status,msg) CHECK_NULL_RET(a,status,msg,NULL);


// Warnings.
#define SIXT_WARNING(msg) (sixt_warning(msg))


/////////////////////////////////////////////////////////////////
// Function declarations.
/////////////////////////////////////////////////////////////////


/** This routine returns a random number. The values are either
    obtained from the Remeis random number server or are created by
    the HEAdas random number generator. The routine is basically a
    wrapper around the respective library routines, either the
    rcl_rand_ndg() or the HEAdas routine HDmtDrand(). The return value
    lies in the interval [0,1). 

    In case the HEAdas random number generator is used, the function
    requires HDmtInit() to be called once before usage for
    initialization. When the HEAdas random number generator is not
    needed any more, it can be realeased with HDmtFree(). Information
    can be found in the HEAdas developer's guide or directly in the
    source files 'headas_rand.h' and 'headas_rand.c'. */
double sixt_get_random_number(int* const status);

/** Initialize the random number generator. */
void sixt_init_rng(const int seed, int* const status);

/** Clean up the random number generator. */
void sixt_destroy_rng();

/** This routine produces two Gaussian distributed random numbers. The
    standard deviation of the Gaussian distribution sigma is assumed
    to be unity. The two numbers are returned via the pointer function
    arguments. */
void sixt_get_gauss_random_numbers(double* const x, 
				   double* const y,
				   int* const status);

/** Returns a random value on the basis of an exponential distribution
    with a given average distance. In the simulation this function is
    used to calculate the temporal differences between individual
    photons from a source. The photons have Poisson statistics. */
double rndexp(const double avg, int* const status);

/** Convert a squence of chars into captial letters. The sequence has
    to be terminated by a '\0' mark. */
void strtoupper(char* const string);

/** Print the given error message for an error occured in the
    specified function. The function name is also part of the
    output. */
void sixt_error(const char* const func, const char* const msg);

/** Print the given warning message. */
void sixt_warning(const char* const msg);

/** Determine the XMLFilename according to the selected mission,
    instrument and mode. */
void sixt_get_XMLFile(char* const filename,
		      const char* const xmlfile,
		      const char* const mission,
		      const char* const instrument,
		      const char* const mode,
		      int* const status);

/** Determine the LAD XMLFilename. */
void sixt_get_LADXMLFile(char* const filename,
			 const char* const xmlfile);

void sixt_add_fits_stdkeywords(fitsfile* const fptr,
			       const int hdunum,
			       char* const telescop,
			       char* const instrume,
			       char* const filter,
			       double mjdref,
			       double timezero,
			       double tstart,
			       double tstop,
			       int* const status);


#endif /* SIXT_H */

