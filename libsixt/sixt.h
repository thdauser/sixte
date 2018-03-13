/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef SIXT_H
#define SIXT_H 1

#include <sixteconfig.h>

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <expat.h>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ape/ape_trad.h"
#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "wcshdr.h"

#include "rmf.h"
#include "parinput.h"

/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/** Maximum length of a filename in cfitsio. */
#define MAXFILENAME FLEN_FILENAME
/** Maximum length of a message string. */
#define MAXMSG (512)

/** Not a Number. */
#define SIXT_NAN (0./0.)


/** MJDREF used in the FITS header of eROSITA event files
    [d]. Corresponds to 2000-01-01T00:00:00. */
extern const double eromjdref;

/** MJDREF used in the FITS header of XMM event files [d]. Corresponds
    to 1998-01-01T00:00:00.00. */
extern const double xmmmjdref;

/////////////////////////////////////////////////////////////////
extern int sixt_argc;
extern char **sixt_argv;

/////////////////////////////////////////////////////////////////
// Macro definitions.
/////////////////////////////////////////////////////////////////


/** Returns the maximum of 2 values. */
#define MAX(a, b) ( (a)>(b) ? (a) : (b) )
/** Returns the minimum of 2 values. */
#define MIN(a, b) ( (a)<(b) ? (a) : (b) )

// Flag for deprecated functions.
// Note: these functions should also make use of SIXT_DEPRECATED().
#define DEPRECATED(func) func __attribute__ ((deprecated))

// Error handling macros.
#define SIXT_ERROR(msg) (sixt_error(__func__, msg))

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_STATUS_VOID(status) \
  if (EXIT_SUCCESS!=status) return;

/** FITS error message for error handling macro retrieving error status */
extern char _fits_err_msg[80];
#define CHECK_STATUS_BREAK_WITH_FITSERROR(status) \
  if (EXIT_SUCCESS!=status){ \
	fits_get_errstatus(status,_fits_err_msg); \
	SIXT_ERROR(_fits_err_msg); \
	status=EXIT_FAILURE; \
	break; \
  }


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

#define CHECK_MALLOC_RET_NULL(a) \
		if (NULL==a) { \
			SIXT_ERROR("memory allocation failed"); \
			return NULL;\
		}

#define CHECK_MALLOC_RET_STATUS(a,b,status) \
		if (NULL==a) { \
			SIXT_ERROR("memory allocation failed"); \
			status=EXIT_FAILURE; \
			return b;\
		}


#define CHECK_MALLOC_RET_NULL_STATUS(a,status) \
		if (NULL==a) { \
			SIXT_ERROR("memory allocation failed"); \
			status=EXIT_FAILURE; \
			return NULL;\
		}

#define CHECK_MALLOC_VOID_STATUS(a,status) \
		if (NULL==a) { \
			SIXT_ERROR("memory allocation failed"); \
			status=EXIT_FAILURE; \
			return;\
		}


#define CHECK_MALLOC_VOID(a) \
		if (NULL==a) { \
			SIXT_ERROR("memory allocation failed"); \
			return;\
		}

// Warnings.
#define SIXT_WARNING(msg) (sixt_warning(msg))
#define SIXT_DEPRECATED(fnc, alt) (sixt_deprecated(fnc, alt))

/////////////////////////////////////////////////////////////////
// Type declarations.
/////////////////////////////////////////////////////////////////

typedef struct {
	/** Telescope keyword */
	char* telescop;

	/** Instrument keyword */
	char* instrume;

	/** Filter keyword */
	char* filter;

	/** Ancillary response file */
	char* ancrfile;

	/** Response file */
	char* respfile;

	/** Extension name */
	char* extname;

	/** Reference MJD*/
	double mjdref;

	/** Time offset */
	double timezero;

	/** Start time */
	double tstart;

	/** Stop time */
	double tstop;

} SixtStdKeywords;

/////////////////////////////////////////////////////////////////
// Function declarations.
/////////////////////////////////////////////////////////////////


/** Return a seed for the random number generator. */
unsigned int getSeed(int seed);

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
void sixt_init_rng(const unsigned int seed, int* const status);

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

/** Print a message that this function is deprecated and propose
    the given alternative (if supplied). */
void sixt_deprecated(const char* const fnc, const char* const alt);

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

/** Determine the eROSITA XMLFilename if different 
 *telescopes are used. */
void sixt_get_eroXMLFile(char *filename,
			const int telescop_index,
			int* const status);

/** Determine a date and a time string for the specified MJDREF offset
    and time given in seconds since MJDREF. */
void sixt_get_date_time(const double mjdref,
			const double t,
			char* const datestr,
			char* const timestr,
			int* const status);

/** Add standard FITS header keywords to the specified file. */
void sixt_add_fits_stdkeywords_obsolete(fitsfile* const fptr,
			       const int hdunum,
			       char* const telescop,
			       char* const instrume,
			       char* const filter,
			       char* const ancrfile,
			       char* const respfile,
			       double mjdref,
			       double timezero,
			       double tstart,
			       double tstop,
			       int* const status);

/** Reads standard header keywords from a FITS file. */
void sixt_read_fits_stdkeywords_obsolete(fitsfile* const ifptr,
			       char* const telescop,
			       char* const instrume,
			       char* const filter,
			       char* const ancrfile,
			       char* const respfile,
			       double *mjdref,
			       double *timezero,
			       double *tstart,
			       double *tstop, 
			       int* const status);

/** Add eROSITA-specific standard FITS header keywords to the
    specified file. */
void sixt_add_fits_erostdkeywords(fitsfile* const fptr, 
				  const int hdunum,
				  char* const creation_date,
				  char* const date_obs,
				  char* const time_obs,
				  char* const date_end,
				  char* const time_end,
				  double tstart,
				  double tstop,
				  double mjdref,
				  double timezero,
				  int ccdnr,
				  int* const status);

/** Determine whether the given value for MJDREF is equivalent to the
    specified reference MJDREF. In order to allow a better
    localization of the problem by the user, an optional description
    can be specified, which is added to the displayed message in case
    of a mismatch. */
void verifyMJDREF(const double refmjdref,
		  const double mjdref,
		  const char* const description,
		  int* const status);

/** Make sure that the value of TIMEZERO is '0.0'. This is required by
    the current implementation. */
void verifyTIMEZERO(const double timezero,
		    int* const status);

/** Determine the signal corresponding to a particular PHA channel
    according to the EBOUNDS table. The input channel must have the
    same offset as in the EBOUNDS table. I.e. if the first channel in
    the EBOUNDS has the number 1, the numbering starts at 1. If the
    first channel has the number 0, the numbering starts at 0.  The
    returned energy is randomized between the lower and the upper bin
    boundary and is given in the same units as the EBOUNDS, (usually
    [keV]). */
float getEBOUNDSEnergy(const long channel,
		       const struct RMF* const rmf, 
		       int* const status);

/** Add standard FITS header keywords to the specified file using info
 *  contained in a SixtStdKeywords structure. */
void sixt_add_fits_stdkeywords(fitsfile* const fptr,
		const int hdunum,
		SixtStdKeywords * keyword_struct,
		int* const status);

/** Reads standard header keywords from a FITS file using a
 *  SixtStdKeywords structure. Does at the same time the
 *  malloc of the different char arrays. */
void sixt_read_fits_stdkeywords(fitsfile* const ifptr,
		SixtStdKeywords* keyword_struct,
		int* const status);

/** Constructor of the SixtStdKeywords structure: returns a pointer to an empty structure of this type */
SixtStdKeywords* newSixtStdKeywords(int* const status);

/** Builds a SixtStdKeywords struct from the individual keywords.
 * 	Does at the same time the malloc of the different char arrays. */
SixtStdKeywords* buildSixtStdKeywords(char* const telescop,
	       char* const instrume,
	       char* const filter,
	       char* const ancrfile,
	       char* const respfile,
	       char* const extname,
	       double mjdref,
	       double timezero,
	       double tstart,
	       double tstop,
	       int* const status);

/** copies a SixtStdKeywords structure (in case many keywords are to be duplicated) **/
SixtStdKeywords* duplicateSixtStdKeywords(const SixtStdKeywords *key, int* const status);


/** Destructor of the SixtStdKeywordsStructure */
void freeSixtStdKeywords(SixtStdKeywords* keyword_struct);

/** convenience function to create a FITS-file or to error out depending on the value of **/
/** clobber **/
int fits_create_file_clobber(fitsfile **fptr, char *filename, int clobber, int *status);

// convenience function: update checksum in current and primary HDU and close the file
void fits_close_file_chksum(fitsfile *fptr,int *status);

// check for obsolete input keywords
void sixt_check_obsolete_keyword(int* status);

#endif /* SIXT_H */

