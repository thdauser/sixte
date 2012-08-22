/** \file pil_error.h
    \brief Declaration of Ape's PIL emulating error codes and error reporting utilities.

    Ape is backward compatible with a subset of the PIL library developed by
    ISDC. Replacements are provided for the most commonly used API functions.
    Less common functions and functions which rely on details of pil's internal
    structures are not provided. Ape's implementation of these functions does
    not always behave identically to PIL's implementation, but the behavior should be
    close enough to provide compatibility.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "pil_error.h"

static const char * PIL_err_table[] = {
  "NULL pointer passed",
  "bad argument(s) passed",
  "no enough memory",
  "file does not exist",
  "error reading file",
  "error writing file",
  "end of string encountered",
  "bad name",
  "bad type",
  "bad mode",
  "bad line",
  "not implemented",
  "file does not exist",
  "file exists",
  "file not readable",
  "file not writable",
  "blank line",
  "comment line",
  "error in line format",
  "not found",
  "PFILES env. variable too long",
#ifdef WIN32
  "PFILES env. missing '|'",
#else
  "PFILES env. missing ';'",
#endif
  "cannot lock/unlock parameter for exclusive access",
  "bogus parameters in command line",
  "no logger function",
  "too many tokens found in the line",
  "not enough tokens found in the line",
  "unmatched quote/doublequote found in the line",
  "line without terminating character",
  "extra space after trailing quote",
  "cannot convert string to boolean",
  "cannot convert string to integer",
  "cannot convert string to real",
  "cannot convert string to a var. vector of integers",
  "cannot convert string to a vector of integers",
  "cannot convert string to a var. vector of reals",
  "cannot convert string to a vector of reals",
  "value is out of range set by min/max fields",
  "value does not match any in enumerated list (min field)",
  "file not found (or has wrong access type)",
  "problem converting string to other type",
  "value is \"undefined\"",
  "non-specific pil error"
};

const char * PIL_err_handler(int r) {
  const char * result = 0;
  if (PIL_ERR_MIN_IDX <= r && PIL_ERR_MAX_IDX >= r)
    result = PIL_err_table[PIL_ERR_MAX_IDX - r];
  return result;
}

/*
 * $Log: pil_error.c,v $
 * Revision 1.4  2006/06/17 03:05:54  peachey
 * Correct spelling in error message.
 *
 * Revision 1.3  2006/06/06 17:47:32  peachey
 * Add some error codes to cover essential ape errors which have no
 * counterpart in pil error code space.
 *
 */
