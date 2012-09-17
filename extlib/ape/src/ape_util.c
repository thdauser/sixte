/** \file ape_util.c
    \brief Implementation of internal utilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_msg.h"
#include "ape/ape_util.h"
#include "ape/ape_test.h"

#ifdef WIN32
#include <windows.h>
static const DWORD s_file_not_found = 0xFFFFFFFF;
#else
#include <unistd.h>
#endif

#ifdef USE_READLINE
#include "readline/history.h"
#include "readline/readline.h"
#endif

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
#define PFILES_DIVIDER "|"
#define QP "\\"
#else
#define PFILES_DIVIDER ";"
#define QP "/"
#endif
static const char s_qp[] = QP;
static const char s_pfiles_divider[] = PFILES_DIVIDER;
static FILE ** get_prompt_stream_ptr(void);

/* Case insensitive string comparison. */
static int ci_strcmp(const char * s1, const char * s2) {
  return ape_util_strcmp(s1, s2, 1);
}

/* Return a copy of a string with leading and trailing white space removed. Client must free the string. */
static char * truncate_white_space(const char * begin) {
  char * output = 0;
  if (0 != begin) {
    const char * end = 0;

    /* Skip leading white space. */
    while (0 != isspace(*begin)) ++begin;

    /* Skip trailing white space. */
    for (end = begin + strlen(begin); end != begin && 0 != isspace(*(end - 1)); --end) {}

    /* Create output buffer, using calloc so that the contents are initialized to 0. */
    output = (char *) calloc(end - begin + 1, sizeof(char));

    if (0 != output) {
      /* Copy truncated input to output buffer. */
      strncpy(output, begin, end - begin);
    }
  }
  return output;
}

/* Strictly convert input string to output bool, ignoring leading/trailing white space, and flagging every possible error. */
static int string2bool(const char * input, char * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != input) {
      /* Truncate white space from input string. */
      char * copy = truncate_white_space(input);

      /* Compare input string (case insensitive) to various boolean expressions. */
      if (0 == ci_strcmp(copy, "true") || 0 == ci_strcmp(copy, "t") ||
        0 == ci_strcmp(copy, "yes") || 0 == ci_strcmp(copy, "y")) {
        *result = 1;
      } else if (0 == ci_strcmp(copy, "false") || 0 == ci_strcmp(copy, "f") ||
        0 == ci_strcmp(copy, "no") || 0 == ci_strcmp(copy, "n")) {
        *result = 0;
      } else {
        /* String does not contain a valid bool expression. */
        status = eConversionError;
      }

      /* Clean up. */
      free(copy); copy = 0;
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

/* Strictly convert input string to output double, ignoring leading/trailing white space, and flagging every possible error. */
static int string2double(const char * input, double * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0.;
    if (0 != input) {
      char * remainder = 0;
      /* Skip leading white space. */
      while (0 != isspace(*input)) ++input;

      if ('\0' == *input || 0 == ci_strcmp(input, "indef") || 0 == ci_strcmp(input, "none") || 0 == ci_strcmp(input, "undef") ||
        0 == ci_strcmp(input, "undefined")) {
        status = eUndefinedValue;
      } else if (0 == ci_strcmp(input, "inf") || 0 == ci_strcmp(input, "infinity") || 0 == ci_strcmp(input, "nan")) {
        status = eNan;
      } else {
        errno = 0;
        *result = strtod(input, &remainder);
        if (ERANGE == errno) {
          if (0. == *result) status = eUnderflow;
          else status = eOverflow;
        } else if (0 != remainder && '\0' != *remainder) {
          /* Ignore trailing whitespace; the C standard is not specific about whether it's returned or eaten by the conversion. */
          while (0 != isspace(*remainder)) ++remainder;
          if ('\0' != *remainder) status = eStringRemainder;
          /* else conversion succeeded. */
        } else {
          /* Conversion succeeded! */
        }
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

/* Strictly convert input string to output long, ignoring leading/trailing white space, and flagging every possible error. */
static int string2long(const char * input, long * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != input) {
      char * remainder = 0;
      errno = 0;
      /* Skip leading white space. */
      while (0 != isspace(*input)) ++input;

      if ('\0' == *input || 0 == ci_strcmp(input, "indef") || 0 == ci_strcmp(input, "none") || 0 == ci_strcmp(input, "undef") ||
        0 == ci_strcmp(input, "undefined")) {
        status = eUndefinedValue;
      } else if (0 == ci_strcmp(input, "inf") || 0 == ci_strcmp(input, "infinity") || 0 == ci_strcmp(input, "nan")) {
        *result = LONG_MAX;
      } else {
        *result = strtol(input, &remainder, 0);
        if (ERANGE == errno) {
          if (LONG_MIN == *result) status = eUnderflow;
          else if (LONG_MAX == *result) status = eOverflow;
          else status = eConversionError;
        } else if (0 != remainder && '\0' != *remainder) {
          /* Ignore trailing whitespace; the C standard is not specific about whether it's returned or eaten by the conversion. */
          while (0 != isspace(*remainder)) ++remainder;
          if ('\0' != *remainder) status = eStringRemainder;
          /* else conversion succeeded. */
        } else {
          /* Conversion succeeded! */
        }
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

static char * copy_range(const char * begin, const char * end) {
  char * result = 0;
  if (0 != begin && begin <= end) {
    size_t len = end - begin;
    result = (char *) calloc(len + 1, sizeof(char));
    if (0 != result) {
      strncpy(result, begin, len);
    }
  }
  return result;
}

static char * copy_string(const char * orig) {
  return copy_range(orig, orig + strlen(orig));
}

#define MAX_NUM_ATEXIT_FUNC 32
typedef void (*at_exit_func_ptr_type)(void);
static at_exit_func_ptr_type s_atexit[MAX_NUM_ATEXIT_FUNC];
static size_t s_atexit_idx = MAX_NUM_ATEXIT_FUNC;
static char s_atexit_registered = 0;

static void ape_all_atexit(void) {
  size_t idx = 0;
  /* Call all ape registered functions. */
  for (idx = s_atexit_idx; idx != MAX_NUM_ATEXIT_FUNC; ++idx) {
    if (0 != s_atexit[idx]) {
      s_atexit[idx]();
    }
  }
}

int ape_util_atexit(void (*func)(void)) {
  int status = eOK;

  /* Check input. */
  if (0 == func) {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* The first time this is called, register ape_all_atexit. */
    if (0 == s_atexit_registered) {
      memset(s_atexit, '\0', sizeof(s_atexit) / sizeof(char));
      if (0 == atexit(&ape_all_atexit)) {
        s_atexit_registered = 1;
      } else {
        status = eAtExitError;
      }
    }
  }

  if (eOK == status) {
    /* Add functions in reverse order so that at exit they will be called in reverse order. */
    if (0 < s_atexit_idx) {
      size_t idx = 0;
      /* See if this function was already added to the container. */
      for (; idx != MAX_NUM_ATEXIT_FUNC && func != s_atexit[idx]; ++idx) {}
      if (MAX_NUM_ATEXIT_FUNC == idx) {
        /* Function was not already added, so add it. */
        --s_atexit_idx;
        s_atexit[s_atexit_idx] = func;
      }
    } else {
      status = eTooManyAtExitFunc;
    }
  }

  return status;
}

static int cat_range(const char * s1, const char * s2begin, const char * s2end, char ** result) {
  int status = eOK;
  const char * s2 = s2begin;
  if (0 != result) {
    /* Initialize output. */
    *result = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    size_t len = 1u; /* For 0 termination. */
    if (0 != s1) len += strlen(s1);
    if (0 != s2) len += strlen(s2);

    /* Allocate output string, setting it to 0. */
    *result = (char *) calloc(len, sizeof(char));

    if (0 != *result) {
      /* Copy input(s) to output. */
      if (0 != s1) strcat(*result, s1);
      if (0 != s2) {
        if (s2begin <= s2end) strncat(*result, s2, s2end - s2begin);
        else strcat(*result, s2);
      }
    } else {
      status = eDynAllocFailed;
    }
  }

  return status;
}

int ape_util_cat_string(const char * s1, const char * s2, char ** result) {
  return cat_range(s1, s2, 0, result);
}

int ape_util_append_file_name(const char * dir_name, const char * file_name, char ** full_file_name) {
  int status = eOK;
  const size_t qp_len = strlen(s_qp);
  size_t full_file_len = qp_len + 1;

  if (0 != full_file_name) {
    *full_file_name = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Compute length needed for full name. */
    if (0 != dir_name) full_file_len += strlen(dir_name);
    if (0 != file_name) full_file_len += strlen(file_name);

    /* Allocate full file name buffer. */
    *full_file_name = calloc(full_file_len, sizeof(const char *));
    if (0 != *full_file_name) {
      strcat(*full_file_name, dir_name);
      strcat(*full_file_name, s_qp);
      strcat(*full_file_name, file_name);
    } else {
      status = eDynAllocFailed;
    }
  }

  return status;
}

int ape_util_copy_range(const char * begin, const char * end, char ** result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 == begin || 0 == end) {
      status = eNullPointer;
    } else if (begin > end) {
      status = eInvalidArgument;
    } else {
      *result = copy_range(begin, end);
      if (0 == *result) {
        status = eDynAllocFailed;
      }
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_util_copy_string(const char * orig, char ** result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != orig) {
      *result = copy_string(orig);
      if (0 == *result) {
        status = eDynAllocFailed;
      }
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

void ape_util_free_string_array(char ** string_array) {
  if (0 != string_array) {
    char ** sp = string_array;
    /* Go to one past the last pointer in the array. */
    for (; 0 != *sp; ++sp) {}
    /* Free the strings in reverse order. */
    for (; string_array != sp; --sp) {
      free(*(sp - 1)); *(sp - 1) = 0;
    }
    /* Finally, free the pointer to the strings. */
    free(string_array);
  }
}

int ape_util_getenv(const char * name, char ** value, const char * def_value) {
  int status = eOK;

  /* Check arguments. */
  if (0 == name) {
    status = eNullPointer;
    ape_msg_debug("ape_util_getenv was called with name == 0.\n");
  } else if (0 == *name) {
    status = eInvalidArgument;
    ape_msg_debug("ape_util_getenv was called with name == \"\".\n");
  }
  if (0 == value) {
    status = eNullPointer;
    ape_msg_debug("ape_util_getenv was called with value == 0.\n");
  }
  if (eOK == status) {
    const char * tmp_value = getenv(name);
    if (0 != tmp_value) {
      *value = copy_string(tmp_value);
    } else {
      status = eVarNotSet;
      if (0 != def_value) {
        *value = copy_string(def_value);
      } else {
        *value = 0;
      }
    }
  }
  return status;
}

int ape_util_interpret_env(void) {
  char * value = 0;
  int status = ape_util_getenv("APEDEBUG", &value, 0);
  if (eOK == status) {
    char debug_mode = ape_msg_get_debug_mode();
    if (0 == debug_mode) {
      ape_msg_debug_enable(1);
      ape_msg_debug("APEDEBUG environment variable set: Ape debugging mode enabled.\n");
    }
  } else if (eVarNotSet == status) {
    /* That is OK, this is just an optional environment variable. */
    status = eOK;
  }
  free(value); value = 0;
  return status;
}

int ape_util_expand_env_var(const char * input, char ** output) {
  int status = eOK;

  /* Check inputs. */
  if (0 != output) *output = 0;
  else status = eNullPointer;
  if (eOK == status && 0 == input) status = eNullPointer;

  if (eOK == status) {
    static const char * begin_var[] = { "$ENV{", "$ENV(", "${", "$(" };
    static const char * end_var[] = { "}", ")", "}", ")" };
    static const size_t num_delim = sizeof(begin_var) / sizeof(begin_var[0]);
    const char * s_begin = input;
    const char * s_end = s_begin;
    char in_quote = 0;
    char * result = 0;
    while (eOK == status && '\0' != *s_end) {
      size_t delim_idx = 0;
      size_t begin_length = 0;

      /* Check for escaped characters. */
      if ('\\' == *s_end) {
        /* Skip this character and the next one. */
        ++s_end;
        if ('\0' != *s_end) ++s_end;
        continue;
      }

      /* Check for quoted characters. */
      if ('\'' == *s_end) {
        /* Toggle quote state. */
        in_quote = 0 == in_quote ? 1 : 0;
        ++s_end;
        continue;
      }

      /* Don't expand environment variable names inside single quotes. Note double quotes are OK. */
      if (0 != in_quote) {
        ++s_end;
        continue;
      }

      /* Find delimiter which indicates beginning of an environment variable name. */
      for (delim_idx = 0; delim_idx != num_delim; ++delim_idx) {
        if (0 == strncmp(s_end, begin_var[delim_idx], strlen(begin_var[delim_idx]))) {
          begin_length = strlen(begin_var[delim_idx]);
          break;
        }
      }

      if (0 != begin_length) {
        const char * n_begin = s_end + begin_length;
        const char * n_end = n_begin;
        char * name = 0;

        /* Any text which preceded this name should be concatenated to the output "as is". */
        if (s_begin < s_end) {
          result = *output;
          status = cat_range(result, s_begin, s_end, output);
          if (eOK == status) free(result);
          result = 0;
        }

        /* Now handle the named environment variable. */

        /* Look for the end of the variable name. */
        while ('\0' != *n_end && 0 != strncmp(n_end, end_var[delim_idx], strlen(end_var[delim_idx]))) ++n_end;

        if ('\0' == *n_end) {
          /* Reached end of string without finding termination. Do not expand variable, just write out as is. */
          result = *output;
          status = cat_range(result, s_begin, n_end, output);
          if (eOK == status) free(result);
          result = 0;
        } else if (n_begin < n_end) {
          /* Make a string containing just the name of the variable. */
          name = copy_range(n_begin, n_end);
          ++n_end;
          if (0 != name) {
            /* Get the value of the environment variable from the operating system. */
            const char * value = getenv(name);
            if (0 != value) {
              /* Concatenate the value to the output. */
              result = *output;
              status = ape_util_cat_string(result, value, output);
              if (eOK == status) free(result);
              result = 0;
            }
          } else {
            status = eDynAllocFailed;
          }
          /* Clean up. */
          free(name); name = 0;
        }

        /* Continue working through input after the end of this variable name. */
        s_begin = n_end;
        s_end = s_begin;

        /* s_begin and s_end were incremented as needed above, so continue the loop rather than
           falling through and incrementing them again. */
        continue;
      }
      ++s_end;
    }
    /* Any text which remains should be concatenated to the output "as is". */
    result = *output;
    status = cat_range(result, s_begin, s_end, output);
    if (eOK == status) free(result);
    result = 0;
  }
  return status;
}

int ape_util_parse_pfiles(const char * pfiles, char ** loc_pfiles, char ** sys_pfiles) {
  int status = eOK;
  /* Check input for validity. */
  if (0 == pfiles) {
    status = eNullPointer;
    ape_msg_debug("ape_util_parse_pfiles was called with pfiles == 0.\n");
  }

  /* Both loc_pfiles and sys_pfiles are optional. If non-0, initialize them to 0. */
  if (0 != loc_pfiles) *loc_pfiles = 0;
  if (0 != sys_pfiles) *sys_pfiles = 0;

  /* Parse out pfiles. */
  if (eOK == status) {
    /* Find the divider, if it is present in pfiles at all. */
    char * divider = strstr(pfiles, s_pfiles_divider);
    if (0 != loc_pfiles) {
      /* If divider was found, local pfiles is the part before the divider, otherwise, it's the whole thing. */
      if (0 != divider) *loc_pfiles = copy_range(pfiles, divider);
      else *loc_pfiles = copy_string(pfiles);
    }
    if (0 != sys_pfiles) {
      /* If divider was found, system pfiles is the part after the divider, otherwise, it's the whole thing. */
      if (0 != divider) *sys_pfiles = copy_string(divider + strlen(s_pfiles_divider));
      else *sys_pfiles = copy_string(pfiles);
    }
  }
  return status;
}
/* Convert a string to a boolean expression if at all possible, including converting numeric expressions (0 false, non-0 true).
   However, in that case, flag the conversion by returning a non-0 status.
*/
int ape_util_s2b(const char * input, char * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != input) {
      /* First attempt to convert string properly to a bool in one of the allowed forms (e.g. true, t, no, n etc.) */
      status = string2bool(input, result);
      if (eOK != status) {
        /* Indirect conversion: attempt to convert to a double cleanly. */
        double tmp_result = 0.;
        int double_status = string2double(input, &tmp_result);
        if (eOK == double_status) {
          /* Value converted as a double; interpret as false if 0 and true otherwise. */
          *result = (DBL_EPSILON < fabs(tmp_result)) ? 1 : 0;
          status = eTypeMismatch;
        } else {
          /* Converting as a double failed; let stand the result of the call to string2bool. */
        }
      } else {
        /* Conversion succeeded! */
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

/* Convert a string to a double expression.
*/
int ape_util_s2d(const char * input, double * result) {
  return string2double(input, result);
}

/* Convert a string to a float expression.
*/
int ape_util_s2f(const char * input, float * result) {
  int status = eOK;
  double double_result = 0.;
  if (0 == input || 0 == result) status = eNullPointer;
  if (eOK == status) {
    status = string2double(input, &double_result);
    *result = (float) double_result;
  }
  if (eOK == status) {
    if (FLT_MAX < fabs(double_result)) {
      status = eOverflow;
    }
  }

  return status;
}

/* Convert a string to a int expression if at all possible, including converting boolean or floating point expressions.
   However, in those cases, flag the conversion by returning a non-0 status.
*/
int ape_util_s2i(const char * input, int * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != input) {
      /* First attempt to convert to a double, in order to detect overflows and underflows more accurately. */
      double double_result = 0.;
      int double_status = string2double(input, &double_result);
      if (eNan == double_status) {
        *result = INT_MAX;
        status = double_status;
      } else if (eUndefinedValue == double_status) {
        status = double_status;
      } else if (INT_MAX < double_result) {
        *result = INT_MAX;
        status = eOverflow;
      } else if (INT_MIN > double_result) {
        *result = INT_MIN;
        status = eUnderflow;
      } else if (eUnderflow == double_status) {
        /* Floating point underflow means abs value is not measurably different from 0. This is considered
           a type mismatch, because the most likely case is e.g. 1.e-99999. */
        *result = 0;
        status = eTypeMismatch;
      } else if (eOverflow == double_status) {
        *result = INT_MAX;
        status = eOverflow;
      } else {
        long long_result = 0l;
        status = string2long(input, &long_result);
        if (eOK != status && eOK == double_status) {
          /* Value was cleanly converted to a double, but not to a long. Already handled above cases where
             double_result cannot be expressed by a long, so it is safe to use the double result, but flag
             this as a type mismatch. */
          *result = (int) double_result;
          status = eTypeMismatch;
        } else {
          *result = long_result;
        }
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

/* Convert a string to a long expression if at all possible, including converting boolean or floating point expressions.
   However, in those cases, flag the conversion by returning a non-0 status.
*/
int ape_util_s2l(const char * input, long * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0l;
    if (0 != input) {
      /* First attempt to convert to a double, in order to detect overflows and underflows more accurately. */
      double double_result = 0.;
      int double_status = string2double(input, &double_result);
      if (eNan == double_status) {
        *result = LONG_MAX;
        status = double_status;
      } else if (eUndefinedValue == double_status) {
        status = double_status;
      } else if (LONG_MAX < double_result) {
        *result = LONG_MAX;
        status = eOverflow;
      } else if (LONG_MIN > double_result) {
        *result = LONG_MIN;
        status = eUnderflow;
      } else if (eUnderflow == double_status) {
        /* Floating point underflow means abs value is not measurably different from 0. This is considered
           a type mismatch, because the most likely case is e.g. 1.e-99999. */
        *result = 0l;
        status = eTypeMismatch;
      } else if (eOverflow == double_status) {
        *result = LONG_MAX;
        status = eOverflow;
      } else {
        status = string2long(input, result);
        if (eOK != status && eOK == double_status) {
          /* Value was cleanly converted to a double, but not to a long. Already handled above cases where
             double_result cannot be expressed by a long, so it is safe to use the double result, but flag
             this as a type mismatch. */
          *result = (long) double_result;
          status = eTypeMismatch;
        }
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

/* Convert a string to a int expression if at all possible, including converting boolean or floating point expressions.
   However, in those cases, flag the conversion by returning a non-0 status.
*/
int ape_util_s2sh(const char * input, short * result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != input) {
      /* First attempt to convert to a double, in order to detect overflows and underflows more accurately. */
      double double_result = 0.;
      int double_status = string2double(input, &double_result);
      if (eNan == double_status) {
        *result = SHRT_MAX;
        status = double_status;
      } else if (eUndefinedValue == double_status) {
        status = double_status;
      } else if (SHRT_MAX < double_result) {
        *result = SHRT_MAX;
        status = eOverflow;
      } else if (SHRT_MIN > double_result) {
        *result = SHRT_MIN;
        status = eUnderflow;
      } else if (eUnderflow == double_status) {
        /* Floating point underflow means abs value is not measurably different from 0. This is considered
           a type mismatch, because the most likely case is e.g. 1.e-99999. */
        *result = 0;
        status = eTypeMismatch;
      } else if (eOverflow == double_status) {
        *result = SHRT_MAX;
        status = eOverflow;
      } else {
        long long_result = 0l;
        status = string2long(input, &long_result);
        if (eOK != status && eOK == double_status) {
          /* Value was cleanly converted to a double, but not to a long. Already handled above cases where
             double_result cannot be expressed by a long, so it is safe to use the double result, but flag
             this as a type mismatch. */
          *result = (int) double_result;
          status = eTypeMismatch;
        } else {
          *result = long_result;
        }
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_util_strcmp(const char * s1, const char * s2, char case_insensitive) {
  int status = 0;
  if (s1 != s2) {
    if (0 == s1) {
      status = -1;
    } else if (0 == s2) {
      status = +1;
    } else if (0 != case_insensitive) {
      /* Scan through both strings until a case-insensitive discrepancy is found or the end of one string is reached. */
      for (; '\0' != *s1 && '\0' != *s2 && toupper(*s1) == toupper(*s2); ++s1, ++s2) {}
      status = toupper(*s1) - toupper(*s2);
    } else {
      /* For case-sensitive comparison, just use native strcmp. */
      status = strcmp(s1, s2);
    }
  }
  return status;
}

int ape_util_cmp_string_array(const char ** s1, const char ** s2, char case_insensitive) {
  int status = 0;
  if (s1 != s2) {
    if (0 == s1) {
      status = -1;
    } else if (0 == s2) {
      status = 1;
    } else {
      /* Compare each string in turn, stopping when either string null terminates or when
         a discrepancy is found. */
      for (; 0 == status && 0 != *s1 && 0 != *s2; ++s1, ++s2) {
        status = ape_util_strcmp(*s1, *s2, case_insensitive);
      }
      /* If status is 0, the loop above exited because one of the pointers is 0. Therefore,
         in this case, do one last check in case both pointers are not 0; ape_util_strcmp
         can handle that. */
      if (0 == status) status = ape_util_strcmp(*s1, *s2, case_insensitive);
    }
  }
  return status;
}

/* Utility to look for a string in a range of pointers to strings. */
int ape_util_find_string(char ** begin, char ** end, const char * input, char *** found, char case_insensitive) {
  int status = eOK;

  /* Check arguments. */
  if (0 == begin || 0 == end || 0 == input || 0 == found) status = eNullPointer;

  if (eOK == status) {
    /* Initialize output. */
    *found = 0;
  }

  if (eOK == status) {
    char ** itor = 0;
    /* Check whether input matches one of the strings in the range. */
    for (itor = begin; itor < end && 0 == *found; ++itor) {
      status = ape_util_strcmp(input, *itor, case_insensitive);
      if (eOK == status) *found = itor;
    }
  }

  return status;
}

static int check_file_access(const char * file_name, const char * open_mode) {
  int status = eOK;
  if (0 == file_name || 0 == open_mode) status = eNullPointer;
  if (eOK == status) {
#ifdef WIN32
    /* Check existence no matter what. */
    DWORD file_attributes = GetFileAttributes(file_name);
    char file_exists = (s_file_not_found != file_attributes) ? 1 : 0;

    /* Confirm file exists if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'E') || 0 != strchr(open_mode, 'e'))) {
      if (0 == file_exists) status = eFileNotAccessible;
    }

    /* Confirm file does not exist if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'N') || 0 != strchr(open_mode, 'n'))) {
      if (0 != file_exists) status = eFileNotAccessible;
    }

    /* Confirm file is readable if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'R') || 0 != strchr(open_mode, 'r'))) {
      /* On Windows, assume file exists and is readable are equivalent. TODO: Is there a better way? */
      if (0 == file_exists) status = eFileNotAccessible;
    }

    /* Confirm file is writable if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'W') || 0 != strchr(open_mode, 'w'))) {
      if (s_file_not_found == file_attributes || (0 != (FILE_ATTRIBUTE_READONLY & file_attributes))) status = eFileNotAccessible;

      if (eFileNotAccessible == status) {
        /* It might be that the file does not exist, but that it would be possible to create it. */
        FILE * fp = 0;
        fp = fopen(file_name, "w");
        if (0 != fp) {
          /* Succeeded in creating the file. */
          status = eOK;
          fclose(fp); fp = 0;
        }
        /* Clean up in the event the file did not exist prior to creating it. */
        if (0 == file_exists) remove(file_name);
      }
    }
#else
    /* Check existence no matter what. */
    char file_exists = (0 == access(file_name, F_OK)) ? 1 : 0;

    /* Confirm file exists if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'E') || 0 != strchr(open_mode, 'e'))) {
      if (0 == file_exists) status = eFileNotAccessible;
    }

    /* Confirm file does not exist if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'N') || 0 != strchr(open_mode, 'n'))) {
      if (0 != file_exists) status = eFileNotAccessible;
    }

    /* Confirm file is readable if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'R') || 0 != strchr(open_mode, 'r'))) {
      if (0 != access(file_name, R_OK)) status = eFileNotAccessible;
    }

    /* Confirm file is writable if requested. */
    if (eOK == status && (0 != strchr(open_mode, 'W') || 0 != strchr(open_mode, 'w'))) {
      if (0 != access(file_name, W_OK)) status = eFileNotAccessible;

      if (eFileNotAccessible == status) {
        /* It might be that the file does not exist, but that it would be possible to create it. */
        FILE * fp = 0;
        fp = fopen(file_name, "w");
        if (0 != fp) {
          /* Succeeded in creating the file. */
          status = eOK;
          fclose(fp); fp = 0;
        }
        /* Clean up in the event the file did not exist prior to creating it. */
        if (0 == file_exists) remove(file_name);
      }
    }
#endif
  }

  return status;
}

static int (*s_check_file_func)(const char *, const char *) = &check_file_access;

int ape_util_check_file_access(const char * file_name, const char * access) {
  int status = eOK;
  if (0 == file_name || 0 == access) status = eNullPointer;
/*  if (eOK == status && '\0' == *access) return eOK; */
  if (eOK == status) status = (*s_check_file_func)(file_name, access);
  return status;
}

int ape_util_get_file_check_func(int (**func)(const char *, const char *)) {
  int status = eOK;
  if (0 == func) status = eNullPointer;

  /* Get the pointer to the function currently being used to check file access. */
  if (eOK == status) *func = s_check_file_func;

  return eOK;
}

int ape_util_set_file_check_func(int (*func)(const char *, const char *)) {
  if (0 == func) s_check_file_func = check_file_access;
  else s_check_file_func = func;
  return eOK;
}

void ape_util_sleep(int sleep_time) {
#ifdef WIN32
  Sleep(sleep_time * 1000);
#else
  sleep(sleep_time);
#endif
}

#ifdef USE_READLINE
/* Clean up after readline's allocations in the hopes of reducing memory sloppiness. */
static void free_readline_memory(void) {
  /* Clean up readline's history list. */
  clear_history();

#if 0
  /* This is not compiled currently because it is not known whether funmap and rl_line_buffer
     are always available. However in the hopes one day of not leaking any memory, keeping the
     code here for possible future use. */
  /* The type of the pointer is not correct, but the ape call does the right thing and char * will work so just cast it. */
  ape_util_free_string_array((char **) funmap);

  /* Free up readline's main buffer. */
  free(rl_line_buffer);
#endif

}
#endif

#ifndef USE_READLINE
#define APE_BUF_SIZE 8192

/* Note: it would be nice if the following function worked in the readline case, but it seems that
   mixing readline calls with direct reads from stdin causes subsequent calls to readline to break.
   This is why ape_par_redirect_prompt_stream does not accept NULL prompt streams. */
static char * read_from_stdin(void) {
  char * text = 0;
  char buf[APE_BUF_SIZE] = "";
  size_t num_read = 0;
  const size_t buf_size = sizeof(buf)/sizeof(buf[0]);

  /* Start with a blank buffer. */
  memset(buf, '\0', buf_size);

  /* Get the input. */
  if (0 != fgets(buf, (int) buf_size, stdin)) {

    /* Determine the number of characters which were actually read. */
    num_read = strlen(buf);

    /* Blank out the newline, if any. */
    if ('\n' == buf[num_read - 1]) buf[num_read - 1] = '\0';

    /* Copy the buffer to the output string. */
    ape_util_copy_string(buf, &text);
  }

  return text;
}

#undef APE_BUF_SIZE
#endif

static char * get_text(const char * prompt) {
  char * text = 0;

#ifndef USE_READLINE
  FILE * prompt_stream = *get_prompt_stream_ptr();
  if (0 != prompt_stream) {

    /* Issue the prompt. */
    fprintf(prompt_stream, "%s", prompt);

    /* Flush buffer so that prompt appears right away even if stream is buffered. */
    fflush(prompt_stream);
  }

  text = read_from_stdin();

#else
  /* Clean up after readline's global allocations. */
  ape_util_atexit(&free_readline_memory);

  /* Use readline to get the user's input. */
  /* Note: it would be nice if the following block worked in the readline case, but it seems that
     mixing readline calls with direct reads from stdin causes subsequent calls to readline to break.
     This is why ape_par_redirect_prompt_stream does not accept NULL prompt streams. */
  /* if (0 != prompt_stream) {
       text = readline(prompt);
     } else {
       text = read_from_stdin();
     }
  */
  text = readline(prompt);

  /* Store significant input in the history list. */
  if (0 != text && '\0' != *text) {
    add_history(text);
  }

#endif

  return text;
}

int ape_util_get_text(const char * prompt, char ** text) {
  int status = eOK;
  if (0 != text && 0 != prompt) {
    *text = get_text(prompt);
    if (0 == *text) status = eInputFailure;
  } else {
    status = eNullPointer;
  }
  return status;
}

static FILE ** get_prompt_stream_ptr(void) {
  static FILE * s_prompt_stream = 0;
  static char s_stream_init = 0;
  if (0 == s_stream_init) {
    s_stream_init = 1;
    s_prompt_stream = stdout;
  }
  return &s_prompt_stream;
}

int ape_util_set_prompt_stream(FILE * stream) {
  int status = eOK;

  /* See below for why this is not allowed. */
  if (0 == stream) {
    status = eNullPointer;
    ape_msg_debug("ape_util_set_prompt_stream was called with null stream.\n");
  }

  if (eOK == status) {
#ifdef USE_READLINE
    /* Note: if rl_outstream is NULL, calling readline() produces a segmentation fault, thus we don't allow it ever. */
    rl_outstream = stream;
#endif

    *get_prompt_stream_ptr() = stream;
  }

  return status;
}

static void check_s2b(const char * input, int expected_output, int expected_status, int output, int status) {
  if (expected_status != status)
    ape_test_failed("ape_util_s2b(\"%s\", &bool_result) returned status %d, not %d as expected.\n",
      0 != input ? input : "0", status, expected_status);
  if (expected_output != output)
    ape_test_failed("ape_util_s2b(\"%s\", &bool_result) converted to %d, not %d as expected.\n",
      0 != input ? input : "0", output, expected_output);
}

static void check_s2d(const char * input, double expected_output, int expected_status, double output, int status) {
  if (expected_status != status)
    ape_test_failed("ape_util_s2d(\"%s\", &double_result) returned status %d, not %d as expected.\n",
      0 != input ? input : "0", status, expected_status);
  if (expected_output != output)
    ape_test_failed("ape_util_s2d(\"%s\", &double_result) converted to %g, not %g as expected.\n",
      0 != input ? input : "0", output, expected_output);
}

static void check_s2i(const char * input, int expected_output, int expected_status, int output, int status) {
  if (expected_status != status)
    ape_test_failed("ape_util_s2l(\"%s\", &int_result) returned status %d, not %d as expected.\n",
      0 != input ? input : "0", status, expected_status);
  if (expected_output != output)
    ape_test_failed("ape_util_s2l(\"%s\", &int_result) converted to %ld, not %ld as expected.\n",
      0 != input ? input : "0", output, expected_output);
}

static void check_s2l(const char * input, long expected_output, int expected_status, long output, int status) {
  if (expected_status != status)
    ape_test_failed("ape_util_s2l(\"%s\", &long_result) returned status %d, not %d as expected.\n",
      0 != input ? input : "0", status, expected_status);
  if (expected_output != output)
    ape_test_failed("ape_util_s2l(\"%s\", &long_result) converted to %ld, not %ld as expected.\n",
      0 != input ? input : "0", output, expected_output);
}

static void check_s2sh(const char * input, short expected_output, int expected_status, short output, int status) {
  if (expected_status != status)
    ape_test_failed("ape_util_s2sh(\"%s\", &short_result) returned status %d, not %d as expected.\n",
      0 != input ? input : "0", status, expected_status);
  if (expected_output != output)
    ape_test_failed("ape_util_s2sh(\"%s\", &short_result) converted to %ld, not %ld as expected.\n",
      0 != input ? input : "0", output, expected_output);
}

static int dont_check_access(const char * file_name, const char * access) { return eOK; }

static void test_cmp_check_ptr(const char * hint, int (* result)(const char *, const char *),
  int (* expected_result)(const char *, const char *), int status, int expected_status) {
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if (result != expected_result) {
    ape_test_failed("%s%spointer to file access checking function did not have expected value.\n", hint, join);
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

/* Test conversions to/from strings from/to numbers. */
static void test_conversion(void) {
  const char * input = 0;
  const char * input2 = 0;
  char input3[] = "fobble";
  int expected_int = 0;
  long expected_long = 0l;
  short expected_short = 0;
  int expected_status = 0;
  int int_result = 0;
  long long_result = 0l;
  short short_result = 0;
  char * expected_string = 0;
  char * string_result = 0;
  char expected_bool = 0;
  char bool_result = 0;
  double expected_double = 0.;
  double double_result = 0.;
  int status = 0;

  /* Test ci_strcmp. */
  /* Test results when one or more null pointers are passed. */
  expected_status = 0;
  input = 0;
  status = ci_strcmp(input, input);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(0, 0) returned %d, not %d as expected.\n", status, expected_status);

  expected_status = 1;
  input = "fobble";
  input2 = 0;
  status = ci_strcmp(input, input2);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", 0) returned %d, not %d as expected.\n", input, status, expected_status);

  expected_status = -1;
  input = 0;
  input2 = "fobble";
  status = ci_strcmp(input, input2);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(0, \"%s\") returned %d, not %d as expected.\n", input, status, expected_status);

  /* Test comparisons involving 2 strings. */
  expected_status = 0;
  input = "fobble";
  status = ci_strcmp(input, input);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", \"%s\") returned %d, not %d as expected.\n", input, input, status, expected_status);

  expected_status = 0;
  input = "fobble";
  /* Note: input3 = "fobble", but a separate copy so that the pointers will not be equal. */
  status = ci_strcmp(input, input3);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", \"%s\") returned %d, not %d as expected.\n", input, input3, status, expected_status);

  expected_status = 0;
  input = "fobble";
  input2 = "FoBblE";
  status = ci_strcmp(input, input2);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", \"%s\") returned %d, not %d as expected.\n", input, input2, status, expected_status);

  expected_status = -1 * ' ';
  input = "fobble";
  input2 = "fobble ";
  status = ci_strcmp(input, input2);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", \"%s\") returned %d, not %d as expected.\n", input, input2, status, expected_status);

  expected_status = ' ';
  input = "fobble ";
  input2 = "fobble";
  status = ci_strcmp(input, input2);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", \"%s\") returned %d, not %d as expected.\n", input, input2, status, expected_status);

  expected_status = 1;
  input = "fobble";
  input2 = "eobble";
  status = ci_strcmp(input, input2);
  if (expected_status != status)
    ape_test_failed("ci_strcmp(\"%s\", \"%s\") returned %d, not %d as expected.\n", input, input2, status, expected_status);

  /* Test truncate_white_space. */
  expected_string = 0;
  input = 0;
  string_result = truncate_white_space(input);
  if (string_result != expected_string)
    ape_test_failed("truncate_white_space(0) returned \"%s\", not 0 as expected.\n", string_result);
  free(string_result); string_result = 0;

  expected_string = "fobble";
  input = " \tfobble \v\t ";
  string_result = truncate_white_space(input);
  if (0 != strcmp(string_result, expected_string))
    ape_test_failed("truncate_white_space(\"s\") returned \"%s\", not \"%s\" as expected.\n", input, string_result,
      expected_string);
  free(string_result); string_result = 0;

  expected_string = "";
  input = " \t \v\t ";
  string_result = truncate_white_space(input);
  if (0 != strcmp(string_result, expected_string))
    ape_test_failed("truncate_white_space(\"s\") returned \"%s\", not \"%s\" as expected.\n", input, string_result,
      expected_string);
  free(string_result); string_result = 0;

  /* Test ape_util_s2i. */
  /* Conversions involving one or more null pointers. */
  input = 0;
  expected_int = 0;
  expected_status = eNullPointer;
  int_result = 0;
  status = ape_util_s2i(input, 0);
  check_s2i(input, expected_int, expected_status, int_result, status);

  input = 0;
  expected_int = 0;
  expected_status = eNullPointer;
  int_result = 45;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  input = " 45\t";
  expected_int = 45;
  expected_status = eNullPointer;
  int_result = 45;
  status = ape_util_s2i(input, 0);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which fails due to trailing characters. */
  input = "45 fish";
  expected_int = 45;
  expected_status = eStringRemainder;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which fails due to mismatched type. */
  input = "45.";
  expected_int = 45;
  expected_status = eTypeMismatch;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which fails due to mismatched type, but also has trailing characters. There is no way to distinguish
     this from a purely trailing character failure.*/
  input = "45. fish";
  expected_int = 45;
  expected_status = eStringRemainder;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which fails due to overflow. */
  input = "45000000000000000000000000000000000000000000000";
  expected_int = INT_MAX;
  expected_status = eOverflow;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which fails due to underflow. */
  input = "-45000000000000000000000000000000000000000000000";
  expected_int = INT_MIN;
  expected_status = eUnderflow;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which fails due to just plain non numeric value. */
  input = "the quick brown fox jumps over the lazy dog";
  expected_int = 0;
  expected_status = eStringRemainder;
  int_result = 45;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which overflows due to a converted NAN. */
  input = "NAN";
  expected_int = INT_MAX;
  expected_status = eNan;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Straightforward conversion which succeeds. */
  input = "45";
  expected_int = 45;
  expected_status = eOK;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which succeeds despite leading and trailing white space, and hexadecimal representation. */
  input = " 0x2d\t";
  expected_int = 45;
  expected_status = eOK;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Conversion which succeeds despite leading white space, and exponential notation (with error about type conversion). */
  input = " 1.024e3\t";
  expected_int = 1024;
  expected_status = eTypeMismatch;
  int_result = 0;
  status = ape_util_s2i(input, &int_result);
  check_s2i(input, expected_int, expected_status, int_result, status);

  /* Test ape_util_s2l. */
  /* Conversions involving one or more null pointers. */
  input = 0;
  expected_long = 0l;
  expected_status = eNullPointer;
  long_result = 0l;
  status = ape_util_s2l(input, 0);
  check_s2l(input, expected_long, expected_status, long_result, status);

  input = 0;
  expected_long = 0l;
  expected_status = eNullPointer;
  long_result = 45l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  input = " 45\t";
  expected_long = 45l;
  expected_status = eNullPointer;
  long_result = 45l;
  status = ape_util_s2l(input, 0);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which fails due to trailing characters. */
  input = "45 fish";
  expected_long = 45l;
  expected_status = eStringRemainder;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which fails due to mismatched type. */
  input = "45.";
  expected_long = 45l;
  expected_status = eTypeMismatch;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which fails due to mismatched type, but also has trailing characters. There is no way to distinguish
     this from a purely trailing character failure.*/
  input = "45. fish";
  expected_long = 45l;
  expected_status = eStringRemainder;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which fails due to overflow. */
  input = "45000000000000000000000000000000000000000000000";
  expected_long = LONG_MAX;
  expected_status = eOverflow;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which fails due to underflow. */
  input = "-45000000000000000000000000000000000000000000000";
  expected_long = LONG_MIN;
  expected_status = eUnderflow;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which fails due to just plain non numeric value. */
  input = "the quick brown fox jumps over the lazy dog";
  expected_long = 0l;
  expected_status = eStringRemainder;
  long_result = 45l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which overflows due to a converted NAN. */
  input = "NAN";
  expected_long = LONG_MAX;
  expected_status = eNan;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Straightforward conversion which succeeds. */
  input = "45";
  expected_long = 45l;
  expected_status = eOK;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which succeeds despite leading and trailing white space, and hexadecimal representation. */
  input = " 0x2d\t";
  expected_long = 45l;
  expected_status = eOK;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Conversion which succeeds despite leading white space, and exponential notation (with error about type conversion). */
  input = " 1.024e3\t";
  expected_long = 1024l;
  expected_status = eTypeMismatch;
  long_result = 0l;
  status = ape_util_s2l(input, &long_result);
  check_s2l(input, expected_long, expected_status, long_result, status);

  /* Test ape_util_s2sh. */
  /* Conversions involving one or more null pointers. */
  input = 0;
  expected_short = 0;
  expected_status = eNullPointer;
  short_result = 0;
  status = ape_util_s2sh(input, 0);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  input = 0;
  expected_short = 0;
  expected_status = eNullPointer;
  short_result = 45;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  input = " 45\t";
  expected_short = 45;
  expected_status = eNullPointer;
  short_result = 45;
  status = ape_util_s2sh(input, 0);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which fails due to trailing characters. */
  input = "45 fish";
  expected_short = 45;
  expected_status = eStringRemainder;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which fails due to mismatched type. */
  input = "45.";
  expected_short = 45;
  expected_status = eTypeMismatch;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which fails due to mismatched type, but also has trailing characters. There is no way to distinguish
     this from a purely trailing character failure.*/
  input = "45. fish";
  expected_short = 45;
  expected_status = eStringRemainder;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which fails due to overflow. */
  input = "45000000000000000000000000000000000000000000000";
  expected_short = SHRT_MAX;
  expected_status = eOverflow;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which fails due to underflow. */
  input = "-45000000000000000000000000000000000000000000000";
  expected_short = SHRT_MIN;
  expected_status = eUnderflow;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which fails due to just plain non numeric value. */
  input = "the quick brown fox jumps over the lazy dog";
  expected_short = 0;
  expected_status = eStringRemainder;
  short_result = 45;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which overflows due to a converted NAN. */
  input = "NAN";
  expected_short = SHRT_MAX;
  expected_status = eNan;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Straightforward conversion which succeeds. */
  input = "45";
  expected_short = 45;
  expected_status = eOK;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which succeeds despite leading and trailing white space, and hexadecimal representation. */
  input = " 0x2d\t";
  expected_short = 45;
  expected_status = eOK;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Conversion which succeeds despite leading white space, and exponential notation (with error about type conversion). */
  input = " 1.024e3\t";
  expected_short = 1024;
  expected_status = eTypeMismatch;
  short_result = 0;
  status = ape_util_s2sh(input, &short_result);
  check_s2sh(input, expected_short, expected_status, short_result, status);

  /* Test ape_util_s2b. */
  /* Conversions involving one or more null pointers. */
  input = 0;
  expected_bool = 0;
  expected_status = eNullPointer;
  bool_result = 0;
  status = ape_util_s2b(input, 0);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = 0;
  expected_bool = 0;
  expected_status = eNullPointer;
  bool_result = 45;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = " truE\t";
  expected_bool = 1;
  expected_status = eNullPointer;
  bool_result = 1;
  status = ape_util_s2b(input, 0);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  /* Conversions which fail due to trailing characters. */
  input = " \ttruehood";
  expected_bool = 0;
  expected_status = eConversionError;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = "ye";
  expected_bool = 0;
  expected_status = eConversionError;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  /* Conversion which throws error but converts correctly due to numeric input. */
  input = " \t123. ";
  expected_bool = 1;
  expected_status = eTypeMismatch;
  bool_result = 123;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  /* Conversion which throws error but converts correctly due to numeric input. */
  input = "0";
  expected_bool = 0;
  expected_status = eTypeMismatch;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  /* Conversion which throws error but converts correctly due to numeric input. */
  input = "1";
  expected_bool = 1;
  expected_status = eTypeMismatch;
  bool_result = 0;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  /* Conversions which succeed. */
  input = " \tF ";
  expected_bool = 0;
  expected_status = eOK;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = " \tFalSe ";
  expected_bool = 0;
  expected_status = eOK;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = " \vnO";
  expected_bool = 0;
  expected_status = eOK;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = "n";
  expected_bool = 0;
  expected_status = eOK;
  bool_result = 1;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = " t ";
  expected_bool = 1;
  expected_status = eOK;
  bool_result = 0;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = "true";
  expected_bool = 1;
  expected_status = eOK;
  bool_result = 0;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = "yEs";
  expected_bool = 1;
  expected_status = eOK;
  bool_result = 0;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  input = "y";
  expected_bool = 1;
  expected_status = eOK;
  bool_result = 0;
  status = ape_util_s2b(input, &bool_result);
  check_s2b(input, expected_bool, expected_status, bool_result, status);

  /* Test ape_util_s2d. */
  /* Conversions involving one or more null pointers. */
  input = 0;
  expected_double = 0.;
  expected_status = eNullPointer;
  double_result = 0.;
  status = ape_util_s2d(input, 0);
  check_s2d(input, expected_double, expected_status, double_result, status);

  input = "1.e-3 \t\v";
  expected_double = 1.e-3;
  expected_status = eNullPointer;
  double_result = 1.e-3;
  status = ape_util_s2d(input, 0);
  check_s2d(input, expected_double, expected_status, double_result, status);

  input = 0;
  expected_double = 0.;
  expected_status = eNullPointer;
  double_result = 1.e-3;
  status = ape_util_s2d(input, &double_result);
  check_s2d(input, expected_double, expected_status, double_result, status);

  /* Conversion which fails due to trailing characters. */
  input = "1.e-3e7";
  expected_double = 1.e-3;
  expected_status = eStringRemainder;
  double_result = 1.e3;
  status = ape_util_s2d(input, &double_result);
  check_s2d(input, expected_double, expected_status, double_result, status);

  /* Conversions which test overflow/underflow. */
  input = "1.e99999";
  string2double("1.e99999", &expected_double);
  expected_status = eOverflow;
  double_result = 1.e-3;
  status = ape_util_s2d(input, &double_result);
  check_s2d(input, expected_double, expected_status, double_result, status);

  input = "1.e-99999";
  expected_double = 0.;
  expected_status = eUnderflow;
  double_result = 1.e30;
  status = ape_util_s2d(input, &double_result);
  check_s2d(input, expected_double, expected_status, double_result, status);

  /* Check behavior of ape_util_getenv. */
  { const char * name = 0;
    char ** value = 0;
    const char * def_value = 0;
    status = ape_util_getenv(name, value, def_value);
    if (eNullPointer != status) {
      ape_test_failed("ape_util_getenv(0, 0, 0) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  { const char * name = 0;
    char * value = 0;
    const char * def_value = 0;
    status = ape_util_getenv(name, &value, def_value);
    if (eNullPointer != status) {
      ape_test_failed("ape_util_getenv(0, &value, 0) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
    if (0 != value) {
      ape_test_failed("ape_util_getenv(0, &value, 0) returned value %s, not 0 as expected.\n", value);
    }
  }
  { const char * name = "";
    char * value = 0;
    const char * def_value = 0;
    status = ape_util_getenv(name, &value, def_value);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, 0) returned status %d, not %d as expected.\n", name, status,
        eInvalidArgument);
    }
  }
  { const char * name = "THIS_ENV_VAR_SHOULD_NOT_BE_SET";
    char ** value = 0;
    const char * def_value = 0;
    status = ape_util_getenv(name, value, def_value);
    if (eNullPointer != status) {
      ape_test_failed("ape_util_getenv(\"%s\", 0, 0) returned status %d, not %d as expected.\n", name, status, eNullPointer);
    }
  }
  { const char * name = "THIS_ENV_VAR_SHOULD_NOT_BE_SET";
    char * value = 0;
    const char * def_value = 0;
    status = ape_util_getenv(name, &value, def_value);
    if (eVarNotSet != status) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, 0) returned status %d, not %d as expected.\n", name, status, eVarNotSet);
    }
    if (0 != value) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, 0) returned value \"%s\", not 0 as expected.\n", name, value);
    }
  }
  { const char * name = "THIS_ENV_VAR_SHOULD_NOT_BE_SET";
    char * value = 0;
    const char * def_value = "non-empty default value";
    status = ape_util_getenv(name, &value, def_value);
    if (eVarNotSet != status) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, \"%s\") returned status %d, not %d as expected.\n",
        name, def_value, status, eVarNotSet);
    }
    if (0 == value) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, \"%s\") unexpectedly returned value == 0.\n", name, def_value);
    } else if (0 != strcmp(value, def_value)) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, \"%s\") returned value \"%s\", not the default value as expected.\n",
        name, def_value, value);
    }
    free(value); value = 0;
  }
  /* Finally, some success cases. */
  { const char * name = "PATH";
    char * value = 0;
    const char * def_value = 0;
    status = ape_util_getenv(name, &value, def_value);
    if (eOK != status) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, 0) returned status %d, not %d as expected.\n",
        name, status, eOK);
    }
    if (0 == value) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, 0) unexpectedly returned value == 0.\n", name);
    }
    free(value); value = 0;
  }
  { const char * name = "PATH";
    char * value = 0;
    const char * def_value = "PATH variable is not set";
    status = ape_util_getenv(name, &value, def_value);
    if (eOK != status) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, \"%s\") returned status %d, not %d as expected.\n",
        name, def_value, status, eVarNotSet);
    }
    if (0 == value) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, \"%s\") unexpectedly returned value == 0.\n", name, def_value);
    } else if (0 == strcmp(value, def_value)) {
      ape_test_failed("ape_util_getenv(\"%s\", &value, \"%s\") unexpected returned default value \"%s\".\n",
        name, def_value, value);
    }
    free(value); value = 0;
  }

  /* Test ape_util_parse_pfiles. */
  /* Test error cases. */
  { const char * pfiles = 0;
    char * loc_pfiles = 0;
    char * sys_pfiles = 0;
    status = ape_util_parse_pfiles(pfiles, &loc_pfiles, &sys_pfiles);
    if (eNullPointer != status) {
      ape_test_failed("ape_util_parse_pfiles(0, 0, 0) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  /* Test case in which nothing is actually parsed because both outputs are 0. */
  { const char * pfiles = ".";
    char ** loc_pfiles = 0;
    char ** sys_pfiles = 0;
    /* Note this has no effect, but it is legal. */
    status = ape_util_parse_pfiles(pfiles, loc_pfiles, sys_pfiles);
    if (eOK != status) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", 0, 0) returned status %d, not %d as expected.\n",
        pfiles, status, eNullPointer);
    }
  }
  /* Get local part only. */
  { const char * pfiles = ".";
    char * loc_pfiles = 0;
    char ** sys_pfiles = 0;
    status = ape_util_parse_pfiles(pfiles, &loc_pfiles, sys_pfiles);
    if (eOK != status) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, 0) returned status %d, not %d as expected.\n",
        pfiles, status, eNullPointer);
    }
    if (0 == loc_pfiles) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, 0) unexpectedly returned loc_pfiles == 0.\n",
        pfiles);
    } else if (0 != strcmp(loc_pfiles, ".")) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, 0) returned loc_pfiles == \"%s\", not \".\", as expected.\n",
        pfiles, loc_pfiles);
    }
    free(loc_pfiles); loc_pfiles = 0;
  }
  /* Get system part only. */
  { const char * pfiles = ".";
    char ** loc_pfiles = 0;
    char * sys_pfiles = 0;
    status = ape_util_parse_pfiles(pfiles, loc_pfiles, &sys_pfiles);
    if (eOK != status) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", 0, &sys_pfiles) returned status %d, not %d as expected.\n",
        pfiles, status, eNullPointer);
    }
    if (0 == sys_pfiles) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", 0, &sys_pfiles) unexpectedly returned sys_pfiles == 0.\n",
        pfiles);
    } else if (0 != strcmp(sys_pfiles, ".")) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", 0, &sys_pfiles) returned sys_pfiles == \"%s\", not \".\", as expected.\n",
        pfiles, sys_pfiles);
    }
    free(sys_pfiles); sys_pfiles = 0;
  }
  /* Get both parts in a few non-trivial cases. */
  { const char * pfiles = "." PFILES_DIVIDER "..";
    char * loc_pfiles = 0;
    char * sys_pfiles = 0;
    status = ape_util_parse_pfiles(pfiles, &loc_pfiles, &sys_pfiles);
    if (eOK != status) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", 0, &sys_pfiles) returned status %d, not %d as expected.\n",
        pfiles, status, eNullPointer);
    }
    if (0 == loc_pfiles) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, &sys_pfiles) unexpectedly returned loc_pfiles == 0.\n",
        pfiles);
    } else if (0 != strcmp(loc_pfiles, ".")) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, &sys_pfiles) returned loc_pfiles == \"%s\", not "
        "\"..\", as expected.\n", pfiles, loc_pfiles);
    }
    if (0 == sys_pfiles) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, &sys_pfiles) unexpectedly returned sys_pfiles == 0.\n",
        pfiles);
    } else if (0 != strcmp(sys_pfiles, "..")) {
      ape_test_failed("ape_util_parse_pfiles(\"%s\", &loc_pfiles, &sys_pfiles) returned sys_pfiles == \"%s\", not "
        "\"..\", as expected.\n", pfiles, sys_pfiles);
    }
    free(sys_pfiles); sys_pfiles = 0;
    free(loc_pfiles); sys_pfiles = 0;
  }

  /* Test ape_util_cat_string. */
  { const char * s1 = "can";
    const char * s2 = "dle";
    const char * expected = "candle";
    char * result = 0;
    status = ape_util_cat_string(s1, s2, &result);
    ape_test_cmp_string("ape_util_cat_string(s1, s2, &result)", result, expected, status, eOK);
    free(result);
  }
  { const char * s1 = "can";
    const char * s2 = 0;
    const char * expected = s1;
    char * result = 0;
    status = ape_util_cat_string(s1, s2, &result);
    ape_test_cmp_string("ape_util_cat_string(s1, 0, &result)", result, expected, status, eOK);
    free(result);
  }
  { const char * s1 = 0;
    const char * s2 = "dle";
    const char * expected = s2;
    char * result = 0;
    status = ape_util_cat_string(s1, s2, &result);
    ape_test_cmp_string("ape_util_cat_string(0, s2, &result)", result, expected, status, eOK);
    free(result);
  }
  { const char * s1 = "can";
    const char * s2 = "dle";
    const char * expected = 0;
    char ** result = 0;
    status = ape_util_cat_string(s1, s2, result);
    ape_test_cmp_string("ape_util_cat_string(s1, s2, 0)", expected, expected, status, eNullPointer);
  }

  /* Test ape_util_strcmp. */
  { const char * string1 = "ab";
    const char * string2 = "Ae";
    status = ape_util_strcmp(0, 0, 0);
    if (0 != status) ape_test_failed("ape_util_strcmp(0, 0, 0) returned %d, not 0 as expected.\n", status);

    status = ape_util_strcmp(string1, 0, 1);
    if (-1 == status) ape_test_failed("ape_util_strcmp(string1, 0, 1) returned %d, not 0 as expected.\n", status);

    status = ape_util_strcmp(0, string2, 1);
    if (1 == status) ape_test_failed("ape_util_strcmp(0, string2, 1) returned %d, not 0 as expected.\n", status);

    /* Case insensitive test should return same as native strcmp. */
    status = ape_util_strcmp(string1, string2, 0);
    expected_status = strcmp(string1, string2);
    if (expected_status != status)
      ape_test_failed("ape_util_strcmp(string1, string2, 0) returned %d, not %d as expected.\n", status, expected_status);

    status = ape_util_strcmp(string1, string2, 1);
    if (0 == status) ape_test_failed("ape_util_strcmp(string1, string2, 1) returned %d, not 0 as expected.\n", status);
  }

  /* Test ape_util_cmp_string_array. */
  { const char * array1[] = { "a", "B", "c", 0 };
    const char * array2[] = { "a", "b", "e", 0 };
    status = ape_util_cmp_string_array(0, 0, 0);
    if (0 != status) ape_test_failed("ape_util_cmp_string_array(0, 0, 0) returned %d, not 0 as expected.\n", status);

    status = ape_util_cmp_string_array(array1, 0, 1);
    if (0 >= status) ape_test_failed("ape_util_cmp_string_array(array1, 0, 0) returned %d, not < 0 as expected.\n", status);

    status = ape_util_cmp_string_array(0, array2, 1);
    if (0 <= status) ape_test_failed("ape_util_cmp_string_array(0, array2, 0) returned %d, not > 0 as expected.\n", status);

    status = ape_util_cmp_string_array(array1, array2, 0);
    expected_status = strcmp(array1[1], array2[1]);
    if (expected_status != status)
      ape_test_failed("ape_util_cmp_string_array(array1, array2, 0) returned %d, not %d as expected.\n", status, expected_status);

    expected_status = 'c' - 'e';
    status = ape_util_cmp_string_array(array1, array2, 1);
    if (expected_status != status)
      ape_test_failed("ape_util_cmp_string_array(array1, array2, 1) returned %d, not %d as expected.\n", status, expected_status);
  }

  /* Test ape_util_find_string. */
  { char * array[] = { "eV", "keV", "MeV" };
    char ** begin = array;
    char ** end = array + sizeof(array)/sizeof(array[0]);
    char ** found = 0;
    const char * input = "GeV";
    int expected_status = eOK;
    size_t idx = 0;

    /* Search for a string that should not be found, case insensitively. */
    status = ape_util_find_string(begin, end, input, &found, 1);
    if (0 != found) {
      if (0 != *found) {
        ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned \"%s\", not 0 as expected.\n",
          input, *found);
      } else {
        ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned null string unexpectedly.\n", input);
      }
    }
    if (eOK == status) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) unexpectedly returned status %d.\n",
        input, status);
    }

    /* Search for a string that should not be found, case sensitively. */
    input = "ev";
    status = ape_util_find_string(begin, end, input, &found, 0);
    if (0 != found) {
      if (0 != *found) {
        ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) returned \"%s\", not 0 as expected.\n",
          input, *found);
      } else {
        ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) returned null string unexpectedly.\n", input);
      }
    }
    if (eOK == status) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) unexpectedly returned status %d.\n",
        input, status);
    }

    /* Search for a string that should be found, case insensitively. */
    input = "ev";
    status = ape_util_find_string(begin, end, input, &found, 1);
    if (0 == found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned null pointer unexpectedly.\n", input);
    } else if (0 == *found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned null string unexpectedly.\n", input);
    } else if (array[idx] != *found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned \"%s\", not \"%s\" as expected.\n",
        input, *found, array[idx]);
    }
    if (expected_status != status) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned status %d, not %d as expected.\n",
        input, status, expected_status);
    }

    /* Search for a string that should be found, case sensitively. */
    input = "keV";
    idx = 1;
    status = ape_util_find_string(begin, end, input, &found, 0);
    if (0 == found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) returned null pointer unexpectedly.\n", input);
    } else if (0 == *found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) returned null string unexpectedly.\n", input);
    } else if (array[idx] != *found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) returned \"%s\", not \"%s\" as expected.\n",
        input, *found, array[idx]);
    }
    if (expected_status != status) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 0) returned status %d, not %d as expected.\n",
        input, status, expected_status);
    }

    /* Search for a string that should be found at the last position. */
    input = "MeV";
    idx = 2;
    status = ape_util_find_string(begin, end, input, &found, 1);
    if (0 == found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned null pointer unexpectedly.\n", input);
    } else if (0 == *found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned null string unexpectedly.\n", input);
    } else if (array[idx] != *found) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned \"%s\", not \"%s\" as expected.\n",
        input, *found, array[idx]);
    }
    if (expected_status != status) {
      ape_test_failed("ape_util_find_string(begin, end, \"%s\", &found, 1) returned status %d, not %d as expected.\n",
        input, status, expected_status);
    }
  }

  /* Test ape_util_check_file_access and related functions. */
  { const char * file_name = "ape_test.par";
    int expected_status = eOK;
    int (*check_func)(const char *, const char *) = 0;

    /* Check file which can be read. */
    status = ape_util_check_file_access(file_name, "r");
    ape_test_cmp_long("ape_util_check_file_access(\"ape_test.par\", \"r\")", 0, 0, status, expected_status);

    file_name = "non-existent-file";
    expected_status = eFileNotAccessible;
    /* Check file which can't be read. */
    status = ape_util_check_file_access(file_name, "R");
    ape_test_cmp_long("ape_util_check_file_access(\"non-existent-file\", \"R\")", 0, 0, status, expected_status);

    file_name = "non-existent-file";
    expected_status = eOK;
    /* Check file which can be written. */
    status = ape_util_check_file_access(file_name, "w");
    ape_test_cmp_long("ape_util_check_file_access(\"non-existent-file\", \"w\")", 0, 0, status, expected_status);

    file_name = "non-existent-directory/spud";
    expected_status = eFileNotAccessible;
    /* Check file which can't be written. */
    status = ape_util_check_file_access(file_name, "W");
    ape_test_cmp_long("ape_util_check_file_access(\"non-existent-directory/spud\", \"W\")", 0, 0, status, expected_status);

    file_name = "ape_test.par";
    expected_status = eOK;
    /* Check file which exists. */
    status = ape_util_check_file_access(file_name, "");
    ape_test_cmp_long("ape_util_check_file_access(\"ape_test.par\", \"\")", 0, 0, status, expected_status);

    file_name = "non-existent-file";
    expected_status = eFileNotAccessible;
    /* Check a file that does not exist for existence. */
    status = ape_util_check_file_access(file_name, "e");
    ape_test_cmp_long("ape_util_check_file_access(\"non-existent-file\", \"\")", 0, 0,
      status, expected_status);

    /* Test file checker getter to make sure it returns correct result. */
    status = ape_util_get_file_check_func(&check_func);
    test_cmp_check_ptr("Before changing file check function, ape_util_get_file_check_func(&check_func)",
      check_func, check_file_access, status, eOK);

    /* Set file checker to function which returns eOK no matter what. */
    status = ape_util_set_file_check_func(&dont_check_access);
    ape_test_cmp_long("ape_util_set_file_check_func(&dont_check_access)", 0, 0, status, eOK);

    /* Test file checker getter to make sure this worked. */
    status = ape_util_get_file_check_func(&check_func);
    test_cmp_check_ptr("After changing file check function, ape_util_get_file_check_func(&check_func)",
      check_func, dont_check_access, status, eOK);

    /* Repeat failure cases from above; they should "succeed" in that they do not check access, but return true always. */
    file_name = "non-existent-file";
    expected_status = eOK;
    status = ape_util_check_file_access(file_name, "r");
    ape_test_cmp_long("After ape_util_set_file_check_func, ape_util_check_file_access(\"non-existent-file\", \"r\")", 0, 0,
      status, expected_status);

    file_name = "/spud";
    expected_status = eOK;
    status = ape_util_check_file_access(file_name, "w");
    ape_test_cmp_long("After ape_util_set_file_check_func, ape_util_check_file_access(\"/spud\", \"w\")", 0, 0,
      status, expected_status);

    file_name = "ape_test.par";
    expected_status = eOK;
    status = ape_util_check_file_access(file_name, "");
    ape_test_cmp_long("After ape_util_set_file_check_func, ape_util_check_file_access(\"ape_test.par\", \"\")", 0, 0,
      status, expected_status);

    file_name = "non-existent-file";
    expected_status = eOK;
    status = ape_util_check_file_access(file_name, "");
    ape_test_cmp_long("After ape_util_set_file_check_func, ape_util_check_file_access(\"non-existent-file\", \"\")", 0, 0,
      status, expected_status);

    /* Restore default checking behavior. */
    status = ape_util_set_file_check_func(0);
    ape_test_cmp_long("ape_util_set_file_check_func(0)", 0, 0, status, eOK);

    /* Check whether restore worked. */
    status = ape_util_get_file_check_func(&check_func);
    test_cmp_check_ptr("After restoring file check function, ape_util_get_file_check_func(&check_func)",
      check_func, check_file_access, status, eOK);
  }

  /* Test ape_util_expand_env_var. */
  { const char * input = 0;
    char * output = 0;
    const char * expected = 0;
    char * expanded = 0;

    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eNullPointer);
    free(output); output = 0;

    input = "";
    output = 0;
    expected = "";
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "No environment variables at all.";
    expected = input;
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "$ENV{NOTANENVVARNAME}";
    expected = "";
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "$ENV{NOTANENVVARNAME)";
    expected = input;
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "$ENV(PATH)";
    expected = getenv("PATH");
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "Singly quoted empty env variable '$ENV{NOTANENVVARNAME}'.";
    expected = input;
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "Doubly quoted empty env variable \"$ENV(NOTANENVVARNAME)\".";
    expected = "Doubly quoted empty env variable \"\".";
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "Escaped env variable \\$ENV{NOTANENVVARNAME}.";
    expected = input;
    status = ape_util_expand_env_var(input, &output);
    ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
    free(output); output = 0;

    input = "$ENV{PATH}$(PATH)\\$ENV{PATH}";
    { const char * value = getenv("PATH");
      expanded = 0;
      value = getenv("PATH");
      status = ape_util_cat_string(value, value, &expanded);
      if (eOK == status) {
        value = expanded;
        status = ape_util_cat_string(value, "\\$ENV{PATH}", &expanded);
        if (value != expanded) free((char *) value);
      }
      value = 0;
    }

    expected = expanded;
    if (eOK == status) {
      status = ape_util_expand_env_var(input, &output);
      ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
      free(output); output = 0;
    } else {
      ape_test_failed("Unable to set up to test ape_util_expand_env_var(\"$ENV{PATH}$(PATH)\\$ENV{PATH}\")", status, eOK);
    }
    free(expanded); expanded = 0;

    input = "\\$ENV{PATH}$ENV{PATH}$ENV(PATH)";
    { const char * value = getenv("PATH");
      expanded = 0;
      value = getenv("PATH");
      status = ape_util_cat_string(value, value, &expanded);
      if (eOK == status) {
        value = expanded;
        status = ape_util_cat_string("\\$ENV{PATH}", value, &expanded);
        if (value != expanded) free((char *) value);
      }
      value = 0;
    }

    expected = expanded;
    if (eOK == status) {
      status = ape_util_expand_env_var(input, &output);
      ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
      free(output); output = 0;
    } else {
      ape_test_failed("Unable to set up to test ape_util_expand_env_var(\"\\$ENV{PATH}$ENV{PATH}$ENV(PATH)\")", status, eOK);
    }
    free(expanded); expanded = 0;

    input = "\\${PATH}$ENV{PATH}${PATH}";
    { const char * value = getenv("PATH");
      expanded = 0;
      value = getenv("PATH");
      status = ape_util_cat_string(value, value, &expanded);
      if (eOK == status) {
        value = expanded;
        status = ape_util_cat_string("\\${PATH}", value, &expanded);
        if (value != expanded) free((char *) value);
      }
      value = 0;
    }

    expected = expanded;
    if (eOK == status) {
      status = ape_util_expand_env_var(input, &output);
      ape_test_cmp_string("ape_util_expand_env_var", output, expected, status, eOK);
      free(output); output = 0;
    } else {
      ape_test_failed("Unable to set up to test ape_util_expand_env_var(\"\\${PATH}$ENV{PATH}${PATH}\")", status, eOK);
    }
    free(expanded); expanded = 0;
  }

}

/* Unit test of utilities. */
void ape_util_test(void) {
  test_conversion();
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_util.c,v $
 * Revision 1.56  2011/02/01 18:01:47  jpeachey
 * Add support and tests for functions to convert string values into int
 * and short types, and to convert a parameter value into short type.
 *
 * Revision 1.55  2010/11/23 22:05:56  jpeachey
 * Replace tabs with spaces to make code line up correctly.
 *
 * Revision 1.54  2010/11/23 21:35:53  jpeachey
 * Cast to int to avoid warnings about converting size_t to int.
 *
 * Revision 1.53  2010/11/23 19:22:21  jpeachey
 * Move read_from_stdin static function from ape_par to ape_util. This is
 * only used in the readline-free version.
 *
 * Revision 1.52  2010/11/12 20:53:03  irby
 * Moved some internal static functions from ape_par to ape_util, and allowed
 * new custom_get_text routine to be used instead of ape's standard getter in
 * order to handle xpi call-back.
 *
 * Revision 1.51  2009/06/22 20:44:22  peachey
 * Use Windows's Sleep function, which takes milliseconds instead of
 * seconds, but otherwise behaves like Unix sleep.
 *
 * Revision 1.50  2009/06/11 19:09:55  peachey
 * Add ape_util_sleep function, needed for unit test of ape_io with
 * new behavior for hidden parameters (system time stamp overrules local
 * hidden parameters.)
 *
 * Revision 1.49  2007/11/12 19:45:19  peachey
 * Add ape_util_interpret_env function for setting internal state of Ape
 * based on environment variables. Remove spaces at end of line.
 *
 * Revision 1.48  2007/10/10 20:16:40  peachey
 * Add and test ape_util_get_file_check, for retrieving pointer to current
 * file checking function.
 *
 * Revision 1.47  2007/10/10 14:50:07  peachey
 * Modify Windows-specific file access checking to give consistentbehavior with Unix version.
 *
 * Revision 1.46  2007/10/09 18:20:59  peachey
 * Handle checking for existent and non-existent files (fe and fn types).
 *
 * Revision 1.45  2007/09/17 17:07:31  peachey
 * Changed Windows-specific file check to be consistent with Unix version.
 *
 * Revision 1.44  2007/09/17 16:36:32  peachey
 * Use access(...) function instead of fopen to check file access.
 *
 * Revision 1.43  2006/12/22 20:30:52  peachey
 * Tweak unit test for unwritable file check by including an invalid
 * directory in the path.
 *
 * Revision 1.42  2006/12/21 21:44:29  peachey
 * Add ape_util_find_string, to find the first match to a string
 * in an array of strings.
 *
 * Revision 1.41  2006/12/01 19:28:48  peachey
 * Fix a bug: pass the address of temporary variable.
 *
 * Revision 1.40  2006/11/30 20:38:35  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.39  2006/11/24 19:48:00  peachey
 * Add ape_util_s2f.
 *
 * Revision 1.38  2006/11/13 17:33:31  peachey
 * Generalize environment variable expansion to recognize ${var}, $(var) and $ENV(var)
 * in addition to $ENV{var}.
 *
 * Revision 1.37  2006/11/08 22:51:26  peachey
 * Revert to v 1.35; the changes had unintended consequences for the headas file checker.
 *
 * Revision 1.36  2006/11/07 15:04:30  peachey
 * Handle blank file name correctly.
 *
 * Revision 1.35  2006/06/23 20:45:59  peachey
 * Correct a bug which crept in during refactoring ape_par_check_file_access.
 * A file with a blank access mode is always considered to exist.
 *
 * Revision 1.34  2006/06/23 19:03:41  peachey
 * Removed code for NaNs which is no longer needed.
 *
 * Revision 1.33  2006/06/23 19:01:25  peachey
 * Handle Nan by returning status == eNan instead of messing with NaN
 * representations. Make sure test code for too-long constants works on 64 bit.
 *
 * Revision 1.32  2006/06/23 03:16:15  peachey
 * Rationalize file access checks to return eOK if file is accessible.
 *
 * Revision 1.31  2006/06/23 02:48:03  peachey
 * Use finite rather than isinf; the former seems to be more portable.
 *
 * Revision 1.30  2006/06/22 20:25:41  peachey
 * Use $ENV{VAR_NAME} instead of simply $VAR_NAME for environment variables.
 *
 * Revision 1.29  2006/06/20 17:20:24  peachey
 * Correct bug which was preventing env variables names with _ from being handled.
 *
 * Revision 1.28  2006/06/15 13:37:09  peachey
 * Add infrastructure to support expansion of environment variables.
 *
 * Revision 1.27  2006/06/07 06:43:00  peachey
 * Further tweak the file checker to work correctly for blank fopen modes.
 *
 * Revision 1.26  2006/06/07 05:54:04  peachey
 * Reverse the output of ape_util_file_check_access: 1 means access is OK,
 * 0 means it isn't.
 *
 * Revision 1.25  2006/06/07 02:43:48  peachey
 * Prevent checking file access if access type is not defined (return eOK).
 *
 * Revision 1.24  2006/06/05 01:30:48  peachey
 * Include file access checking facilities, including ability for client
 * to supply a custom access checking function.
 *
 * Revision 1.23  2006/05/31 03:57:10  peachey
 * Remove ape_strtod and ape_strtol which are no longer needed.
 *
 * Revision 1.22  2006/05/31 03:53:59  peachey
 * Set status to eUndefinedValue if numeric parameters have special values
 * indef, undef, none, etc.
 *
 * Revision 1.21  2006/05/18 03:17:07  peachey
 * Add ape_util_append_file_name for constructing full file names
 * from directory name + file name. Also improve handling of indef, nan etc.
 * in numeric conversions.
 *
 * Revision 1.20  2006/05/01 19:13:20  peachey
 * Improve conditional definitions for pfiles/syspfiles dividers.
 *
 * Revision 1.19  2006/05/01 17:16:55  peachey
 * Fix unsolved problem with unit test in which strcmp returned different
 * results for the same inputs on OSX.
 *
 * Revision 1.18  2006/05/01 15:13:14  peachey
 * Use native strcmp for case-sensitive string comparisons in ape_util_strcmp.
 *
 * Revision 1.17  2006/04/28 02:50:59  peachey
 * Fix sloppy implementation of ape_util_strcmp and add similar
 * comparison for arrays of strings, ape_util_cmp_string_array.
 *
 * Revision 1.16  2006/04/26 01:30:02  peachey
 * Add ape_util_free_string, a utility for freeing arrays of char *.
 *
 * Revision 1.15  2006/04/22 01:35:07  peachey
 * Add public function ape_util_copy_range.
 *
 * Revision 1.14  2006/04/21 14:26:13  peachey
 * Add ape_util_strcmp, for comparing strings with optional case insensitivity.
 *
 * Revision 1.13  2006/04/19 15:44:15  peachey
 * Add atexit facilities: public function ape_util_atexit plus static
 * function ape_all_atexit. The latter is the only ape function registered with
 * the system atexit function, but it then calls all other functions which
 * were registered using ape_util_atexit.
 *
 * Revision 1.12  2006/04/14 03:45:17  peachey
 * Change ape_util_copy_range so that it tolerates copying a null pointer
 * to a null output pointer.
 *
 * Revision 1.11  2006/04/13 18:43:39  peachey
 * Add and test ape_util_cat_string, for concatenating strings.
 *
 * Revision 1.10  2006/04/12 18:13:42  peachey
 * Add WIN32 version of test for non-trivial PFILES.
 *
 * Revision 1.9  2006/04/12 17:59:20  peachey
 * Add ape_util_parse_pfiles, for splitting into local and system parts of PFILES.
 *
 * Revision 1.8  2006/04/12 14:17:24  peachey
 * Add and test ape_util_getenv.
 *
 * Revision 1.7  2006/04/11 20:39:41  peachey
 * Make nan still more portable.
 *
 * Revision 1.6  2006/04/11 20:30:31  peachey
 * Make nan more portable.
 *
 * Revision 1.5  2006/04/11 17:20:24  peachey
 * Refine the conversion process for obtaining a long, in order to detect overflow
 * and underflow situations more accurately.
 *
 * Revision 1.4  2006/04/11 16:43:41  peachey
 * Add Windows-specific code to handle isnan/_isnan and isinf.
 *
 * Revision 1.3  2006/04/11 16:38:57  peachey
 * Use 0 not 0. to initialize long variable.
 *
 * Revision 1.2  2006/04/11 01:03:19  peachey
 * Add ape_util_s2d and its unit test. Improve handling of overflows to silence osx warning.
 *
 * Revision 1.1  2006/04/10 21:10:56  peachey
 * New module for general utilities, such as string conversions.
 *
*/
