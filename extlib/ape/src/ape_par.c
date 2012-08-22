/** \file ape_par.h
    \brief Declaration of parameter facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_list.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_test.h"
#include "ape/ape_util.h"

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

struct ApePar {
  char * comment;
  ApeList * field; /* A list of lists of char *. */
  int prompt_style;
  char cmd_line; /* Flag indicating the parameter was set on command line. */
  char modified; /* Flag indicating the parameter was modified since being read from file. */
};

enum {
  eTokenPrefix,
  eTokenLeadingQuote,
  eTokenValue,
  eTokenTrailingQuote,
  eTokenSuffix,
  eTokenNull,
  eNumTokens
};

static int compare_par(const ApePar * ape_par, const char * correct_line);
static void test_one_par(const char * line, const char ** correct_field);
static void test_par_check(const char * par_text, int expected_status, char check_value);
static int parse_par(const char * line, ApePar * ape_par);
static int parse_field(const char * begin, const char * end, char ** token_list);
static int find_field(const ApePar * ape_par, ParFieldId field_id, char *** token_list, int * quote);
static int get_field(const ApePar * ape_par, ParFieldId field_id, char ** field_string);
/* TODO: reimplement this to support vectors. */
/* static int get_field_array(const ApePar * ape_par, ParFieldId field_id, char *** field_array); */
static int set_field(const ApePar * ape_par, ParFieldId field_id, const char * field_text);
static const char * find_unquoted(const char * s, const char * delim);
static int collapse_escape(const char * input, const char * sequence, char ** output);
static int expand_escape(const char * begin, const char * end, const char * escape_me, char ** output);
static int quote_field(const ApePar * ape_par, ParFieldId field_id);
static int is_blank(const char * s);
static ApeParGetTextFunc custom_get_text = 0;

enum { eNoQuote, eSingleQuote, eDoubleQuote };

int ape_par_clone(ApePar * orig_par, ApePar ** clone_par) {
  int status = eOK;
  char * orig_text = 0;

  if (0 != clone_par) {
    /* Initialize output. */
    *clone_par = 0;
  } else {
    /* Output pointer is not valid. */
    status = eNullPointer;
  }

  if (eOK == status && 0 == orig_par) {
    /* Input pointer is not valid. */
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Get text of original parameter. */
    orig_text = ape_par_get_line(orig_par);
    if (0 == orig_text) {
      status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* Create the clone. */
    status = ape_par_create(clone_par, orig_text);
  }

  if (eOK == status) {
    /* Copy the prompt style. */
    (*clone_par)->prompt_style = orig_par->prompt_style;
  }

  free(orig_text); orig_text = 0;

  return status;
}

/* Look at the current character to determine if it is a quote, taking into account the escape character. */
static char detect_quote(const char ** cp) {
  char result = 0;
  if (0 != cp && 0 != *cp) {
    if ('\\' == **cp) {
      /* Ignore the next character after backslash (escape) by going on to the next character unless it is 0. */
      if ('\0' != (*cp)[1]) ++(*cp);
    } else if ('\'' == **cp) {
      /* An unescaped single quote was detected. */
      result = eSingleQuote;
    } else if ('"' == **cp) {
      /* An unescaped double quote was detected. */
      result = eDoubleQuote;
    }
  }

  return result;
}

static const char * find_comment(const char * line) {
  char in_quote = eNoQuote;
  int field_index = 1;
  char found_quote = '\0';

  /* Skip leading space. */
  for(; 0 != isspace(*line); ++line) {}

  /* A quote either must start here... */
  if ('#' == *line) return line;

  /* ... or after all required fields were found. */
  for(; '\0' != *line; ++line) {
    if ('\\' == *line) {
      /* Skip this character and the next, if any. */
      if ('\0' != *(line + 1)) ++line;
      continue;
    }
    found_quote = detect_quote(&line);
    if (eNoQuote != found_quote) {
      /* If not already inside a quote, start one of the correct type. */
      if (eNoQuote == in_quote) in_quote = found_quote;
      /* If in a quote of the correct type, close the quote. */
      else if (in_quote == found_quote) in_quote = eNoQuote;
    } else if (eNoQuote == in_quote) {
      /* An unquoted comma indicates a new field was found.
         Don't look for # unless all required fields were found. */
      if (',' == *line) ++field_index;
      else if (eEndOfField > field_index) continue;
      else if ('#' == *line) break;
    }
  }
  return line;
}

int ape_par_create(ApePar ** ape_par_ptr, const char * line) {
  int status = eOK;
  if (0 == ape_par_ptr || 0 == line) {
    if (0 == ape_par_ptr)
      ape_msg_debug("ape_par_create: %s.\n", "NULL pointer passed for first argument \"ape_par_ptr\"");
    if (0 == line) {
      ape_msg_debug("ape_par_create: %s.\n", "NULL pointer passed for second argument \"line\"");
    } else {
      ape_msg_debug("ape_par_create: %s.\n", "line was");
      ape_msg_debug("ape_par_create: %s.\n", line);
    }
    status = eNullPointer;
  } else {
    /* Allocate space for the parameter and initialize it to 0. */
    *ape_par_ptr = (ApePar *) calloc(1, sizeof(ApePar));

    if (0 == *ape_par_ptr) {
      ape_msg_debug("ape_par_create: %s.\n", "allocation failed for new ApePar");
      ape_msg_debug("ape_par_create: %s.\n", "line was");
      ape_msg_debug("ape_par_create: %s.\n", line);
      status = eDynAllocFailed;
    } else {
      const char * comment = find_comment(line);
      char * no_comment = 0;
      status = ape_util_copy_range(line, comment, &no_comment);
      if (eOK == status) {
        status = ape_util_copy_string(comment, &(*ape_par_ptr)->comment);
      }
      if (eOK == status) {
        parse_par(no_comment, *ape_par_ptr);
      }
      free(no_comment); no_comment = 0;
      (*ape_par_ptr)->prompt_style = eDefaultPrompt;
    }
  }

  return status;
}

void ape_par_destroy(ApePar * ape_par) {
  if (0 != ape_par) {
    ape_par->prompt_style = eDefaultPrompt;
    ape_par->modified = 0;
    ape_par->cmd_line = 0;
    if (0 != ape_par->field) {
      ApeListIterator field_itor = 0;
      /* Loop backwards over the list, destroying contents of each element. */
      for (field_itor = ape_list_end(ape_par->field); field_itor != ape_list_begin(ape_par->field); ) {
        char ** token_list = 0;

        /* Decrement the iterator here, because iteration is in reverse order. */
        field_itor = ape_list_prev(field_itor);

        /* Get the list of tokens for this field. */
        token_list = (char **) ape_list_get(field_itor);

        ape_util_free_string_array(token_list); token_list = 0;
      }

      /* Destroy list of fields. */
      ape_list_destroy(ape_par->field);
      ape_par->field = 0;
    }

    /* Free the comment string. */
    free(ape_par->comment); ape_par->comment = 0;

    /* Free the memory held by the parameter. */
    free(ape_par);
  }
}

char * ape_par_get_line(const ApePar * ape_par) {
  char * line = 0;
  if (0 == ape_par) {
    ape_msg_debug("ape_par_get_line: %s.\n", "NULL pointer passed for argument \"ape_par\"");
  } else if (0 == ape_par->field) {
    ape_msg_debug("ape_par_get_line: %s.\n", "parameter's list of fields is empty");
  } else {
    ApeListIterator field_itor = 0;
    ApeList * field = ape_par->field;
    ApeListIterator end = ape_list_end(field);
    size_t line_size = 1; /* For terminating NULL. */

    /* Loop once just to get the total size of the line. */
    for (field_itor = ape_list_begin(field); field_itor != end; field_itor = ape_list_next(field_itor)) {
      char ** list = (char **) ape_list_get(field_itor);
      for (; 0 != list && 0 != *list; ++list)
        line_size += strlen(*list);
      /* Leave room for comma delimiter. */
      ++line_size;
    }

    /* Add room for the comment. */
    line_size += strlen(ape_par->comment);

    /* Create output buffer. */
    line = (char *) calloc(line_size, sizeof(char));
    if (0 == line) {
      ape_msg_debug("ape_par_get_line: %s.\n", "allocation failed for output line buffer");
    } else {
      /* Loop again to concatenate all the parameters to the output line. */
      for (field_itor = ape_list_begin(field); field_itor != end; field_itor = ape_list_next(field_itor)) {
        char ** list = (char **) ape_list_get(field_itor);
        for (; 0 != list && 0 != *list; ++list)
          strcat(line, *list);
        /* At end of field, insert terminating comma. */
        strcat(line, ",");
      }

      /* An extra comma was concatenated to the output string, so zero it out. */
      line[strlen(line) - 1] = '\0';

      /* Append comment. */
      strcat(line, ape_par->comment);
    }
  }
  return line;
}

int ape_par_get_name(const ApePar * ape_par, char ** name) {
  int status = eOK;

  /* Check arguments and initialize where possible. */
  if (0 == ape_par || 0 == name) status = eNullPointer;
  if (eOK == status) *name = 0;

  /* Get the name field. */
  if (eOK == status) status = ape_par_get_field(ape_par, eName, name);

  /* Check whether name is blank. */
  if (eOK == status && 0 != is_blank(*name)) {
    status = eUnnamedPar;
    free(*name); *name = 0;
  }

  return status;
}
 
int ape_par_get_bool(const ApePar * par, char * result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0;
    status = ape_par_get_field(par, eValue, &str_result);
    if (eOK == status) {
      status = ape_util_s2b(str_result, result);
    }
    free(str_result); str_result = 0;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_double(const ApePar * par, double * result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0.;
    status = ape_par_get_field(par, eValue, &str_result);
    if (eOK == status) {
      status = ape_util_s2d(str_result, result);
    }
    free(str_result); str_result = 0;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_float(const ApePar * par, float * result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0.;
    status = ape_par_get_field(par, eValue, &str_result);
    if (eOK == status) {
      status = ape_util_s2f(str_result, result);
    }
    free(str_result); str_result = 0;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_file_name(const ApePar * par, char ** result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0;
    status = ape_par_get_field(par, eValue, &str_result);
    *result = str_result;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_int(const ApePar * par, int * result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0;
    status = ape_par_get_field(par, eValue, &str_result);
    if (eOK == status) {
      status = ape_util_s2i(str_result, result);
    }
    free(str_result); str_result = 0;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_long(const ApePar * par, long * result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0;
    status = ape_par_get_field(par, eValue, &str_result);
    if (eOK == status) {
      status = ape_util_s2l(str_result, result);
    }
    free(str_result); str_result = 0;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_short(const ApePar * par, short * result) {
  int status = eOK;
  if (0 != result) {
    char * str_result = 0;
    int check_status = eOK;
    *result = 0;
    status = ape_par_get_field(par, eValue, &str_result);
    if (eOK == status) {
      status = ape_util_s2sh(str_result, result);
    }
    free(str_result); str_result = 0;
    check_status = ape_par_check(par, 1);
    status = eOK == status ? check_status : status;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_string(const ApePar * par, char ** result) {
  int status = eOK;
  char * min_string = 0;

  /* Check arguments. */
  if (0 == par || 0 == result) status = eNullPointer;

  /* TODO Fix fthedit and other tools that depend on enumerated strings being upper-cased, then just call
     ape_par_get_string_case. */
  if (eOK == status) {
    /* Get minimum string just in order to learn whether this is an enumerated parameter. */
    status = ape_par_get_min_string(par, &min_string);
    free(min_string); min_string = 0;
  }

  if (eOK == status) {
    /* For non-enumerated strings, use case supplied by the user. */
    status = ape_par_get_string_case(par, result, eDefaultCase);
  } else if (eRangeEnum == status) {
    /* For enumerated range, convert to upper case. */
    status = ape_par_get_string_case(par, result, eUpperCase | eEnumCase);
  }

  return status;
}

int ape_par_get_string_case(const ApePar * par, char ** result, char case_code) {
  int status = eOK;
  char * str_result = 0;
  int check_status = eOK;

  /* Check arguments. */
  if (0 == par || 0 == result) status = eNullPointer;

  if (eOK == status) {
    /* Initialize output, then attempt to get the current value of the parameter. */
    *result = 0;
    status = ape_par_get_field(par, eValue, &str_result);
  }

  if (eOK == status) {
    if (0 != (eEnumCase & case_code)) {
      /* Try to match current value to one of the enumerated values. If a match is found,
         use the enumerated value and not the string the user entered. */
      char ** range = 0;
      char ** found = 0;
      status = ape_par_get_enum_string(par, &range);
      if (eOK == status) {
        char ** begin = range;
        char ** end = begin;
        /* Find end of range array. */
        for (; 0 != *end; ++end);
        status = ape_util_find_string(begin, end, str_result, &found, 1);
        if (eOK == status) {
          free(str_result); str_result = 0;
          status = ape_util_copy_string(*found, &str_result);
        }
      } else if (eRangeNoEnum == status) {
        /* If this is not an enumerated parameter, don't worry about it. */
        status = eOK;
      }
      ape_util_free_string_array(range); range = 0;

      /* Done with enumerated case handling, so mask off the enumerated case flag. */
      case_code &= ~eEnumCase;
    }
  }
  if (eOK == status) {
    /* Handle remaining switches; they are mutually exclusive. */
    switch (case_code) {
      case eUnknownCase:
      case eDefaultCase:
      case eEnumCase: /* Note this should not happen, as this bit was masked off above. */
        break;
      case eLowerCase: {
          char * cp = 0;
          status = eOK;
          for (cp = str_result; '\0' != *cp; ++cp) *cp = tolower(*cp);
        }
        break;
      case eUpperCase: {
          char * cp = 0;
          status = eOK;
          for (cp = str_result; '\0' != *cp; ++cp) *cp = toupper(*cp);
        }
        break;
      default:
        status = eInvalidArgument;
        break;
    }
  }

  /* No matter what, return the best effort at retrieving the string, in case the caller decides to ignore some error. */
  *result = str_result;

  check_status = ape_par_check(par, 1);
  status = eOK == status ? check_status : status;

  return status;
}

int ape_par_get_comment(const ApePar * par, char ** comment) {
  int status = eOK;

  /* Check arguments. */
  if (0 == par || 0 == comment) status = eNullPointer;

  if (eOK == status) {
    /* Get original comment from parameter object, or empty string. */
    const char * orig_comment = 0 != par->comment ? par->comment : "";

    /* Initialize output to null. */
    *comment = 0;

    /* Copy original comment to output string. */
    status = ape_util_copy_string(orig_comment, comment);
  }

  return status;
}

int ape_par_get_field(const ApePar * par, ParFieldId id, char ** result) {
  int status = eOK;
  if (0 != result) {
    *result = 0;
    if (0 != par) {
      status = get_field(par, id, result);
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_get_field_array(const ApePar * par, ParFieldId id, char *** result) {
  int status = eOK;

  /* Check arguments. */
  if (0 != result) {
    *result = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == par) {
    status = eNullPointer;
  }

  if (eOK == status && (id < eName || id > ePrompt)) {
    status = eInvalidArgument;
  }

  if (eOK == status) {
/* TODO: reimplement this to support vectors.
    status = get_field_array(par, id, result);
*/
  }

  return status;
}

enum { eNoMatch = -1 };

static int encode_msg(const char * msg, const char * code_char, size_t num_char,
  const size_t ** known_code, const char ** known_msg, size_t num_code, char * encoded_msg) {
  int status = eOK;

  /* Resultant codes found in the input string. */
  size_t * found_code = (size_t *) calloc(num_char, sizeof(size_t));

  /* Flag indicating presence of multiple known characters. */
  char found_multiple = 0;

  /* Flag indicating presence of an unknown character. */
  char found_unknown = 0;

  size_t code_idx = 0;
  size_t known_idx = 0;

  /* Iterate over the string containing input. */
  for (; '\0' != *msg; ++msg) {
    /* Ignore white space. */
    if (0 != isspace(*msg)) continue;

    /* Iterate over known code characters; flag any which are found in the string. */
    for (code_idx = 0; num_char != code_idx; ++code_idx) {
      if (tolower(*msg) == code_char[code_idx]) {
        /* Check if this code character was already found one or more times. */
        if (0 != found_code[code_idx]) found_multiple = 1;

        /* Flag this code as found. */
        found_code[code_idx] = 1;
        break;
      }
    }

    /* See if the loop above exited because we exhausted the known characters without finding a match.
       If so, this character is not in the set of understood codes, so flag it as "unknown". */
    if (num_char == code_idx) found_unknown = 1;
  }

  /* At this point the string has been parsed and found_code contains the result; compare it against known codes. */
  for (known_idx = 0; num_code != known_idx; ++known_idx) {
    for (code_idx = 0; num_char != code_idx; ++code_idx) {
      if (found_code[code_idx] != known_code[known_idx][code_idx]) break;
    }
    /* If (inside) loop above exited because it reached the end, a match was found, so break out of the loop. */
    if (num_char == code_idx) break;
  }

  /* If (outside) loop above exited because it reached the end, *no* match was found, so the code is not
     a known code. */
  if (num_code == known_idx) {
    status = eNoMatch;
  } else {
    /* Although there may have been problems, msg was successfully encoded, so return the correct canonical string. */
    strcat(encoded_msg, known_msg[known_idx]);

    /* Check for lesser errors, which the client may choose to ignore. */
    if (found_unknown || found_multiple) status = eFormatError;
  }

  /* Clean up. */
  free(found_code); found_code = 0;

  return status;
}

/* Parse the input string and determine the canonical "code" which is in standard form. */
static int string2mode(const char * code_string, char * mode_string) {
  /* eNumCodes is the number of reserved characters, which are known to be possible pieces of the mode string.
     eNumModes is the number of known modes which can result from combinations of these codes. */
  enum { eNumCodes = 4, eNumModes = 6 };

  int status = eOK;

  /* known_code is the set of known characters which may appear in the mode string. */
  const char known_code[eNumCodes] = { 'a', 'h', 'l', 'q' };

  /* mode_code is the set of known valid combinations of characters from known_code. */
  const size_t mode_code[][eNumCodes] = {
    { 1, 0, 0, 0 }, /* a */
    { 1, 0, 1, 0 }, /* al */
    { 0, 1, 0, 0 }, /* h */
    { 0, 1, 1, 0 }, /* hl */
    { 0, 0, 0, 1 }, /* q */
    { 0, 0, 1, 1 }, /* ql */
  };

  /* known_mode is the set of canonical mode strings, corresponding to the characters in mode_code. */
  const char * known_mode[eNumModes] = { "a", "al", "h", "hl", "q", "ql" };

  size_t idx = 0;
  const size_t ** mode_code_ptr = 0;

  memset(mode_string, '\0', APE_PAR_MODE_CODE_LEN);

  mode_code_ptr = (const size_t **) calloc(eNumModes, sizeof(const size_t *));
  if (0 == mode_code_ptr) {
    status = eDynAllocFailed;
  }

  if (eOK == status) {
    for (; idx != eNumModes; ++idx) {
      mode_code_ptr[idx] = mode_code[idx];
    }
    status = encode_msg(code_string, known_code, eNumCodes, mode_code_ptr, known_mode, eNumModes, mode_string);
    if (eNoMatch == status) status = eUnknownMode;
  }

  free((void *) mode_code_ptr); mode_code_ptr = 0;

  return status;
}

int ape_par_get_mode(const ApePar * par, char * mode_string) {
  int status = eOK;

  /* Check arguments. */
  if (0 != mode_string) {
    memset(mode_string, '\0', APE_PAR_MODE_CODE_LEN);
  } else {
    status = eNullPointer;
  }
  if (eOK == status && 0 == par) status = eNullPointer;

  if (eOK == status) {
    char * code_string = 0;

    /* Get the actual code string. */
    status = ape_par_get_field(par, eMode, &code_string);
    if (eOK == status) {
      status = string2mode(code_string, mode_string);
    }
    free(code_string); code_string = 0;
  }
  return status;
}

int ape_par_get_eff_mode(const ApePar * par, const char * auto_string, char * mode_string) {
  int status = eOK;
  int prompt_style = eDefaultPrompt;

  /* Start by getting the innate mode of the parameter. */
  status = ape_par_get_mode(par, mode_string);

  if (eOK == status) {
    if ('a' == *mode_string) {
      char auto_mode[APE_PAR_MODE_CODE_LEN] = "";

      /* Parse the string given to determine the automatic mode. */
      status = string2mode(auto_string, auto_mode);
      if (eOK == status) {
        /* Substitute the mode parameter value for the "a" automatic mode code. */
        mode_string[0] = auto_mode[0];

        /* Substitute the mode's value for the "learn" field, but this only if the mode parameter
           turns *on* learning. */
        if ('l' == auto_mode[1]) mode_string[1] = auto_mode[1];
      }
    }
  }

  if (eOK == status) {
    /* Get the mode of the parameter, considering the default mode as appropriate. */
    status = ape_par_get_prompt_style(par, &prompt_style);
  }

  if (eOK == status && eDefaultPrompt == prompt_style) {
    status = ape_par_get_default_prompt_style(&prompt_style);
  }

  if (eOK == status) {
    /* If hidden parameter prompting is enabled for this parameter, effective mode is 'q' for all parameters. */
    if (0 != (eQueryHidden & prompt_style)) mode_string[0] = 'q';

    /* If this parameter was set on the command line, effective mode is 'h'. This trumps hidden parameter prompting. */
    if (0 != par->cmd_line) mode_string[0] = 'h';

    /* If multi-query mode is enabled for this parameter, effective mode is 'q' for all parameters. This trumps
       command-line prompt disabling. */
    if (0 != (eMultiQuery & prompt_style)) mode_string[0] = 'q';

  }

  return status;
}

int ape_par_get_min_string(const ApePar * par, char ** min) {
  int status = eOK;
  const char * found = 0;

  if (0 != min) {
    *min = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == par) {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Get the min field as a string. */
    status = ape_par_get_field(par, eMin, min);
  }

  if (eOK == status) {
    /* Check for existence of an unquoted pipe, which is used to indicate enumerated parameter. */
    found = find_unquoted(*min, "|");
    if (0 != found && '\0' != *found) status = eRangeEnum;
  }

  return status;
}

int ape_par_get_enum_string(const ApePar * par, char *** range) {
  int status = eOK;
  char * min = 0;
  const char * end = 0;
  size_t num_pipe = 0u;

  if (0 != range) {
    *range = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == par) {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Get the min field as a string. */
    status = ape_par_get_field(par, eMin, &min);
  }

  if (eOK == status) {
    /* Look for all unquoted pipes, at least one of which is needed to indicate an enumerated parameter range. */
    end = find_unquoted(min, "|");
    while (0 != end && '\0' != *end) {
      ++num_pipe;
      end = find_unquoted(end + 1, "|");
    }
  }

  if (eOK == status) {
    /* Make sure enough space is allocated. Number of enumerated strings = num_pipe + 1.
       Number of pointers needed = Number of enumerated strings + 1 for terminating 0.
       Therefore allocate num_pipe + 2 pointers, initializing them to 0. */
    *range = (char **) calloc(num_pipe + 2, sizeof(char *));
    if (0 == *range) {
      status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    size_t idx = 0;
    const char * begin = min;

    /* Skip leading white space. */
    while (0 != isspace(*begin)) ++begin;

    /* For sure there are some pipes, so find them again, this time copying strings to the output as we go. */
    end = find_unquoted(begin, "|");
    for (; eOK == status && 0 != end && '\0' != *end; ++idx) {
      const char * true_end = end;

      /* Skip trailing white space. */
      while (end > begin && 0 != isspace(*(end - 1))) --end;

      /* Copy the string to the output vector. */
      status = ape_util_copy_range(begin, end, *range + idx);

      /* Go on to next string. */
      begin = true_end + 1;

      /* Skip leading white space. */
      while (0 != isspace(*begin)) ++begin;

      end = find_unquoted(begin, "|");
    }

    if (eOK == status) {
      /* Skip trailing white space. */
      while (end > begin && 0 != isspace(*(end - 1))) --end;

      /* Copy the final string to the output vector. */
      status = ape_util_copy_range(begin, end, *range + idx);
    }
  }

  /* Clean up. */
  free(min); min = 0;

  /* Final status reset: flag non-enumerated ranges. */
  if (eOK == status && 0u == num_pipe) status = eRangeNoEnum;

  return status;
}

int ape_par_get_type(const ApePar * par, char * type_string) {
  int status = eOK;

  /* Check arguments. */
  if (0 != type_string) {
    memset(type_string, '\0', APE_PAR_TYPE_CODE_LEN);
  } else {
    status = eNullPointer;
  }

  if (0 == par) {
    status = eNullPointer;
  }

  if (eOK == status) {
    char * code_string = 0;

    /* Get the actual code string. */
    status = ape_par_get_field(par, eType, &code_string);
    if (eOK == status) {
      /* eNumCodes is the number of reserved characters, which are known to be possible pieces of the mode string.
         eNumTypes is the number of known types which can result from combinations of these codes. */
      enum { eNumCodes = 10, eNumTypes = 22 };
      /* known_code is the set of known characters which may appear in the type string. */
      /* TODO: support fully type 'g', used by xselect? */
      const char known_code[eNumCodes] = { 'b', 'd', 'e', 'f', 'g', 'i', 'n', 'r', 's', 'w' };

      /* type_code is the set of known valid combinations of characters from known_code. */
      const size_t type_code[eNumTypes][eNumCodes] = {
       /* b  d  e  f  g  i  n  r  s  w */
        { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /* b */
        { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /* d */
        { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /* f */
        { 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 }, /* fe */
        { 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /* fn */
        { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0 }, /* fr */
        { 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 }, /* fw */
        { 0, 0, 1, 1, 0, 0, 1, 0, 0, 0 }, /* fen */
        { 0, 0, 1, 1, 0, 0, 0, 1, 0, 0 }, /* fer */
        { 0, 0, 1, 1, 0, 0, 0, 0, 0, 1 }, /* few */
        { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0 }, /* fnr */
        { 0, 0, 0, 1, 0, 0, 1, 0, 0, 1 }, /* fnw */
        { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1 }, /* frw */
        { 0, 0, 1, 1, 0, 0, 1, 1, 0, 0 }, /* fenr */
        { 0, 0, 1, 1, 0, 0, 1, 0, 0, 1 }, /* fenw */
        { 0, 0, 1, 1, 0, 0, 0, 1, 0, 1 }, /* ferw */
        { 0, 0, 0, 1, 0, 0, 1, 1, 0, 1 }, /* fnrw */
        { 0, 0, 1, 1, 0, 0, 1, 1, 0, 1 }, /* fenrw */
        { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /* g */
        { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /* i */
        { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /* r */
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }  /* s */
      };

      /* known_type is the set of canonical type strings, corresponding to the characters in type_code. */
      const char * known_type[eNumTypes] = {
        "b", "d", "f", "fe", "fn", "fr", "fw",
        "fen","fer", "few", "fnr", "fnw", "frw",
        "fenr", "fenw", "ferw", "fnrw", "fenrw",
        "g", "i", "r", "s"
      };

      size_t idx = 0;
      const size_t ** type_code_ptr = (const size_t **) calloc(eNumTypes, sizeof(const size_t *));
      if (0 == type_code_ptr) {
        status = eDynAllocFailed;
      }

      if (eOK == status) {
        for (; idx != eNumTypes; ++idx) {
          type_code_ptr[idx] = type_code[idx];
        }
        status = encode_msg(code_string, known_code, eNumCodes, type_code_ptr, known_type, eNumTypes, type_string);
        memset((void *) type_code_ptr, '\0', eNumTypes * sizeof(const size_t *));
        if (eNoMatch == status) status = eUnknownType;
      }

      free((void *) type_code_ptr); type_code_ptr = 0;
    }
    free(code_string); code_string = 0;
  }
  return status;
}

int ape_par_set_field(ApePar * par, ParFieldId id, const char * field_text) {
  int status = eOK;
  if (0 == par || 0 == field_text) {
    status = eNullPointer;
  }
  if (eOK == status && (eName > id || ePrompt < id)) {
    status = eInvalidArgument;
  }
  if (eOK == status) {
    status = set_field(par, id, field_text);
  }
  return status;
}

int ape_par_set_value_string(ApePar * par, const char * value) {

  int status = eOK;
  const char * begin = value;
  const char * end = value + strlen(value);
  char * value_copy = 0;
  char par_type[APE_PAR_TYPE_CODE_LEN] = "";

  /* Prepare input; its leading and trailing white space is not significant. */
  while (begin < end && 0 != isspace(*begin)) ++begin;
  while (begin < end && 0 != isspace(*(end - 1))) --end;

  status = ape_util_copy_range(begin, end, &value_copy);
  if (eOK == status) {
    status = ape_par_set_field(par, eValue, value_copy);
    free(value_copy); value_copy = 0;
  }

  if (eOK == status) {
    par->modified = 1;
  }

  /* Determine type of this parameter. */
  if (eOK == status) status = ape_par_get_type(par, par_type);

  /* For string parameters only, re-escape the quotes around the value. This will add double quotes if they are not present. */
  if (eOK == status && ('f' == *par_type || 'g' == *par_type || 's' == *par_type)) status = quote_field(par, eValue);

  return status;
}

int ape_par_set_comment(ApePar * par, const char * comment) {
  int status = eOK;
  char * current_comment = 0;

  /* Check arguments. */
  if (0 == par || 0 == comment) {
    status = eNullPointer;
  }

  /* Clear out current comment value, being somewhat paranoid. */
  current_comment = par->comment;
  par->comment = 0;
  free(current_comment); current_comment = 0;

  /* Copy the comment given by the caller into the parameter's comment. */
  status = ape_util_copy_string(comment, &par->comment);

  return status;
}

static int s_default_prompt_style = eDefaultPrompt;

int ape_par_get_default_prompt_style(int * prompt_style) {
  int status = eOK;
  if (0 != prompt_style) {
    *prompt_style = s_default_prompt_style;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_set_default_prompt_style(int prompt_style) {
  int status = eOK;
  s_default_prompt_style = prompt_style;
  return status;
}

int ape_par_get_prompt_style(const ApePar * par, int * prompt_style) {
  int status = eOK;
  if (0 != par && 0 != prompt_style) {
    *prompt_style = par->prompt_style;
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_par_set_prompt_style(ApePar * par, int prompt_style) {
  int status = eOK;
  if (0 != par) {
    par->prompt_style = prompt_style;
  } else {
    status = eNullPointer;
  }
  return status;
}

/* \brief Determine whether the given error is due to an invalid parameter value.
   Such errors are sometimes recoverable, e.g. in the context of interactive prompting.
   \param par The parameter that has the status described by status.
   \param value_status The status.
*/
static int is_recoverable(ApePar * par, int value_status) {
  int status = eOK;
  switch (value_status) {
    case eOK:
    case eStringRemainder:
    case eTypeMismatch:
    case eConversionError:
    case eOverflow:
    case eUnderflow:
    case eValueBelowMin:
    case eValueAboveMax:
    case eInvalidChoice:
    case eFileNotAccessible:
      status = eOK;
      break;
    default:
      status = value_status;
      break;
  }

  return status;
}

/* \brief Describe the status in the context of the given parameter fields. Note that this
          is intended to be called only by ape_par_check.
   \param field The unpacked fields of the parameter: note there is no check for null parameter.
   \param status The status being described.
*/
static void describe_status(char ** field, int * field_status, int status) {
  const char * name = 0 != field[eName] ? field[eName] : "no name";
  const char * type = 0 != field[eType] ? field[eType] : "";
  const char * type_description = "unknown";
  const char * mode = 0 != field[eMode] ? field[eMode] : "unknown mode";
  const char * value = 0 != field[eValue] ? field[eValue] : "";
  const char * min = 0 != field[eMin] ? field[eMin] : "";
  const char * max = 0 != field[eMax] ? field[eMax] : "";
  char treat_as_string = 0;
  const char * fmt = "";

  if (eOK == status) return;

  /* Note that this function deliberately does not check status as usual because it
     is doing the best it can to report errors. */

  /* If possible, get the type of the parameter, for use in description. */
  if ('\0' != *type) {
    switch (*type) {
      case 'b':
        type_description = "boolean";
        break;
      case 'd':
        type_description = "double precision floating point";
        break;
      case 'f':
        type_description = "file";
        treat_as_string = 1;
        break;
      case 'g':
        type_description = "range";
        break;
      case 'i':
        type_description = "integer";
        break;
      case 'r':
        type_description = "floating point";
        break;
      case 's':
        type_description = "string";
        treat_as_string = 1;
        break;
      default:
        type_description = type;
        treat_as_string = 1;
        break;
    }
  }

  switch (status) {
    case eNullPointer:
      ape_msg_error("Ape: Program logic error: null pointer unexpectedly passed.\n");
      break;
    case eDynAllocFailed:
      ape_msg_error("Ape: Run-time error: failed to allocate dynamic memory.\n");
      break;
    case eInvalidArgument:
      ape_msg_error("Ape: Program logic error: invalid argument unexpectedly passed.\n");
      break;
    case eStringRemainder:
      ape_msg_error("Parameter %s: trailing text in value \"%s\". Cannot convert to %s type.\n", name, value, type_description);
      break;
    case eTypeMismatch:
      ape_msg_error("Parameter %s: value \"%s\" is not of %s type.\n", name, value, type_description);
      break;
    case eConversionError:
      ape_msg_error("Parameter %s: cannot convert value \"%s\" cleanly to %s type.\n", name, value, type_description);
      break;
    case eOverflow:
      ape_msg_error("Parameter %s: overflow converting value \"%s\" to %s type.\n", name, value, type_description);
      break;
    case eUnderflow:
      ape_msg_error("Parameter %s: underflow converting value \"%s\" to %s type.\n", name, value, type_description);
      break;
    case eVarNotSet:
      ape_msg_error("Parameter %s: unset environment variable in value \"%s\".\n", name, value);
      break;
    case eFieldNotFound: {
        size_t idx = 0;
        size_t num_missing = 0;
        ape_msg_error("Parameter %s: missing field(s)", name);
        for (idx = 0; idx != eEndOfField; ++idx) {
          if (eFieldNotFound == field_status[idx]) {
            ++num_missing;
            if (num_missing > 1) ape_msg_error(",");
            switch (idx) {
              case eName:
                ape_msg_error(" name");
                break;
              case eType:
                ape_msg_error(" type");
                break;
              case eMode:
                ape_msg_error(" mode");
                break;
              case eValue:
                ape_msg_error(" value");
                break;
              case eMin:
                ape_msg_error(" minimum");
                break;
              case eMax:
                ape_msg_error(" maximum");
                break;
              case ePrompt:
                ape_msg_error(" prompt");
                break;
              default:
                break;
            }
          }
        }
        ape_msg_error(".\n");
      }
      break;
    case eUnknownMode:
      ape_msg_error("Parameter %s: invalid mode \"%s\".\n", name, mode);
      break;
    case eUnknownType:
      ape_msg_error("Parameter %s: unknown type \"%s\".\n", name, field[eType]);
      break;
    case eFormatError:
      ape_msg_error("Parameter %s: invalid format.\n", name);
      break;
    case eInputFailure:
      ape_msg_error("Parameter %s: end of input encountered unexpectedly.\n", name);
      break;
    case eTooManyFields:
      ape_msg_error("Parameter %s: too many fields.\n", name);
      break;
    case eMinConversionError:
      ape_msg_error("Parameter %s: unable to convert minimum value \"%s\" to %s type.\n", name, min, type_description);
      break;
    case eMaxConversionError:
      ape_msg_error("Parameter %s: unable to convert maximum value \"%s\" to %s type.\n", name, max, type_description);
      break;
    case eValueBelowMin:
      if (0 == treat_as_string) fmt = "Parameter %s: value %s is less than minimum value %s.\n";
      else fmt = "Parameter %s: as a string, value \"%s\" is less than minimum value \"%s\".\n";
      ape_msg_error(fmt, name, value, min);
      break;
    case eValueAboveMax:
      if (0 == treat_as_string) fmt = "Parameter %s: value %s is greater than maximum value %s.\n";
      else fmt = "Parameter %s: as a string, value \"%s\" is less than maximum value \"%s\".\n";
      ape_msg_error(fmt, name, value, max);
      break;
    case eInvalidRange:
      if (0 == treat_as_string) fmt = "Parameter %s: invalid range (%s, %s).\n";
      else fmt = "Parameter %s: invalid range (\"%s\", \"%s\").\n";
      ape_msg_error(fmt, name, min, max);
      break;
    case eInvalidChoice:
      if (0 == treat_as_string) fmt = "Parameter %s: invalid choice %s; choose from %s.\n";
      else fmt = "Parameter %s: invalid choice \"%s\"; choose from \"%s\".\n";
      ape_msg_error(fmt, name, value, min);
      break;
    case eInvalidName:
      if ('\0' == *name)
        ape_msg_error("Ape: blank parameter name in parameter file.\n");
      else
        ape_msg_error("Ape: invalid parameter name \"%s\" in parameter file.\n", name);
      break;
    case eFileNotAccessible:
      if ('\0' != type[1])
        ape_msg_error("Parameter %s: file \"%s\" is not accessible in mode %s.\n", name, value, type + 1);
      else
        ape_msg_error("Parameter %s: file \"%s\" is not accessible.\n", name, value);
      break;
    case eNan:
       /* Note that by the time ape_par_check calls this function, this status should have
          been converted to a different status based on context, e.g. eInvalidRange. Thus
          this message SHOULD NOT ever be displayed. If this message is displayed, it is necessary
          to check the logic of ape_par_check. */
       ape_msg_debug("Parameter %s: value \"%s\" is interpreted as infinite/NaN.\n", name, value);  
      break;
    case eUndefinedValue:
       /* Note that by the time ape_par_check calls this function, this status should have
          been converted to a different status based on context, e.g. eInvalidRange. Thus
          this message SHOULD NOT ever be displayed. If this message is displayed, it is necessary
          to check the logic of ape_par_check. */
       ape_msg_debug("Parameter %s: value \"%s\" is interpreted as undefined.\n", name, value);  
      break;
    default:
      ape_msg_error("Parameter %s: unexpected Ape error code %d.\n", name, status);
      break;
  }
}

int ape_par_prompt(ApePar * par) {
  int status = eOK;
  char * prompt = 0;
  char * new_prompt = 0;
  char * prompt_base = 0;
  char * value = 0;
  char * min = 0;
  char * max = 0;
  int range_status = eOK;
  int prompt_style = eDefaultPrompt;
  char reprompt = 0;

  /* Check arguments. */
  if (0 == par) status = eNullPointer;

  /* Get the prompt style for the parameter. */
  if (eOK == status) status = ape_par_get_prompt_style(par, &prompt_style);

  /* Suppress prompt if prompt style or default prompt style is "no prompt". */
  if (eOK == status && 0 != (eNoPrompt & (s_default_prompt_style | par->prompt_style))) return status;

  if (eOK == status) {
    /* Get the string containing the minimum value for the parameter. */
    range_status = ape_par_get_min_string(par, &min);
    if (eOK == range_status || eRangeEnum == range_status) status = eOK;
  }

  if (eOK == status) {
    /* Get the string containing the maximum value for the parameter. */
    status = ape_par_get_field(par, eMax, &max);
  }

  if (eOK == status) {
    /* Get the prompt from the parameter. */
    status = ape_par_get_field(par, ePrompt, &prompt);
  }

  if (eOK == status && ('\0' != *min || '\0' != *max)) {
    /* Append opening parenthesis. */
    status = ape_util_cat_string(prompt, " (", &new_prompt);
    if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
  }

  if (eOK == status && '\0' != *min) {
    /* Append minimum. */
    status = ape_util_cat_string(prompt, min, &new_prompt);
    if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
  }

  if (eOK == status && eOK == range_status && ('\0' != *min || '\0' != *max)) {
    /* Append colon. */
    status = ape_util_cat_string(prompt, ":", &new_prompt);
    if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
  }

  if (eOK == status && eOK == range_status && '\0' != *max) {
    /* Append maximum. */
    status = ape_util_cat_string(prompt, max, &new_prompt);
    if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
  }

  if (eOK == status && ('\0' != *min || '\0' != *max)) {
    /* Append closing parenthesis. */
    status = ape_util_cat_string(prompt, ") ", &new_prompt);
    if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
  }

  if (eOK == status) {
    /* Append [ to the prompt. */
    status = ape_util_cat_string(prompt, "[", &new_prompt);
    if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
  }

  if (eOK == status) {
    /* Save a copy of the prompt as constructed thus-far. */
    status = ape_util_copy_string(prompt, &prompt_base);
  }

  do {
    /* Remain optimistic; assume reprompting will not be necessary. */
    reprompt = 0;

    if (eOK == status) {
      /* Get the string containing the value for the parameter. */
      status = ape_par_get_field(par, eValue, &value);
    }

    if (eOK == status) {
      /* Append value to the prompt base. */
      status = ape_util_cat_string(prompt_base, value, &new_prompt);
      if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
    }

    if (eOK == status) {
      /* Append final ] and space to the prompt. */
      status = ape_util_cat_string(prompt, "] ", &new_prompt);
      if (eOK == status) { free(prompt); prompt = new_prompt; new_prompt = 0; }
    }

    if (eOK == status) {
      /* Get the value from the user, using readline to make it more user-friendly. */
      char * new_value = 0;
      if (0 != custom_get_text) {
        char * name = 0;
        status = ape_par_get_field(par, eName, &name);
        if (eOK == status) {
          status = custom_get_text(prompt, name, &new_value);
        }
        free(name);
      } else {
        status = ape_util_get_text(prompt, &new_value);
      }

      if (eOK == status) {
        /* A non-empty string is assigned to the parameter. Note that in this context, empty means truly
           empty. If new_value contains only whitespace that means blank out the new_value. */
        if ('\0' != *new_value)
          status = ape_par_set_value_string(par, new_value);

        /* Check parameter for validity. */
        if (eOK == status) {
          /* Check whether parameter prompting was successful. */
          status = ape_par_check(par, 1);

          /* Problems with the user's input cause reprompting. */
          if (eOK != status) {
            if (eOK == is_recoverable(par, status)) {
              reprompt = 1;
              status = eOK;
            }
          }
        }
      }
      free(new_value); new_value = 0;
    }
    free(value); value = 0;
    free(prompt); prompt = 0;
  } while (0 != reprompt);

  free(prompt_base); prompt_base = 0;
  free(max); max = 0;
  free(min); min = 0;

  return status;
}

int ape_par_query(ApePar * par, const char * auto_string) {
  int status = eOK;
  char eff_mode[APE_PAR_MODE_CODE_LEN] = "";

  if (0 == par) status = eNullPointer;

  if (eOK == status) {
    status = ape_par_get_eff_mode(par, auto_string, eff_mode);
  }

  if (eOK == status && 'q' == eff_mode[0]) {
    status = ape_par_prompt(par);
  }

  return status;
}

static void par_error_msg(const ApePar * par, char * first_error) {
  if ('\0' != *first_error) {
    char * line = 0;
    *first_error = '\0';
    if (0 != (line = ape_par_get_line(par))) {
      ape_msg_debug("\nThe following parameter:\n    \"%s\"\nhas the following error(s):\n", line);
    }
    free(line); line = 0;
  }
}

static void par_note_msg(const ApePar * par) {
  char * line = 0;
  if (0 != (line = ape_par_get_line(par))) {
    ape_msg_debug("\nNote: in the following parameter:\n    \"%s\"\n", line);
  }
  free(line); line = 0;
}

/* TODO Is this struct really useful with the current breakdown of responsibilities? */
/** \brief The RangedParValue structure is a purely internal struct used to pass unpacked
    parameter value, min and max to the check_range function.
*/
typedef struct RangedParValue {
  const ApePar * par;
  const char * value_string;
  int * value_status;
  const char * min_string;
  int * min_status;
  const char * max_string;
  int * max_status;
  char ** enum_range;
  const char * type_string;
  char type;
} RangedParValue;

/** \brief Helper function for check_range function. This converts a parameter string to
    the appropriate C primitive type, and returns a status code indicating whether the
    string may be converted to the appropriate C primitive type without error.
    \param type Code giving the type of this parameter.
    \param string_value The string being evaluated.
*/
static int check_as(char type, const char * string_value) {
  int status = eOK;

  switch (type) {
    case 'b': {
        char value = 0;
        status = ape_util_s2b(string_value, &value);
      }
      break;
    case 'd':
    case 'r': {
        double value = 0.;
        status = ape_util_s2d(string_value, &value);
      }
      break;
    case 'f':
    case 'g':
    case 's':
      break;
    case 'i': {
        long value = 0;
        status = ape_util_s2l(string_value, &value);
      }
      break;
    default:
      status = eUnknownType;
      break;
  }

  return status;
}

/** \brief Helper function for check_range function. This converts two parameter strings to
    the appropriate C primitive type, and returns -1 if the first is < the second, 0 if they are ==
    and +1 if the first is > the second.
    \param type Code giving the type of this parameter.
    \param s1 The first string being compared.
    \param s2 The second string being compared.
    \param result The result of the comparison.
*/
static int compare_as(char type, const char * s1, const char * s2, int * result) {
  int status = eOK;

  *result = 0;

  switch (type) {
    case 'b': {
        char value1 = 0;
        char value2 = 0;
        status = ape_util_s2b(s1, &value1);
        if (eOK == status) {
          status = ape_util_s2b(s2, &value2);
        }
        if (eOK == status) {
          *result = value1 - value2;
        }
      }
      break;
    case 'd':
    case 'r': {
        double value1 = 0;
        double value2 = 0;
        status = ape_util_s2d(s1, &value1);
        if (eOK == status) {
          status = ape_util_s2d(s2, &value2);
        }
        if (eOK == status) {
          if (value1 == value2 ||
            ((value1 != 0.) && (DBL_EPSILON >= fabs((value2 - value1) / value1))) ||
            ((value2 != 0.) && (DBL_EPSILON >= fabs((value1 - value2) / value2)))
          ) *result = 0;
          else if (value1 < value2) *result = -1;
          else *result = +1;
        }
      }
      break;
    case 'f':
    case 'g':
    case 's': {
        int comparison = ape_util_strcmp(s1, s2, 1);
        *result = comparison < 0 ? -1 : (comparison > 0 ? 1 : 0);
      }
      break;
    case 'i': {
        long value1 = 0;
        long value2 = 0;
        status = ape_util_s2l(s1, &value1);
        if (eOK == status) {
          status = ape_util_s2l(s2, &value2);
        }
        if (eOK == status) {
          int comparison = value1 - value2;
          *result = comparison < 0 ? -1 : (comparison > 0 ? 1 : 0);
        }
      }
      break;
    default:
      status = eUnknownType;
      break;
  }

  return status;
}

/** \brief Helper function for ape_par_check. This checks whether the value for the parameter
    can be cleanly converted to the correct primitive C type for the parameter type. In addition,
    the value is checked using min/max parameter values.

    This is a long complex function. The first section checks value, min and max to see if they
    may safely be converted/compared at all, both for normal range and enumerated ranges. Then
    independent checks for the appropriate range (is value >= min, etc.) are performed. Finally,
    the min and max are compared to make sure they are mutually compatible.
    \param par A struct containing the unpacked value, min, max etc. as strings.
    \param first_error Flag indicating whether an error is the first one, used just for cosmetic output.
    \param check_value Flag indicating that the parameter value should be checked.
*/
static int check_range(RangedParValue * par, char * first_error, char check_value) {
  int status = eOK; /* Global status, used to report errors in the test process itself. */
  int * range_status = 0; /* Array of status values for each element in enumerated value list. */
  size_t range_len = 0; /* Number of elements in enumerated value list. */
  int min_status = eOK; /* Used just for checking value >= min. */
  int max_status = eOK; /* Used just for checking value <= max. */
  int enum_status = eOK; /* Used just for checking whether value is one of the enumerated values. */

  if (0 != check_value) {
    /* Check whether value field may be converted to the correct type for this parameter. */
    if (eOK == status && eOK == *par->value_status) {
      *par->value_status = check_as(par->type, par->value_string);
      if (eUndefinedValue == *par->value_status) {
        par_note_msg(par->par);
        ape_msg_debug("the value \"%s\" is considered an \"undefined\" parameter of type \"%s\".\n",
          par->value_string, par->type_string);
        *par->value_status = eOK;
        check_value = 0; /* Disable further range checking below. */
      } else if (eNan == *par->value_status) {
        par_note_msg(par->par);
        ape_msg_debug("the value \"%s\" is considered an NaN parameter of type \"%s\".\n",
          par->value_string, par->type_string);
        *par->value_status = eOK;
        check_value = 0; /* Disable further range checking below. */
      } else if (eOK != *par->value_status) {
        par_error_msg(par->par, first_error);
        ape_msg_debug("  o Problem converting value field \"%s\" to parameter type \"%s\".\n",
          par->value_string, par->type_string);
        /* Note that check_as status is returned directly. */
      }
    }
  }

  /* Check whether minimum field may be converted to the correct type for this parameter. */
  if (eOK == status && eOK == *par->min_status) {
    *par->min_status = check_as(par->type, par->min_string);
    if (eOK != *par->min_status) {
      par_error_msg(par->par, first_error);
      if (eUndefinedValue == *par->min_status) {
        ape_msg_debug("  o Minimum field \"%s\" may not contain strings meaning \"undefined\".\n", par->min_string);
        /* Change status to reflect that the range is invalid. */
        *par->min_status = eInvalidRange;
      } else if (eNan == *par->min_status) {
        ape_msg_debug("  o Minimum field \"%s\" may not contain strings meaning \"NaN\".\n", par->min_string);
        /* Change status to reflect that the range is invalid. */
        *par->min_status = eInvalidRange;
      } else {
        ape_msg_debug("  o Problem converting minimum field \"%s\" to parameter type \"%s\".\n", par->min_string, par->type_string);
        /* Change status: not so much detail about the min field. */
        *par->min_status = eMinConversionError;
      }
    }
  }

  /* Check whether maximum field may be converted to the correct type for this parameter. */
  if (eOK == status && eOK == *par->max_status) {
    *par->max_status = check_as(par->type, par->max_string);
    if (eOK != *par->max_status) {
      par_error_msg(par->par, first_error);
      if (eUndefinedValue == *par->max_status) {
        ape_msg_debug("  o Maximum field \"%s\" may not contain strings meaning \"undefined\".\n", par->max_string);
        /* Change status to reflect that the range is invalid. */
        *par->max_status = eInvalidRange;
      } else if (eNan == *par->max_status) {
        ape_msg_debug("  o Maximum field \"%s\" may not contain strings meaning \"NaN\".\n", par->max_string);
        /* Change status to reflect that the range is invalid. */
        *par->max_status = eInvalidRange;
      } else {
        ape_msg_debug("  o Problem converting maximum field \"%s\" to parameter type \"%s\".\n", par->max_string, par->type_string);
        /* Change status: not so much detail about the max field. */
        *par->max_status = eMaxConversionError;
      }
    }
  }

  /* Convert enumerated range, if present. */
  if (eOK == status && eOK == *par->value_status && eRangeEnum == *par->min_status) {
    /* Compute one more than the number of enumerated choices. */
    while (0 != par->enum_range[range_len]) ++range_len;

    if (0 != range_len) {
      /* Allocate a buffer to hold a flag indicating status for each of the enumerated choices. */
      range_status = (int *) calloc(range_len, sizeof(int));
      if (0 == range_status) status = eDynAllocFailed;

      if (eOK == status) {
        /* Convert each enumerated value into the appropriate type value. */
        size_t idx = 0;
        for (idx = 0; idx != range_len; ++idx) {
          range_status[idx] = check_as(par->type, par->enum_range[idx]);
        }
      }
    }
  }

  if (0 != check_value) {
    /* Checks of value against min, max and/or enumerated ranges are done separately, and results
       are stored separately. This is so that *all* current errors can be reported, not just the
       first one. The separate status values from each test are then combined to give the status
       of the value, min and max fields. */
    /* Check whether value is >= minimum. */
    if (eOK == status && eOK == *par->value_status && eOK == *par->min_status) {
      int comparison = 0;
      min_status = compare_as(par->type, par->value_string, par->min_string, &comparison);
      if (eOK == min_status && 0 > comparison) {
        par_error_msg(par->par, first_error);
        ape_msg_debug("  o Value field \"%s\" < minimum field value \"%s\".\n", par->value_string, par->min_string);
        min_status = eValueBelowMin;
      }
    }

    /* Check whether value is <= maximum. */
    if (eOK == status && eOK == *par->value_status && eOK == *par->max_status) {
      int comparison = 0;
      max_status = compare_as(par->type, par->value_string, par->max_string, &comparison);
      if (eOK == max_status && 0 < comparison) {
        par_error_msg(par->par, first_error);
        ape_msg_debug("  o Value field \"%s\" > maximum field value \"%s\".\n", par->value_string, par->max_string);
        max_status = eValueAboveMax;
      }
    }

    /* Check whether value is one of the enumerated range values. */
    /* TODO Can this be changed to use ape_util_find_string? */
    if (eOK == status && eOK == *par->value_status && eRangeEnum == *par->min_status) {
      size_t idx = 0;
      char found = 0;
      for (idx = 0; 0 == found && idx != range_len; ++idx) {
        int comparison = 0;
        if (eOK != range_status[idx]) continue;
        range_status[idx] = compare_as(par->type, par->value_string, par->enum_range[idx], &comparison);
        if (eOK == range_status[idx] && 0 == comparison) found = 1;
      }
      if (0 == found) {
        par_error_msg(par->par, first_error);
        ape_msg_debug("  o Value field \"%s\" does not match one of the values \"%s\".\n", par->value_string, par->min_string);
        enum_status = eInvalidChoice;
      }
    }

    /* Use the separate results from the three range checks above to set the value's status. */
    if (eOK == status && eOK == *par->value_status) {
      if (eOK != min_status) *par->value_status = min_status;
      else if (eOK != max_status) *par->value_status = max_status;
      else if (eOK != enum_status) *par->value_status = enum_status;
    }

    /* Check for file access using the remainder of the field string as the access type to check. */
    if (eOK == status && eOK == *par->value_status && '\0' != *par->value_string && 'f' == par->type) {
      *par->value_status = ape_util_check_file_access(par->value_string, par->type_string + 1);
      if (eOK != *par->value_status)
        ape_msg_debug("  o File \"%s\" is not accessible in mode \"%s\".\n", par->value_string, par->type_string + 1);
    }

  }

  /* Check whether range itself is valid. */
  if (eOK == status && eOK == *par->min_status && eOK == *par->max_status) {
    int comparison = 0;
    *par->min_status = compare_as(par->type, par->min_string, par->max_string, &comparison);
    if (eOK == *par->min_status && 0 < comparison) {
      par_error_msg(par->par, first_error);
      ape_msg_debug("  o Minimum field value \"%s\" > maximum field value \"%s\".\n", par->min_string, par->max_string);
      *par->min_status = eInvalidRange;
      *par->max_status = eInvalidRange;
    }
  }

  /* Check whether enumerated range is valid. */
  if (eOK == status && eRangeEnum == *par->min_status) {
    size_t idx = 0;
    for (idx = 0; idx != range_len; ++idx) {
      if (eOK != range_status[idx]) {
        par_error_msg(par->par, first_error);
        if (eUndefinedValue == range_status[idx]) {
          ape_msg_debug("  o Enumerated field \"%s\" may not contain strings meaning \"undefined\".\n", par->enum_range[idx]);
        } else if (eNan == range_status[idx]) {
          ape_msg_debug("  o Enumerated field \"%s\" may not contain strings meaning \"NaN\".\n", par->enum_range[idx]);
        } else {
          ape_msg_debug("  o Invalid enumerated field value: \"%s\".\n", par->enum_range[idx]);
        }
        *par->min_status = eInvalidRange;
      }
    }
    if (eOK == *par->max_status) {
      par_error_msg(par->par, first_error);
      ape_msg_debug("  o Minimum field \"%s\" defines an enumerated range, but maximum field \"%s\" is also defined.\n",
        par->min_string, par->max_string);
      *par->max_status = eInvalidRange;
    }
  }

  /* The flag which indicates minimum field contains an enumerated range is not strictly an error, so reset it
     now that the enumerated nature has been handled. */
  if (eRangeEnum == *par->min_status) *par->min_status = eOK;

  /* Clean up. */
  free(range_status); range_status = 0;

  return status;
}

int ape_par_check(const ApePar * par, char check_value) {
  /* Global status indicating whether the check itself is working (no memory problems etc.). */
  int status = eOK;
  size_t num_fields = 0;
  size_t idx = 0;
  /* Individual status for each field, indicating problems with the format of the parameter, but not memory issues etc. */
  int field_status[eEndOfField] = { eFieldNotFound, eFieldNotFound, eFieldNotFound, eFieldNotFound, eFieldNotFound,
    eFieldNotFound, eFieldNotFound };
  char ** field = 0;
  char ** enum_range = 0;
  char first_error = 1;
  char blank_line = 0;

  /* Check arguments. */
  if (0 == par) {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Create pointers to hold each field. */
    field = (char **) calloc(eEndOfField + 1, sizeof(char *));
    if (0 == field) status = eDynAllocFailed;
  }

  if (eOK == status) {
    /* Get number of fields in this parameter. */
    num_fields = ape_list_get_size(par->field);

    /* Iterate over all fields which are present, getting their values. */
    for (idx = 0; eOK == status && idx != num_fields && idx != eEndOfField; ++idx) {
      switch (idx) {
        case eName:
        case eValue:
        case eMax:
        case ePrompt:
          status = ape_par_get_field(par, idx, field + idx);
          if (eOK == status && '\0' != *field[idx]) field_status[idx] = status;
          break;
        case eType: {
            /* Get the type code; interpret status specially and handle case of unknown type. */
            char type[APE_PAR_TYPE_CODE_LEN];
            status = ape_par_get_type(par, type);
            if (eOK == status) {
              status = ape_util_copy_string(type, field + idx);
              field_status[idx] = status;
            } else if (eUnknownType == status) {
              /* Type is unknown, but still get the content of the field so it may be reported. */
              status = ape_par_get_field(par, eType, field + idx);
              field_status[idx] = eUnknownType;
            } else {
              field_status[idx] = status;
            }
          }
          break;
        case eMode: {
            /* Get the mode code; interpret status specially and handle case of unknown mode. */
            char mode[APE_PAR_MODE_CODE_LEN];
            status = ape_par_get_mode(par, mode);
            if (eOK == status) {
              status = ape_util_copy_string(mode, field + idx);
              field_status[idx] = status;
            } else if (eUnknownMode == status) {
              /* Mode is unknown, but still get the content of the field so it may be reported. */
              status = ape_par_get_field(par, eMode, field + idx);
              field_status[idx] = eUnknownMode;
            } else {
              field_status[idx] = status;
            }
          }
          break;
        case eMin:
          /* Get the minimum value; this will return a non-zero status if the field has invalid content. */
          status = ape_par_get_min_string(par, field + idx);
          if (eOK == status && '\0' != *field[idx]) {
            field_status[idx] = status;
          } else if (eRangeEnum == status) {
            field_status[idx] = status;
            status = ape_par_get_enum_string(par, &enum_range);
          }
          break;
        default:
          break;
      }
    }
  }

  /* Describe all errors in the format of this parameter. */
  /* A blank name field is an error unless there are no other fields, indicating a blank line. */
  if (eOK == status && (0 == field[eName] || '\0' == *field[eName])) {
    if (1 == num_fields) {
      blank_line = 1;
    } else {
      field_status[eName] = eInvalidName;
      par_error_msg(par, &first_error);
      ape_msg_debug("  o Name field is blank.\n");
    }
  }

  /* Report problems with the type field. */
  if (eOK == status && eUnknownType == field_status[eType]) {
    static char blank[] = "blank";
    const char * type_string = field[eType];
    if (0 == type_string) type_string = blank;
    par_error_msg(par, &first_error);
    ape_msg_debug("  o Type field \"%s\" is invalid.\n", type_string);
  }

  /* Report problems with the mode field. */
  if (eOK == status && eUnknownMode == field_status[eMode]) {
    static char blank[] = "blank";
    const char * mode_string = field[eMode];
    if (0 == mode_string) mode_string = blank;
    par_error_msg(par, &first_error);
    ape_msg_debug("  o Mode field \"%s\" is invalid.\n", mode_string);
  }

  /* Report problems with the value; most of the checking requires the type to be known. */
  if (eOK == status && eOK == field_status[eType]) {
    /* Pack all relevant information about value, min, max into a structure. */
    RangedParValue par_value;
    par_value.par = par;
    par_value.value_string = field[eValue];
    par_value.value_status = field_status + eValue;
    par_value.min_string = field[eMin];
    par_value.min_status = field_status + eMin;
    par_value.max_string = field[eMax];
    par_value.max_status = field_status + eMax;
    par_value.enum_range = enum_range;
    par_value.type_string = field[eType];
    par_value.type = *field[eType];

    /* Check types of all parts which have a type. */
    status = check_range(&par_value, &first_error, check_value);
  }

  /* Report unexpected numbers of fields, except in the case of a blank line. */
  if (eOK == status && eEndOfField != num_fields && 0 == blank_line) {
    par_error_msg(par, &first_error);
    ape_msg_debug("  o Parameter has %u field%s, not %u as expected.\n", num_fields, num_fields != 1 ? "s" : "", eEndOfField);
  }

  /* Determine final status to report, which should be the status of the most serious error encountered. In
     determining the relative seriousness of problems, the most critical aspect is that errors having to do with
     the value are considere the least serious. This is so that the caller can tell that a status associated with
     a bad value necessarily implies nothing *worse* than that happened. Bad values are expected to be the most
     common and most easily correctable (e.g. repeat a prompt for value out of range). In descending order of
     seriousness, possible problems are:
       Invalid name (null name string)
       Invalid type
       Invalid mode
       Min or max cannot be converted to correct C type for this parameter type
       Min > max
       Min is enum but max is also given
       Too few fields present
       Too many fields present
       Value cannot be converted to correct C type for this parameter type
       Value out of range
       File not found (file type only)
  */
  if (eOK == status) {
    if (eOK != field_status[eName] && eFieldNotFound != field_status[eName]) status = field_status[eName];
    else if (eOK != field_status[eType] && eFieldNotFound != field_status[eType]) status = field_status[eType];
    else if (eOK != field_status[eMode] && eFieldNotFound != field_status[eMode]) status = field_status[eMode];
    else if (eOK != field_status[eMax] && eFieldNotFound != field_status[eMax]) status = field_status[eMax];
    else if (eOK != field_status[eMin] && eFieldNotFound != field_status[eMin]) status = field_status[eMin];
    else if (eEndOfField > num_fields && 0 == blank_line) status = eFieldNotFound;
    else if (eEndOfField < num_fields && 0 == blank_line) status = eTooManyFields;
    else if (eOK != field_status[eValue] && eFieldNotFound != field_status[eValue]) status = field_status[eValue];
  }

  if (0 != field) describe_status(field, field_status, status);

  /* TODO Only print the following extra line if actual final status != 0. */
/*  ape_msg_debug("\n"); */

  /* Clean up. */
  ape_util_free_string_array(enum_range); enum_range = 0;
  ape_util_free_string_array(field); field = 0;

  return status;
}

int ape_par_flag_cmd_line(ApePar * par) {
  int status = eOK;

  if (0 == par) status = eNullPointer;

  if (eOK == status) {
    par->cmd_line = 1;
  }

  return status;
}

int ape_par_redirect_prompt_stream(FILE * stream) {
  return ape_util_set_prompt_stream(stream);
}

int ape_par_register_get_text(ApeParGetTextFunc client_get_text) {
  custom_get_text = client_get_text;
  return eOK;
}

static int test_register_get_text(const char * prompt, const char * name, char ** new_value) {
  *new_value = "Hello World";
  return eInputFailure;
}

void ape_par_test(void) {
  { const char * correct_field[] = { "par0", "s", "a", "value \"0\"", "", "", "Prompt, for parameter par0", 0 };
    test_one_par(" par0    , s , a,  \"value \\\"0\\\"\", ,, \"Prompt, for parameter par0\"", correct_field);
  }
  { const char * correct_field[] = { "par0", "s", "a", "value \"0\"", "", "", 0 };
    test_one_par(" par0    , s , a, \"value \\\"0\\\"\", ,    ", correct_field);
  }
  { const char * correct_field[] = { "", 0 };
    test_one_par("", correct_field);
  }
  { const char * correct_field[] = { "", 0 };
    test_one_par("\"\"", correct_field);
  }
  { const char * correct_field[] = { "", 0 };
    test_one_par(" \t", correct_field);
  }
  { const char * correct_field[] = { "", "", 0 };
    test_one_par(",", correct_field);
  }
  { const char * correct_field[] = { "", "", "", 0 };
    test_one_par(",,", correct_field);
  }
  { const char * correct_field[] = { "", "", "", 0 };
    test_one_par(" ,,", correct_field);
  }
  { const char * correct_field[] = { "", "", "", 0 };
    test_one_par(", ,", correct_field);
  }
  { const char * correct_field[] = { "", "", "", 0 };
    test_one_par(",, ", correct_field);
  }
  { const char * correct_field[] = { "\"value,", 0 };
    test_one_par("\"value, ", correct_field);
  }
  { const char * correct_field[] = { "\"value", "", 0 };
    test_one_par("\\\"value, ", correct_field);
  }
  { const char * correct_field[] = { "value", 0 };
    test_one_par(" value", correct_field);
  }
  { const char * correct_field[] = { "\"value", 0 };
    test_one_par(" \"value", correct_field);
  }
  { const char * correct_field[] = { "value0 \" value1 \" value2", 0 };
    test_one_par("value0 \" value1 \" value2", correct_field);
  }
  { const char * correct_field[] = { "value0\" \" value1 \"\"value2", 0};
    test_one_par("\"value0\" \" value1 \"\"value2\"", correct_field);
  }
  { const char * correct_field[] = { "par0", "s", "a", "value0 value1 value2", "", "", 0 };
    test_one_par(" par0    , s , a, value0 value1 value2, ,    ", correct_field);
  }
  { const char * correct_field[] = { ", ", "s,", ",", ",value0 value1 value2", "", "min0 , min1", 0 };
    test_one_par(" \", \" , \"s,\" , \",\", \",value0 value1 value2\", , \"min0 , min1\"   ", correct_field);
  }
  { const char * correct_field[] = { "\",,,,,,", 0 };
    test_one_par(" \",,,,,,", correct_field);
  }
  /* Make sure single and double quotes both work. */
  { const char * correct_field[] = { "par", "s", "h", "the value", "the min", "the max", "the prompt", 0 };
    test_one_par("par, s, h, the value, \"the min\", 'the max', the prompt", correct_field);
  }
  /* Make sure comments are always handled correctly. */
  { const char * correct_field[] = {
      "par", "s", "h", "value with # embedded", "min with # embedded", "max with # embedded", "prompt", 0
    };
    test_one_par("par, s, h, value with # embedded, min with # embedded, max with # embedded, prompt # with comment",
      correct_field);
  }
  { const char * correct_field[] = {
      "par", "s", "h", "value with # embedded", "min with # embedded", "max with # embedded", "prompt with # embedded", 0
    };
    test_one_par("par, s, h, value with # embedded, min with # embedded, max with # embedded, prompt with \\# embedded",
      correct_field);
  }
  { const char * correct_field[] = {
      "par", "s", "h", "value with # embedded", "min with # embedded", "max with # embedded", "prompt with # embedded", 0
    };
    test_one_par("par, s, h, value with # embedded, min with # embedded, max with # embedded, \"prompt with # embedded\"",
      correct_field);
  }
  { const char * correct_field[] = {
      "par", "s", "h", "value with # embedded", "min with # embedded", "max with # embedded", "prompt with # embedded", 0
    };
    test_one_par("par, s, h, value with # embedded, min with # embedded, max with # embedded, 'prompt with # embedded'",
      correct_field);
  }

  { const char * correct_field[] = { "par 0", "s", "a", "value", "a", "z",
      "Prompt with 'single', \"double\", even unmatched \" quotes", 0 };
    test_one_par(
      " 'par 0'    , \"s\" , 'a', 'value', \"a\", 'z', 'Prompt with \\'single\\', \"double\", even unmatched \" quotes'",
      correct_field
    );
  }
  /* Test whether a single space is preserved when a parameter is created:*/
  { const char * correct_field[] = { "par", "s", "h", " ", "", "", "prompt", 0 };
    test_one_par("par, s, h, \" \", , , prompt",correct_field);
  }

  /* Test ape_par_get_double. */
  { ApePar * par = 0;
    int status = ape_par_create(&par, "double_par, r, a, 3., 1., 5., Unquoted prompt");
    if (0 == par) {
      ape_test_failed("When about to test ape_par_get_double, unexpected problem creating a test parameter object.\n");
    } else {
      double result = 1.;

      /* Test ape_par_get_double with one or more null pointers as inputs. */
      status = ape_par_get_double(0, 0);
      if (eNullPointer != status)
        ape_test_failed("ape_par_get_double(0, 0) returned status %d, not %d as expected.\n", status, eNullPointer);

      status = ape_par_get_double(0, &result);
      ape_test_cmp_double("ape_par_get_double(0, &result)", result, 0., status, eNullPointer);

      status = ape_par_get_double(par, 0);
      ape_test_cmp_double("ape_par_get_double(par, 0)", result, result, status, eNullPointer);

      /* Test valid conversion of field. */
      result = 0;
      status = ape_par_get_double(par, &result);
      ape_test_cmp_double("ape_par_get_double(par, 0)", result, 3., status, eOK);
    }

    /* Clean up. */
    ape_par_destroy(par);
  }

  /* Test ape_par_get_short. */
  { ApePar * par = 0;
    int status = ape_par_create(&par, "short_par, i, a, 40000, , , Unquoted prompt");
    if (0 == par) {
      ape_test_failed("When about to test ape_par_get_short, unexpected problem creating a test parameter object.\n");
    } else {
      short result = SHRT_MAX;

      /* Test ape_par_get_short with one or more null pointers as inputs. */
      status = ape_par_get_short(0, 0);
      if (eNullPointer != status)
        ape_test_failed("ape_par_get_short(0, 0) returned status %d, not %d as expected.\n", status, eNullPointer);

      status = ape_par_get_short(0, &result);
      ape_test_cmp_long("ape_par_get_short(0, &result)", result, 0., status, eNullPointer);

      status = ape_par_get_short(par, 0);
      ape_test_cmp_long("ape_par_get_short(par, 0)", result, result, status, eNullPointer);

      /* Test valid conversion of field, but resulting in an overflow error. */
      result = 0;
      status = ape_par_get_short(par, &result);
      ape_test_cmp_long("ape_par_get_short(par, 0)", result, SHRT_MAX, status, eOverflow);
    }

    /* Clean up. */
    ape_par_destroy(par);
  }

  /* Test ape_par_get_field. */
  { ApePar * par = 0;
    int status = ape_par_create(&par,
      "correct par name, rF, a, \"correct $ENV{NON_EXISTENT_ENV_VAR} \\$ENV{value}\", '$ENV{QUOTED_ENV_VAR}', , Unquoted prompt");
    if (0 == par) {
      ape_test_failed("When about to test ape_par_get_field, unexpected problem creating a test parameter object.\n");
    } else {
      const char * expected[] = { "correct par name", "rF", "a", "correct  $ENV{value}", "$ENV{QUOTED_ENV_VAR}", "",
        "Unquoted prompt" };
      char * result = 0;
      char incorrect[] = "incorrect";
      ParFieldId idx;

      for (idx = eName; idx != eEndOfField; ++idx) {
        /* Test ape_par_get_field with one or more null pointers as inputs. */
        status = ape_par_get_field(0, idx, 0);
        if (eNullPointer != status)
          ape_test_failed("ape_par_get_field(0, %d, 0) returned status %d, not %d as expected.\n", idx, status, eNullPointer);

        result = incorrect;
        status = ape_par_get_field(0, idx, &result);
        if (eNullPointer != status)
          ape_test_failed("ape_par_get_field(0, %d, &result) returned status %d, not %d as expected.\n", idx, status,
            eNullPointer);
        if (0 != result)
          ape_test_failed("ape_par_get_field(0, %d, &result) gave result \"%s\", not 0 as expected.\n", idx, result);

        status = ape_par_get_field(par, idx, 0);
        if (eNullPointer != status)
          ape_test_failed("ape_par_get_field(par, %d, 0) returned status %d, not %d as expected.\n", idx, status, eNullPointer);

        /* Test valid conversion of fields. */
        result = 0;
        status = ape_par_get_field(par, idx, &result);
        if (eOK != status)
          ape_test_failed("ape_par_get_field(par, %d, &result) returned status %d, not %d as expected.\n", idx, status, eOK);
        if (0 != strcmp(expected[idx], result))
          ape_test_failed("ape_par_get_field(par, %d, &result) returned string \"%s\", not \"%s\" as expected.\n",
            idx, result, expected[idx]);
        free(result); result = incorrect;
      }
    }

    /* Clean up. */
    ape_par_destroy(par);
  }
  { ApePar * par = 0;
    int status = ape_par_create(&par, "sa, s, a, opTIon1, OPTion1|OPTION2|option3, , sa enumerated");
    if (0 == par) {
      ape_test_failed("Unable to set up to test ape_par_get_string_case for enumerated string, status is %d.\n", status);
    } else {
      const char * expected = "OPTION1";
      char * result = 0;
      status = ape_par_get_string(par, &result);
      ape_test_cmp_string("ape_par_get_string called for enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "opTIon1";
      status = ape_par_get_string_case(par, &result, eDefaultCase);
      ape_test_cmp_string("ape_par_get_string_case called with eDefaultCase for enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "option1";
      status = ape_par_get_string_case(par, &result, eLowerCase);
      ape_test_cmp_string("ape_par_get_string_case called with eLowerCase for enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "OPTION1";
      status = ape_par_get_string_case(par, &result, eUpperCase);
      ape_test_cmp_string("ape_par_get_string_case called with eUpperCase for enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "OPTion1";
      status = ape_par_get_string_case(par, &result, eEnumCase);
      ape_test_cmp_string("ape_par_get_string_case called with eEnumCase for enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "option1";
      status = ape_par_get_string_case(par, &result, eLowerCase | eEnumCase);
      ape_test_cmp_string("ape_par_get_string_case called with eLowerCase | eEnumCase for enumerated range", result,
        expected, status, eOK);
      free(result); result = 0;
    }

    /* Clean up. */
    ape_par_destroy(par);
  }
  { ApePar * par = 0;
    int status = ape_par_create(&par, "sa, s, a, opTIon1, , , sa enumerated");
    if (0 == par) {
      ape_test_failed("Unable to set up to test ape_par_get_string_case for non-enumerated string, status is %d.\n", status);
    } else {
      const char * expected = "opTIon1";
      char * result = 0;
      status = ape_par_get_string(par, &result);
      ape_test_cmp_string("ape_par_get_string called for non-enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "opTIon1";
      status = ape_par_get_string_case(par, &result, eDefaultCase);
      ape_test_cmp_string("ape_par_get_string_case called with eDefaultCase for non-enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "option1";
      status = ape_par_get_string_case(par, &result, eLowerCase);
      ape_test_cmp_string("ape_par_get_string_case called with eLowerCase for non-enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "OPTION1";
      status = ape_par_get_string_case(par, &result, eUpperCase);
      ape_test_cmp_string("ape_par_get_string_case called with eUpperCase for non-enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "opTIon1";
      status = ape_par_get_string_case(par, &result, eEnumCase);
      ape_test_cmp_string("ape_par_get_string_case called with eEnumCase for non-enumerated range", result, expected, status, eOK);
      free(result); result = 0;

      expected = "option1";
      status = ape_par_get_string_case(par, &result, eLowerCase | eEnumCase);
      ape_test_cmp_string("ape_par_get_string_case called with eLowerCase | eEnumCase for enumerated range", result,
        expected, status, eOK);
      free(result); result = 0;
    }

    /* Clean up. */
    ape_par_destroy(par);
  }


  /* Test ape_par_clone. */
  { const char * par_text = "stringpar   , s, a, \"value\", , , \"Enter a string\"";
    ApePar * orig_par = 0;
    ApePar ** clone_par = 0;
    int status = ape_par_create(&orig_par, par_text);
    if (eOK == status) {
      status = ape_par_clone(orig_par, clone_par);
      ape_test_cmp_ptr("ape_par_clone(orig_par, 0)", clone_par, 0, status, eNullPointer);
    } else {
      ape_test_failed("Setup for ape_par_clone(orig_par, 0) test could not create orig_par.\n");
    }
    ape_par_destroy(orig_par);
  }
  { const char * par_text = "stringpar   , s, a, \"value\", , , \"Enter a string\"";
    ApePar * orig_par = 0;
    ApePar * clone_par = 0;
    ApePar * save_par = 0;
    int status = ape_par_create(&clone_par, par_text);
    save_par = clone_par;
    if (eOK == status) {
      status = ape_par_clone(orig_par, &clone_par);
      ape_test_cmp_ptr("ape_par_clone(0, &clone_par)", clone_par, 0, status, eNullPointer);
    } else {
      ape_test_failed("Setup for ape_par_clone(0, &clone_par) test could not create clone_par.\n");
    }
    ape_par_destroy(save_par);
  }
  { const char * par_text = "stringpar   , s, a, \"value\", , , \"Enter a string\"";
    ApePar * orig_par = 0;
    ApePar * clone_par = 0;
    int status = ape_par_create(&orig_par, par_text);
    if (eOK == status) {
      status = ape_par_clone(orig_par, &clone_par);
      if (eOK == status) {
        char * clone_text = ape_par_get_line(clone_par);
        ape_test_cmp_string("after ape_par_clone(orig_par, &clone_par)", clone_text, par_text, status, eOK);
        free(clone_text); clone_text = 0;
      }
    } else {
      ape_test_failed("Setup for ape_par_clone(orig_par, &clone_par) test could not create orig_par.\n");
    }
    ape_par_destroy(clone_par);
    ape_par_destroy(orig_par);
  }

  /* Test ape_par_set_field. */
  { const char * par_text = "stringpar   , s, a, \"value\", , , \"Enter a string\"";
    char buf[1024] = "";
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    if (eOK == status) {
      char * new_value = "new-value";
      status = ape_par_set_field(par, eValue, new_value);
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \"new-value\")", value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_cmp_long("ape_par_set_field(par, eValue, \"new-value\")", 0, 0, status, eOK);
      }
      new_value = "comma, separated,values";
      status = ape_par_set_field(par, eValue, new_value);
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \"comma, separated,values\")", value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_cmp_long("ape_par_set_field(par, eValue, \"comma, separate,values\")", 0, 0, status, eOK);
      }
      new_value = "\"comma, separated,values\"";
      status = ape_par_set_field(par, eValue, new_value);
      new_value = "\"comma, separated,values\"";
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \\\"comma, separated,values\\\")", value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_cmp_long("ape_par_set_field(par, eValue, \\\"comma, separate,values\\\")", 0, 0, status, eOK);
      }
      new_value = "regfilter(\"tpow.reg\", a.p1, a.p2) ? 1 : 0";
      status = ape_par_set_field(par, eValue, new_value);
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \"regfilter(\"tpow.reg\", a.p1, a.p2) ? 1 : 0\")",
          value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_cmp_long("ape_par_set_field(par, eValue, \"regfilter(\"tpow.reg\", a.p1, a.p2) ? 1 : 0\")", 0, 0, status, eOK);
      }

      /* Test a hideous but realistic string. */
      strcat(buf, "infile = bat_b19_data/ftcoco1.fits[col *; RA_OBJ(D)=(((arctan2(      (0.45");
      strcat(buf, "598377618) * (cos((GLAT_TEMPLATE/57.2957795130823)) *  sin((GLON_TEMPLATE");
      strcat(buf, "/57.2957795130823 - (0.574770433)))) + (-0.88998808748) *  sin((GLAT_TEMP");
      strcat(buf, "LATE/57.2957795130823)), cos((GLAT_TEMPLATE/57.2957795130823)) *  cos(   ");
      strcat(buf, "GLON_TEMPLATE/57.2957795130823 - (0.574770433)))  )+(4.9368292465)+12.56");
      strcat(buf, "06143592) % 6.28318530717959)*57.2957795130823); TUNIT#( physical unit o");
      strcat(buf, "f field) = \"deg\"; DEC_OBJ(D)=(arcsin(min((-(-0.88998808748) *  (cos((GLON");
      strcat(buf, "_TEMPLATE/57.2957795130823)) *  sin((GLON_TEMPLATE/57.2957795130823 -   ");
      strcat(buf, "574770433)))) + (0.45598377618) *  sin((GLAT_TEMPLATE/57.2957795130823   ");
      strcat(buf, ",1.0))*57.2957795130823); TUNIT#( physical unit of field) =  \"deg\";]");
      new_value = buf;
      status = ape_par_set_field(par, eValue, new_value);
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \"hideous expression\")",
          value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_cmp_long("ape_par_set_field(par, eValue, \"hideous expression\")", 0, 0, status, eOK);
      }
      new_value = "\\\"backslash quote\\\"";
      status = ape_par_set_field(par, eValue, new_value);
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \"\\\"backslash quote\\\"\")",
          value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_cmp_long("ape_par_set_field(par, eValue, \"regfilter(\"tpow.reg\", a.p1, a.p2) ? 1 : 0\")", 0, 0, status, eOK);
      }

      /* Test single space */
      new_value = " ";
      status = ape_par_set_field(par, eValue, new_value);
      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("after ape_par_set_field(par, eValue, \" \")",
          value, new_value, status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("ape_par_set_field: assignment of string \" \" to value failed with status %d.\n", status);
      }

    } else {
      ape_test_failed("Setup for ape_par_set_field(par, eValue, \"new_value\") test could not create test par.\n");
    }
    ape_par_destroy(par);
  }

  /* Test ape_par_get_type. */
  { const char * par_text = "testpar,    B , a, , , , \"Test parameter\"";
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    if (eOK == status) {
      char type_string[APE_PAR_TYPE_CODE_LEN + 1] = "";

      /* Success cases: identify all known types. */
      status = ape_par_get_type(par, type_string);
      ape_test_cmp_string("ape_par_get_type for type \"B\"", type_string, "b", status, eOK);

      status = ape_par_set_field(par, eType, " d");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"d\"", type_string, "d", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"d\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "F ");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"F \"", type_string, "f", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"F \" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "r F");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"r F\"", type_string, "fr", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"r F\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "F w");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"F w\"", type_string, "fw", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"F w\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, " i ");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \" i \"", type_string, "i", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \" i \" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "R");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"R\"", type_string, "r", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"R\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "s");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"s\"", type_string, "s", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"s\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "fWr");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"fWr\"", type_string, "frw", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"fWr\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "Ef");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"Ef\"", type_string, "fe", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"Ef\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "Efr");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"Efr\"", type_string, "fer", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"Efr\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "wEf");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"wEf\"", type_string, "few", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"wEf\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "wERf");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"wERf\"", type_string, "ferw", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"wERf\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "nf");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"nf\"", type_string, "fn", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"nf\" to par type failed with status %d.\n", status);
      }

      status = ape_par_set_field(par, eType, "wnf");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"wnf\"", type_string, "fnw", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"wnf\" to par type failed with status %d.\n", status);
      }

      /* Error case: blank type string . */
      status = ape_par_set_field(par, eType, "");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"\"", type_string, "", status, eUnknownType);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"rq\" to par type failed with status %d.\n", status);
      }

      /* Error case: include unknown characters. */
      status = ape_par_set_field(par, eType, "rq");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"rq\"", type_string, "r", status, eFormatError);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"rq\" to par type failed with status %d.\n", status);
      }

      /* Error case: include duplicated but known characters. */
      status = ape_par_set_field(par, eType, "bb");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"bb\"", type_string, "b", status, eFormatError);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"bb\" to par type failed with status %d.\n", status);
      }

      /* Error case: include two known but contradicatory characters. */
      status = ape_par_set_field(par, eType, "si");
      if (eOK == status) {
        status = ape_par_get_type(par, type_string);
        ape_test_cmp_string("ape_par_get_type for type \"si\"", type_string, "", status, eUnknownType);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"si\" to par type failed with status %d.\n", status);
      }
    } else {
      ape_test_failed("Setup for ape_par_get_type(par, type_string) failed with status %d.\n", status);
    }
    ape_par_destroy(par);
  }

  /* Test ape_par_get_mode. */
  { const char * par_text = "testpar,    B , a, , , , \"Test parameter\"";
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    if (eOK == status) {
      char mode_string[APE_PAR_MODE_CODE_LEN] = "";

      /* Success cases: identify all known modes. */
      status = ape_par_get_mode(par, mode_string);
      ape_test_cmp_string("ape_par_get_mode for mode \"A\"", mode_string, "a", status, eOK);

      status = ape_par_set_field(par, eMode, " l h");
      if (eOK == status) {
        status = ape_par_get_mode(par, mode_string);
        ape_test_cmp_string("ape_par_get_mode for mode \" l h\"", mode_string, "hl", status, eOK);
      } else {
        ape_test_failed("ape_par_test: assignment of string \" l h\" to par mode failed with status %d.\n", status);
      }

      /* Error case: blank mode string . */
      status = ape_par_set_field(par, eMode, "");
      if (eOK == status) {
        status = ape_par_get_mode(par, mode_string);
        ape_test_cmp_string("ape_par_get_mode for mode \"\"", mode_string, "", status, eUnknownMode);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"rq\" to par mode failed with status %d.\n", status);
      }

      /* Error case: include unknown characters. */
      status = ape_par_set_field(par, eMode, "az");
      if (eOK == status) {
        status = ape_par_get_mode(par, mode_string);
        ape_test_cmp_string("ape_par_get_mode for mode \"az\"", mode_string, "a", status, eFormatError);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"rq\" to par mode failed with status %d.\n", status);
      }

      /* Error case: include duplicated but known characters. */
      status = ape_par_set_field(par, eMode, "qq");
      if (eOK == status) {
        status = ape_par_get_mode(par, mode_string);
        ape_test_cmp_string("ape_par_get_mode for mode \"qq\"", mode_string, "q", status, eFormatError);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"qq\" to par mode failed with status %d.\n", status);
      }

      /* Error case: include two known but contradicatory characters. */
      status = ape_par_set_field(par, eMode, "qh");
      if (eOK == status) {
        status = ape_par_get_mode(par, mode_string);
        ape_test_cmp_string("ape_par_get_mode for mode \"qh\"", mode_string, "", status, eUnknownMode);
      } else {
        ape_test_failed("ape_par_test: assignment of string \"qh\" to par mode failed with status %d.\n", status);
      }
    } else {
      ape_test_failed("Setup for ape_par_get_mode(par, mode_string) failed with status %d.\n", status);
    }
    ape_par_destroy(par);
  }

  /* Test ape_par_get_field_array. */
#if 0
  { const char * par_text = "testpar,    i , a, 1 2   3   4 5  6, , , \"Test parameter\"";
    const char * expected[] = { "1", "2", "3", "4", "5", "6", 0 };
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    if (eOK == status) {
      char ** actual = 0;
      status = ape_par_get_field_array(par, eValue, &actual);
      ape_test_cmp_string_array("ape_par_get_field_array", actual, expected, status, eOK);
      ape_util_free_string_array(actual);
    } else {
      ape_test_failed("Setup for ape_par_get_field_array test failed with status %d.\n", status);
    }
    ape_par_destroy(par);
  }
#endif

  /* Test ape_par_get_min_string and ape_par_get_enum_string. */
  { const char * par_text = "testpar,    i , a, 1, 1| 2 |3|4, , \"Test parameter\"";
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    if (eOK == status) {
      char * actual = 0;
      const char * expected = "1| 2 |3|4";
      char ** actual_range = 0;
      const char * expected_range[] = { "1", "2", "3", "4", 0 };

      /* Get minimum field as a string, which should work, but return a non-eOK status because the
         minimum field contains an enumerated list. */
      status = ape_par_get_min_string(par, &actual);
      ape_test_cmp_string("ape_par_get_min_string", actual, expected, status, eRangeEnum);
      free(actual); actual = 0;

      /* Get minimum field as a set of enumerated values, which should work and return eOK because
         that is how enumerated range is supposed to be encoded. */
      status = ape_par_get_enum_string(par, &actual_range);
      ape_test_cmp_string_array("ape_par_get_enum_string", actual_range, expected_range, status, eOK);
      ape_util_free_string_array(actual_range); actual_range = 0;

      /* Test a non-enumerated range with a format error. */
      expected = "1 2";
      expected_range[0] = "1 2";
      expected_range[1] = 0;
      /* Note the space makes this illegal for integer parameter, but ape_par_get_enum_string does not see
         this error, because it treats the field as a string. */
      status = ape_par_set_field(par, eMin, "1 2");
      if (eOK == status) {
        status = ape_par_get_min_string(par, &actual);
        ape_test_cmp_string("ape_par_get_min_string", actual, expected, status, eOK);
        free(actual); actual = 0;

        status = ape_par_get_enum_string(par, &actual_range);
        ape_test_cmp_string_array("ape_par_get_enum_string", actual_range, expected_range, status, eRangeNoEnum);
        ape_util_free_string_array(actual_range); actual_range = 0;
      } else {
        ape_test_failed("Assignemt for ape_par_get_min_string test failed with status %d.\n", status);
      }
    } else {
      ape_test_failed("Setup for ape_par_get_min_string test failed with status %d.\n", status);
    }
    ape_par_destroy(par);
  }

  /* Test redirecting the prompt stream to a NULL stream, which should generate an error but have no effect. */
  { int status = ape_par_redirect_prompt_stream(0);
    ape_test_cmp_string("After ape_par_redirect_prompt_stream(0)", "", "", status, eNullPointer);
  }

  /* Test ape_par_prompt. */
  { const char * par_text = "testpar,    i , a, 1, 1| 2, , 'Enter \"test\" parameter'";
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    if (eOK == status) {
      /* Get input, expected value is 2. */
      status = ape_par_prompt(par);
      /* For cosmetic reasons, print a carriage return here. */
      ape_msg_debug("\n");

      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_prompt", value, "2", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("ape_par_prompt returned status %d, not eOK as expected.\n", status);
      }

      /* Get input, expected value is C/R (no text), which does not change the parameter. */
      status = ape_par_prompt(par);
      /* For cosmetic reasons, print a carriage return here. */
      ape_msg_debug("\n");

      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_prompt", value, "2", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("ape_par_prompt returned status %d, not eOK as expected.\n", status);
      }

      /* Get input, expected value is a single space, which blanks out the parameter. */
      status = ape_par_prompt(par);
      /* For cosmetic reasons, print a carriage return here. */
      ape_msg_debug("\n");

      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_prompt", value, "", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("ape_par_prompt returned status %d, not eOK as expected.\n", status);
      }

      /* Get input, expected value is a 1, which resets the parameter to its original value. */
      status = ape_par_prompt(par);
      /* For cosmetic reasons, print a carriage return here. */
      ape_msg_debug("\n");

      if (eOK == status) {
        char * value = 0;
        status = ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_prompt", value, "1", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("ape_par_prompt returned status %d, not eOK as expected.\n", status);
      }

    } else {
      ape_test_failed("Setup for ape_par_prompt test failed with status %d.\n", status);
    }
    ape_par_destroy(par);
  }

  /* Test ape_par_query. */
  { const char * par_text = "testpar,    i , a, 1, 1| 2, , \"Enter test parameter\"";
    ApePar * par = 0;
    int status = ape_par_create(&par, par_text);
    int orig_prompt_style = 0;
    if (eOK == status) {
      /* Get ready to test ape_par_query for hidden parameter. Set value to 1 and set prompt to be something
         clarifying that it's an error for a hidden parameter to be prompted for. */
      status = ape_par_set_field(par, eValue, "1");
      if (eOK == status) {
        status = ape_par_set_field(par, ePrompt, "Mode parameter is h; you should not see this. Enter 2 if you do");
      }
      if (eOK == status) {
        char * value = 0;
        /* Call ape_par_query, which should do nothing, since mode parameter is hidden. */
        status = ape_par_query(par, "h");
        ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_query called for with mode parameter == \"h\"", value, "1", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("Could not set up to test ape_par_query for a hidden parameter.\n");
      }

      /* Get ready to test ape_par_query for queried parameter. Set value to 1 and set prompt to be something
         clarifying the test. */
      status = ape_par_set_field(par, eValue, "1");
      if (eOK == status) {
        status = ape_par_set_field(par, ePrompt,
          "'Mode par is \"q\", you should see this. To test re-prompt, first enter 3, then 2'");
      }
      if (eOK == status) {
        char * value = 0;
        /* Call ape_par_query, which should prompt; if the user cooperates and enters first an invalid entry,
           reprompting should occur. */
        status = ape_par_query(par, "q");
        /* For cosmetic reasons, print a carriage return here. */
        ape_msg_debug("\n");

        ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_query called with mode parameter == \"q\"", value, "2", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("Could not set up to test ape_par_query for a queried parameter.\n");
      }

      /* Get ready to test ape_par_query for ranged parameter. Set value to 1 and set prompt to be something
         clarifying the test. */
      status = ape_par_set_field(par, eValue, "1");
      if (eOK == status) {
        status = ape_par_set_field(par, ePrompt,
          "'Test re-prompts: enter -8000000000, 8000000000, blah, -1, 3, 1., 2'");
      }
      if (eOK == status) {
        status = ape_par_set_field(par, eMin, "0");
      }
      if (eOK == status) {
        status = ape_par_set_field(par, eMax, "2");
      }
      if (eOK == status) {
        char * value = 0;
        /* Call ape_par_query, which should prompt; if the user cooperates and enters first an invalid entry,
           reprompting should occur. */
        status = ape_par_query(par, "q");
        /* For cosmetic reasons, print a carriage return here. */
        ape_msg_debug("\n");

        ape_par_get_field(par, eValue, &value);
        ape_test_cmp_string("ape_par_query called for re-prompt test", value, "2", status, eOK);
        free(value); value = 0;
      } else {
        ape_test_failed("Could not set up to test ape_par_query for re-prompting.\n");
      }

      /* Disable prompts and prompt again, expected value is again 2, with no error despite the fact that
         stdin is at EOF. */
      status = ape_par_get_default_prompt_style(&orig_prompt_style);
      if (eOK == status) {
        status = ape_par_set_default_prompt_style(eNoPrompt);
      }
      if (eOK == status) status = ape_par_set_field(par, eValue, "1");
      if (eOK == status) {
        status = ape_par_set_field(par, ePrompt, "'Prompt style is no prompt. you should not see this. Enter 2 if you do'");
      }
      if (eOK == status) {
        status = ape_par_prompt(par);
        /* For cosmetic reasons, print a carriage return here. */
        ape_msg_debug("\n");

        if (eOK == status) {
          char * value = 0;
          status = ape_par_get_field(par, eValue, &value);
          ape_test_cmp_string("After disabling prompts, ape_par_prompt", value, "1", status, eOK);
          free(value); value = 0;
        } else {
          ape_test_failed("After disabling prompts, ape_par_prompt returned status %d, not eOK as expected.\n", status);
        }
      } else {
        ape_test_failed("ape_par_set_default_prompt_style(NoPrompt) returned status %d, not eOK as expected.\n", status);
      }
      /* Enable prompts and prompt again, which should generate an error because there is no more input. */
      if (eOK == status) {
        status = ape_par_set_default_prompt_style(eDefaultPrompt);
      }
      if (eOK == status) status = ape_par_set_field(par, eValue, "1");
      if (eOK == status) {
        status = ape_par_set_field(par, ePrompt, "'Test EOF condition; hit Ctrl-D'");
      }
      if (eOK == status) {
        status = ape_par_prompt(par);
        /* For cosmetic reasons, print a carriage return here. */
        ape_msg_debug("\n");

        if (eInputFailure == status) {
          char * value = 0;
          status = ape_par_get_field(par, eValue, &value);
          ape_test_cmp_string("After reenabling prompts, ape_par_prompt", value, "1", status, eOK);
          free(value); value = 0;
        } else {
          ape_test_failed("After reenabling prompts, ape_par_prompt returned status %d, not %d as expected.\n",
            status, eInputFailure);
        }
      } else {
        ape_test_failed("ape_par_set_default_prompt_style(eNoPrompt) returned status %d, not eOK as expected.\n", status);
      }
      /* Restore previous status quo. */
      ape_par_set_default_prompt_style(orig_prompt_style);
    } else {
      ape_test_failed("Setup for ape_par_query test failed with status %d.\n", status);
    }
    ape_par_destroy(par);
  }

  /* Test parameter checking without checking value. */
  test_par_check("", eOK, 0);
  test_par_check(" ", eOK, 0);
  test_par_check(" #", eOK, 0);
  test_par_check("parname,  b,  a,   ,    ,    ,", eOK, 0);
  test_par_check("parname,  d, la,   ,    ,  10, Prompt", eOK, 0);
  test_par_check("parname,  i,  q,   ,   0,    , Prompt", eOK, 0);
  test_par_check("parname,  r, ql,   ,   0,  10, Prompt", eOK, 0);
  test_par_check("parname,  i,  h,  1,    ,    , Prompt", eOK, 0);
  test_par_check("parname,  i, hl,  1,    ,  10, Prompt", eOK, 0);
  test_par_check("parname,  i,  a,  1,   0,  10, Prompt", eOK, 0);
  test_par_check("parname,  i,  a,  1,   0,  10, \"Prompt\\\",\\\" with quoted comma\"", eOK, 0);
  test_par_check("       ,  b,  a,   ,    ,    ,", eInvalidName, 0);
  test_par_check("parname,  i,  a,  1,   0,  10", eFieldNotFound, 0);
  test_par_check("parname,  i,  a,  1,   0,  10, \"Prompt\",\" plus extra comma\"", eTooManyFields, 0);
  test_par_check("parname,  Z,  a,  1,   0,  10, Prompt", eUnknownType, 0);
  test_par_check("parname,  i,  Z,  1,   0,  10, Prompt", eUnknownMode, 0);
  test_par_check("parname,  i,  a, 1.,   0,  10, Prompt", eTypeMismatch, 1);
  test_par_check("parname,  i,  a, 1.,   0,  10, Prompt", eOK, 0); /* For now, ape_par_check does not check value. */
  test_par_check("parname,  i,  a,  1,  0.,  10, Prompt", eMinConversionError, 0);
  test_par_check("parname,  i,  a,  1,   0, 10., Prompt", eMaxConversionError, 0);
  test_par_check("parname,  i,  a,one,   0,  10, Prompt", eStringRemainder, 1);
  test_par_check("parname,  i,  a,one,   0,  10, Prompt", eOK, 0); /* For now, ape_par_check does not check value. */
  test_par_check("parname,  i,  a,  1,zero,  10, Prompt", eMinConversionError, 0);
  test_par_check("parname,  i,  a,  1,   0, ten, Prompt", eMaxConversionError, 0);
  test_par_check("parname,  i,  a,  1, 0|1,  10, Prompt", eInvalidRange, 0);
  test_par_check("parname,  i,  a,  1,  10,   0, Prompt", eInvalidRange, 0);
  test_par_check("parname,  i,  a,  1,   0, 0|1, Prompt", eMaxConversionError, 0);
  test_par_check("parname,  i,  a, -1,   0,  10, Prompt", eValueBelowMin, 1);
  test_par_check("parname,  i,  a, -1,   0,  10, Prompt", eOK, 0); /* For now, ape_par_check does not check value. */
  test_par_check("parname,  i,  a, 11,   0,  10, Prompt", eValueAboveMax, 1);
  test_par_check("parname,  i,  a, 11,   0,  10, Prompt", eOK, 0); /* For now, ape_par_check does not check value. */
  test_par_check("parname,  i,  a,  2, 0|1,    , Prompt", eInvalidChoice, 1);
  test_par_check("parname,  i,  a,  2, 0|1,    , Prompt", eOK, 0); /* For now, ape_par_check does not check value. */
  test_par_check("parname,  i,  h,  -1,    ,    , Prompt", eOK, 0); /* Make sure blank minimum not treated as 0. */

  remove("non-existent-file"); /* Just in case. */
  test_par_check("parname,  fw,  h,  non-existent-file,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  fr,  h,  non-existent-file,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  f,  h,  non-existent-file,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  fw,  h,  ape_test.par,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  fr,  h,  ape_test.par,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  f,  h,  ape_test.par,    ,    , Prompt", eOK, 1);
  remove("non-existent-file"); /* Just in case. */
  test_par_check("parname,  fwe,  h,  non-existent-file,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  fre,  h,  non-existent-file,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  fe,  h,  non-existent-file,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  fwe,  h,  ape_test.par,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  fre,  h,  ape_test.par,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  fe,  h,  ape_test.par,    ,    , Prompt", eOK, 1);
  remove("non-existent-file"); /* Just in case. */
  test_par_check("parname,  fwN,  h,  non-existent-file,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  frN,  h,  non-existent-file,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  fN,  h,  non-existent-file,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  fwN,  h,  ape_test.par,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  frN,  h,  ape_test.par,    ,    , Prompt", eFileNotAccessible, 1);
  test_par_check("parname,  fn,  h,  ape_test.par,    ,    , Prompt", eFileNotAccessible, 1);

  /* Test parameter checking including checking value. */
  test_par_check("", eOK, 1);
  test_par_check(" ", eOK, 1);
  test_par_check(" #", eOK, 1);
  test_par_check("parname,  b,  a,   ,    ,    ,", eOK, 1);
  test_par_check("parname,  d, la,   ,    ,  10, Prompt", eOK, 1);
  test_par_check("parname,  i,  q,   ,   0,    , Prompt", eOK, 1);
  test_par_check("parname,  r, ql,   ,   0,  10, Prompt", eOK, 1);
  test_par_check("parname,  i,  h,  1,    ,    , Prompt", eOK, 1);
  test_par_check("parname,  i, hl,  1,    ,  10, Prompt", eOK, 1);
  test_par_check("parname,  i,  a,  1,   0,  10, Prompt", eOK, 1);
  test_par_check("parname,  i,  a,  1,   0,  10, \"Prompt\\\",\\\" with quoted comma\"", eOK, 1);
  test_par_check("       ,  b,  a,   ,    ,    ,", eInvalidName, 1);
  test_par_check("parname,  i,  a,  1,   0,  10", eFieldNotFound, 1);
  test_par_check("parname,  i,  a,  1,   0,  10, \"Prompt\",\" plus extra comma\"", eTooManyFields, 1);
  test_par_check("parname,  Z,  a,  1,   0,  10, Prompt", eUnknownType, 1);
  test_par_check("parname,  i,  Z,  1,   0,  10, Prompt", eUnknownMode, 1);
  test_par_check("parname,  i,  a, 1.,   0,  10, Prompt", eTypeMismatch, 1);
  test_par_check("parname,  i,  a,  1,  0.,  10, Prompt", eMinConversionError, 1);
  test_par_check("parname,  i,  a,  1,   0, 10., Prompt", eMaxConversionError, 1);
  test_par_check("parname,  i,  a,one,   0,  10, Prompt", eStringRemainder, 1);
  test_par_check("parname,  i,  a,  1,zero,  10, Prompt", eMinConversionError, 1);
  test_par_check("parname,  i,  a,  1,   0, ten, Prompt", eMaxConversionError, 1);
  test_par_check("parname,  i,  a,  1, 0|1,  10, Prompt", eInvalidRange, 1);
  test_par_check("parname,  i,  a,  1,  10,   0, Prompt", eInvalidRange, 1);
  test_par_check("parname,  i,  a,  1,   0, 0|1, Prompt", eMaxConversionError, 1);
  test_par_check("parname,  i,  a, -1,   0,  10, Prompt", eValueBelowMin, 1);
  test_par_check("parname,  i,  a, 11,   0,  10, Prompt", eValueAboveMax, 1);
  test_par_check("parname,  i,  a,  2, 0|1,    , Prompt", eInvalidChoice, 1);
  test_par_check("parname,  i,  h,  -1,    ,    , Prompt", eOK, 1); /* Make sure blank minimum not treated as 0. */
  test_par_check("parname,  g,  h,  0-1,    ,    , Prompt", eOK, 1); /* Test xselect's "g" type, which should be tolerated. */
  /* Test behavior with INDEF, with and without a range specified. */
  test_par_check("parname,  r,  h,  INDEF,   ,   , Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INDEF, 0.,   , Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INDEF,   , 5., Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INDEF, 0., 5., Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INDEF, 0.|1., , Prompt", eOK, 1);
  /* Test behavior with INDEF in ranges, which is not allowed. */
  test_par_check("parname,  r,  h,  INDEF, INDEF, 5., Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  1., INDEF, 5., Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  6., INDEF, 5., Prompt", eInvalidRange, 1);

  test_par_check("parname,  r,  h,  INDEF, 0., INDEF, Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  1., 0., INDEF, Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h, -1., 0., INDEF, Prompt", eInvalidRange, 1);

  test_par_check("parname,  r,  h,  INDEF, INDEF, INDEF, Prompt", eInvalidRange, 1);

  test_par_check("parname,  r,  h,  INDEF, 1.|INDEF, , Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  1., 1.|INDEF, , Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  5., 1.|INDEF, , Prompt", eInvalidRange, 1);

  /* Test behavior with INF, with and without a range specified. */
  test_par_check("parname,  r,  h,  INF,   ,   , Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INF, 0.,   , Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INF,   , 5., Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INF, 0., 5., Prompt", eOK, 1);
  test_par_check("parname,  r,  h,  INF, 0.|1., , Prompt", eOK, 1);
  /* Test behavior with INF in ranges, which is not allowed. */
  test_par_check("parname,  r,  h,  INF, INF, 5., Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  1., INF, 5., Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  6., INF, 5., Prompt", eInvalidRange, 1);

  test_par_check("parname,  r,  h,  INF, 0., INF, Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  1., 0., INF, Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h, -1., 0., INF, Prompt", eInvalidRange, 1);

  test_par_check("parname,  r,  h,  INF, INF, INF, Prompt", eInvalidRange, 1);

  test_par_check("parname,  r,  h,  INF, 1.|INF, , Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  1., 1.|INF, , Prompt", eInvalidRange, 1);
  test_par_check("parname,  r,  h,  5., 1.|INF, , Prompt", eInvalidRange, 1);

  /* Tests of ape_par_get_eff_mode. */
  { const char * par_string[] = {
      "a, i,  a, , , , auto",
      "al, i, al, , , , auto learned",
      " h, i,  h, , , , hidden",
      "hl, i, hl, , , , hidden learned",
      " q, i,  q, , , , queried",
      "ql, i, ql, , , , queried learned"
    };
    const char * mode_par[] = { "a", "al", "h", "hl", "q", "ql" };
    const char * eff_mode_expected[][sizeof(mode_par) / sizeof(mode_par[0])] = {
      { "a", "al", "h", "hl", "q", "ql" },
      { "al", "al", "hl", "hl", "ql", "ql" },
      { "h", "h", "h", "h", "h", "h" },
      { "hl", "hl", "hl", "hl", "hl", "hl" },
      { "q", "q", "q", "q", "q", "q" },
      { "ql", "ql", "ql", "ql", "ql", "ql" }
    };
    const char * eff_mode_expected_after_set[][sizeof(mode_par) / sizeof(mode_par[0])] = {
      { "h", "hl", "h", "hl", "h", "hl" },
      { "hl", "hl", "hl", "hl", "hl", "hl" },
      { "h", "h", "h", "h", "h", "h" },
      { "hl", "hl", "hl", "hl", "hl", "hl" },
      { "h", "h", "h", "h", "h", "h" },
      { "hl", "hl", "hl", "hl", "hl", "hl" }
    };
    size_t ii = 0;
    size_t jj = 0;

    for (ii = 0; ii != sizeof(mode_par) / sizeof(mode_par[0]); ++ii) {
      ApePar * par = 0;
      int status = ape_par_create(&par, par_string[ii]);
      if (eOK == status) {
        for (jj = 0; jj != sizeof(mode_par) / sizeof(mode_par[0]); ++jj) {
          char msg[128] = "";
          char eff_mode[APE_PAR_MODE_CODE_LEN] = "";
          status = ape_par_get_eff_mode(par, mode_par[jj], eff_mode);
          sprintf(msg, "ape_par_get_eff_mode(\"%s\") for mode \"%s\"", par_string[ii], mode_par[jj]);
          ape_test_cmp_string(msg, eff_mode, eff_mode_expected[ii][jj], status, eOK);
        }
        status = ape_par_flag_cmd_line(par);
        if (eOK == status) {
          for (jj = 0; jj != sizeof(mode_par) / sizeof(mode_par[0]); ++jj) {
            char msg[128] = "";
            char eff_mode[APE_PAR_MODE_CODE_LEN] = "";
            status = ape_par_get_eff_mode(par, mode_par[jj], eff_mode);
            sprintf(msg, "after ape_par_flag_cmd_line, ape_par_get_eff_mode(\"%s\") for mode \"%s\"", par_string[ii], mode_par[jj]);
            ape_test_cmp_string(msg, eff_mode, eff_mode_expected_after_set[ii][jj], status, eOK);
          }
        } else {
          ape_test_failed("Unable to set parameter value to test ape_par_get_eff_mode (status was &d).\n", status);
        }
      } else {
        ape_test_failed("Unable to set up to test ape_par_get_eff_mode (status was &d).\n", status);
      }
      ape_par_destroy(par); par = 0;
    }
  }
  /* Test collapse_escape. */
  { const char * input = "\\\\ \\ \" \\\" ' \\' \\, \\t \\q \\z \\n";
    const char * expected_output = "\\  \" \" ' ' , \t q \\z \n";

    char * output = 0;
    int status = collapse_escape(input, "'\", tnq\\", &output);
    ape_test_cmp_string("collapse_escape", output, expected_output, status, eOK);
    free(output); output = 0;
  }
  /* Test expand_escape. */
  { const char * input = "\"this\" and 'that', \"&& something else.";
    const char * expected_output = "\"this\" and 'that', \"&& something else.";
    char * output = 0;
    int status = expand_escape(input, input + strlen(input), "", &output);
    ape_test_cmp_string("expand_escape(input, input + strlen(input), \"\" &output)", output, expected_output, status, eOK);
    free(output); output = 0;

    input = "\"this\" and 'that', \"&& something else.";
    expected_output = "\\\"this\\\" and \\'that\\'\\, \\\"&& something else.";
    output = 0;
    status = expand_escape(input, input + strlen(input), "\"',", &output);
    ape_test_cmp_string("expand_escape(input, input + strlen(input), \"\\\"',\", &output)", output, expected_output, status, eOK);
    free(output); output = 0;

    input = "\"this\" and 'that', \"&& something else.";
    expected_output = "\\\"this\\\" and 'that', \\\"&& something else.";
    output = 0;
    status = expand_escape(input, input + strlen(input), "\"", &output);
    ape_test_cmp_string("expand_escape(input, input + strlen(input), \"\\\"\"", output, expected_output, status, eOK);
    free(output); output = 0;
  }
  /* Test quote_field. */
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "sq, s, q, \"quoted\"\\, but \\\"that\\'s not all, , , prompt";
    char * line = 0;
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      const_line = "sq, s, q, \"\\\"quoted\\\", but \\\"that's not all\", , , prompt";
      status = quote_field(ape_par, eValue);
      if (eOK == status) line = ape_par_get_line(ape_par);
      ape_test_cmp_string("quote_field", line, const_line, status, eOK);
    } else {
      ape_test_failed("Could not set up to test quote_field, status is %d.\n", status);
    }
    free(line); line = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "sq, s, q, \"quoted\\, but \\\"that's not all\", , , prompt";
    char * line = 0;
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = quote_field(ape_par, eValue);
      if (eOK == status) line = ape_par_get_line(ape_par);
      ape_test_cmp_string("quote_field", line, const_line, status, eOK);
    } else {
      ape_test_failed("Could not set up to test quote_field, status is %d.\n", status);
    }
    free(line); line = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }
  /* Test ape_par_set_value_string. */
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "sq, s, q, value\\, \\\" with escaped quote, , , prompt";
    char * line = 0;
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      const_line = "sq, s, q, \"value, \\\" with escaped quote\", , , prompt";
      status = ape_par_set_value_string(ape_par, "value, \" with escaped quote");
      if (eOK == status) line = ape_par_get_line(ape_par);
      ape_test_cmp_string("ape_par_set_value_string", line, const_line, status, eOK);
    } else {
      ape_test_failed("Could not set up to test ape_par_set_value_string, status is %d.\n", status);
    }
    free(line); line = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "iq, i, q, , , , prompt";
    char * line = 0;
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      const_line = "iq, i, q, 10, , , prompt";
      status = ape_par_set_value_string(ape_par, "10");
      if (eOK == status) line = ape_par_get_line(ape_par);
      ape_test_cmp_string("ape_par_set_value_string", line, const_line, status, eOK);
    } else {
      ape_test_failed("Could not set up to test ape_par_set_value_string, status is %d.\n", status);
    }
    free(line); line = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }

  /* Test ape_par_get_name. */
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = " par name , s, a, , , , prompt";
    char * name = 0;
    const char * expected_name = "par name";
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_get_name(ape_par, &name);
      ape_test_cmp_string("ape_par_get_name", name, expected_name, status, eOK);
    } else {
      ape_test_failed("Could not set up to test ape_par_get_name for a valid parameter, status is %d.\n", status);
    }
    free(name); name = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = " par name only\t";
    char * name = 0;
    const char * expected_name = "par name only";
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_get_name(ape_par, &name);
      ape_test_cmp_string("ape_par_get_name", name, expected_name, status, eOK);
    } else {
      ape_test_failed("Could not set up to test ape_par_get_name for an invalid but named parameter, status is %d.\n",
        status);
    }
    free(name); name = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = " # comment only";
    char * name = 0;
    const char * expected_name = 0;
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_get_name(ape_par, &name);
      ape_test_cmp_string("ape_par_get_name", name, expected_name, status, eUnnamedPar);
    } else {
      ape_test_failed("Could not set up to test ape_par_get_name for a valid but unnamed parameter, status is %d.\n",
        status);
    }
    free(name); name = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "   \t   ";
    char * name = 0;
    const char * expected_name = 0;
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_get_name(ape_par, &name);
      ape_test_cmp_string("ape_par_get_name", name, expected_name, status, eUnnamedPar);
    } else {
      ape_test_failed("Could not set up to test ape_par_get_name for a blank parameter, status is %d.\n",
        status);
    }
    free(name); name = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }

  /* Test ape_par_get_comment with no comment. */
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "iq, i, q, , , , prompt  ";
    char * comment = 0;
    const char * expected_comment = "";
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_get_comment(ape_par, &comment);
      ape_test_cmp_string("ape_par_get_comment", comment, expected_comment, status, eOK);
    } else {
      ape_test_failed("Could not set up to test ape_par_get_comment for no comment, status is %d.\n", status);
    }
    free(comment); comment = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }

  /* Test ape_par_get_comment with comment. */
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "iq, i, q, , , , prompt  #   comment goes here";
    char * comment = 0;
    const char * expected_comment = "#   comment goes here";
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_get_comment(ape_par, &comment);
      ape_test_cmp_string("ape_par_get_comment", comment, expected_comment, status, eOK);
    } else {
      ape_test_failed("Could not set up to test ape_par_get_comment, status is %d.\n", status);
    }
    free(comment); comment = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }

  /* Test ape_par_set_comment with comment. */
  { ApePar * ape_par = 0;
    int status = eOK;
    const char * const_line = "iq, i, q, , , , prompt  #   comment goes here";
    char * comment = 0;
    const char * expected_comment = "#   new comment went here";
    status = ape_par_create(&ape_par, const_line);
    if (eOK == status) {
      status = ape_par_set_comment(ape_par, expected_comment);
      if (eOK == status) {
        status = ape_par_get_comment(ape_par, &comment);
        ape_test_cmp_string("ape_par_set_comment", comment, expected_comment, status, eOK);
      } else {
        ape_test_failed("ape_par_set_comment failed unexpectedly with status %d.\n", status);
      }
    } else {
      ape_test_failed("Could not set up to test ape_par_set_comment, status is %d.\n", status);
    }
    free(comment); comment = 0;
    ape_par_destroy(ape_par); ape_par = 0;
  }

  /* Test ape_par_register_get_text. */
  {
    char * test_value = 0;
    int status = eOK;
    status = ape_par_register_get_text(&test_register_get_text);
    status = custom_get_text("Test Prompt","Test Name",&test_value);
    if (0 != strcmp(test_value,"Hello World"))
      ape_test_failed("ape_par_register_get_text failed to return \'Hello World\'\n");
    if (eInputFailure != status)
      ape_test_failed("ape_par_register_get_text failed unexpectedly with status %d\n",status);
  }
}

static void test_one_par(const char * line, const char ** correct_field) {
  ApePar * ape_par = 0;
  int status = ape_par_create(&ape_par, line);
  if (0 != status) {
    ape_test_failed("Calling ape_par_create(&ape_par, \"%s\") returned non-0 status.\n", line);
  } else {
    size_t field_index = 0;
    size_t correct_field_index = 0;
    char * field = 0;
    if (0 != compare_par(ape_par, line)) {
      char * par_line = ape_par_get_line(ape_par);
      ape_test_failed("Just after ape_par_create, parameter was %s, not %s as expected.\n", par_line, line);
      free(par_line);
    }
    for (field_index = 0; field_index != eEndOfField; ++field_index) {
      status = get_field(ape_par, field_index, &field);
      if (field != 0 && correct_field[correct_field_index] == 0) {
        ape_test_failed("In line |%s|, field %d was parsed to be |%s| but should have been 0.\n", line, field_index, field);
      } else if (field == 0 && correct_field[correct_field_index] != 0) {
        ape_test_failed("In line |%s|, field %d was parsed to be 0, but should have been |%s|.\n", line, field_index,
          correct_field[correct_field_index]);
      } else if (field != 0 && 0 != strcmp(field, correct_field[correct_field_index])) {
        ape_test_failed("In line |%s|, field %d was parsed to be |%s|, but should have been |%s|.\n", line, field_index, field,
          correct_field[correct_field_index]);
      }
      if (0 != correct_field[correct_field_index]) ++correct_field_index;
      free(field);
    }
    /* Make sure getting line from parameter returns the same as the original line. */
    field = ape_par_get_line(ape_par);
    if (0 != strcmp(line, field)) {
      ape_test_failed("When called for an ApePar created from line |%s|, ape_par_get_line unexpectedly gave |%s|.\n",
        line, field);
    }
    free(field);
  }

  ape_par_destroy(ape_par);
}

static void test_par_check(const char * par_text, int expected_status, char check_value) {
  int status = eOK;
  ApePar * par = 0;
  status = ape_par_create(&par, par_text);
  if (eOK == status) {
    status = ape_par_check(par, check_value);
    if (status != expected_status)
      ape_test_failed("ape_par_check(\"%s\", %d) returned status %d, not %d as expected.\n\n", par_text, check_value,
        status, expected_status);
    else
      ape_msg_debug("ape_par_check(\"%s\", %d) behaved as expected.\n", par_text, check_value);
  } else {
    ape_test_failed("Could not set up for ape_par_check(\"%s\", %d).\n", par_text, check_value);
  }
  ape_par_destroy(par); par = 0;
}

static int compare_par(const ApePar * ape_par, const char * correct_line) {
  int difference_found = 0;
  if (0 == ape_par) {
    ape_test_failed("Test logic error: compare_par was passed null parameter pointer.\n");
  } else {
    char * line = ape_par_get_line(ape_par);
    if (0 == line && 0 == correct_line) {
      /* OK: both NULL. */
    } else if (0 == line && 0 != correct_line) {
      ape_msg_debug("compare_par: %s.\n", "compare_par: parameter's line was NULL, but it should be");
      ape_msg_debug("compare_par: %s.\n", correct_line);
      difference_found = 1;
    } else if (0 != line && 0 == correct_line) {
      ape_msg_debug("compare_par: %s.\n", "parameter's line should be NULL, but it was");
      ape_msg_debug("compare_par: %s.\n", line);
      difference_found = 1;
    } else if (0 != strcmp(line, correct_line)) {
      ape_msg_debug("compare_par: %s.\n", "parameter's line was");
      ape_msg_debug("compare_par: %s.\n", line);
      ape_msg_debug("compare_par: %s.\n", "but it should have been");
      ape_msg_debug("compare_par: %s.\n", correct_line);
      difference_found = 1;
    }
    free(line);
  }
  return difference_found;
}

static int collapse_escape(const char * input, const char * escape, char ** output) {
  int status = eOK;

  /* Check arguments. */
  if (0 != output) *output = 0;
  else status = eNullPointer;

  if (eOK == status && (0 == input || 0 == escape)) status = eNullPointer;

  if (eOK == status) {
    /* Allocate room for output. In case there are no escapes, make it as long as the input; this is an upper limit. */
    *output = (char *) calloc(strlen(input) + 1, sizeof(char));
    if (0 == *output) status = eDynAllocFailed;
  }

  if (eOK == status) {
    const char * in_p = input;
    const char * escape_p = 0;
    char * out_p = *output;
    for (; '\0' != *in_p; ++in_p, ++out_p) {
      if ('\\' == *in_p) {
        /* Check whether the next character is one which needs to be converted. */
        for (escape_p = escape; '\0' != *escape_p && *(in_p + 1) != *escape_p; ++escape_p) {}

        /* If this character is not one of the known escape codes, just copy it. */
        if ('\0' == *escape_p) {
          *out_p = *in_p;
        } else {
          ++in_p;
          /* TODO Add all escape codes. */
          if ('b' == *in_p) *out_p = '\b';
          else if ('n' == *in_p) *out_p = '\n';
          else if ('t' == *in_p) *out_p = '\t';
          else if ('v' == *in_p) *out_p = '\v';
          else *out_p = *in_p;
        }
      } else {
        *out_p = *in_p;
      }
    }
  }

  return status;
}

static int expand_escape(const char * begin, const char * end, const char * escape_me, char ** output) {
  int status = eOK;
  size_t num_escape = 0;
  const char * in_p = begin;

  /* Check arguments. */
  if (0 != output) *output = 0;
  else status = eNullPointer;

  if (eOK == status && (0 == begin || 0 == end || 0 == escape_me)) status = eNullPointer;
  if (eOK == status && begin > end) status = eInvalidArgument;

  if (eOK == status) {
    for (; in_p != end; ++in_p) {
      const char * escape_p = escape_me;
      for (; '\0' != *escape_p && *in_p != *escape_p; ++escape_p) {}
      if ('\0' != *escape_p) ++num_escape;
    }
    *output = (char *) calloc(end - begin + num_escape + 1, sizeof(char));
    if (0 == *output) status = eDynAllocFailed;
  }

  if (eOK == status) {
    char * out_p = *output;
    for (in_p = begin; in_p != end; ++in_p, ++out_p) {
      const char * escape_p = escape_me;
      for (; '\0' != *escape_p && *in_p != *escape_p; ++escape_p) {}
      if ('\0' != *escape_p) {
        *out_p = '\\';
        ++out_p;
      }
      *out_p = *in_p;
    }
  }

  return status;
}

/* Make sure the given field is quoted; add double quotes if none are present. */
static int quote_field(const ApePar * ape_par, ParFieldId field_id) {
  int status = eOK;
  char ** token_list = 0;
  int quote = eNoQuote;

  if (0 == ape_par || 0 == ape_par->field) {
    status = eNullPointer;
  }

  if (eOK == status && (field_id < eName || ePrompt < field_id)) {
    status = eInvalidArgument;
  }

  if (eOK == status) {
    /* Extract the significant part of the field. */
    status = find_field(ape_par, field_id, &token_list, &quote);
  }

  /* If field is not already quoted, it's necessary to re-parse this token with quotes added. */
  if (eOK == status && eNoQuote == quote) {
    char * leading_quote = 0;
    char * trailing_quote = 0;
    char * tmp_field1 = 0;
    char * tmp_field2 = 0;

    /* Don't quote blank values. */
    if ('\0' == *token_list[eTokenValue]) return status;

    /* Handle escaped characters. (Remove one level of backslashes from special escape sequences.) */
    status = collapse_escape(token_list[eTokenValue], "'\",\\", &tmp_field1);

    if (eOK == status) {
      leading_quote = (char *) calloc(2, sizeof(char));
      if (0 != leading_quote) *leading_quote = '"';
      else status = eDynAllocFailed;
    }

    /* Re-escape just double quotes and \\. */
    if (eOK == status) expand_escape(tmp_field1, tmp_field1 + strlen(tmp_field1), "\"\\", &tmp_field2);

    if (eOK == status) {
      trailing_quote = (char *) calloc(2, sizeof(char));
      if (0 != trailing_quote) *trailing_quote = '"';
      else status = eDynAllocFailed;
    }

    if (eOK == status) {
      free(token_list[eTokenTrailingQuote]); token_list[eTokenTrailingQuote] = trailing_quote;
      free(token_list[eTokenValue]); token_list[eTokenValue] = tmp_field2;
      free(token_list[eTokenLeadingQuote]); token_list[eTokenLeadingQuote] = leading_quote;
    } else {
      free(trailing_quote); trailing_quote = 0;
      free(tmp_field2); tmp_field2 = 0;
      free(leading_quote); leading_quote = 0;
    }

    /* Clean up. */
    free(tmp_field1); tmp_field1 = 0;
  }

  return status;
}

/* To give a picture: each line in the par file is assumed to have the form
   field0, field1, field2, ..., where each field is of the form
   prefix-token value-token suffix-token
   where prefix-token = leading white space + opening quote, if any, and
   suffix-token = closing quote + trailing white space, if any.
*/
static int parse_par(const char * line, ApePar * ape_par) {
  int status = eOK;
  const char * begin = line;
  const char * end = 0;
  ApeList * field = 0;
  char ** token_list = 0;

  /* Clear out the current list held by the parameter. */
  ape_list_destroy(ape_par->field);

  /* Create a new list of fields for the parameter. */
  ape_par->field = ape_list_create();
  if (0 == ape_par->field) {
    ape_msg_debug("parse_par: %s.\n", "allocation failed for list of fields for line");
    ape_msg_debug("parse_par: %s.\n", line);
    status = eDynAllocFailed;
    return status;
  }

  /* For convenience, make a local shortcut for the list, and an iterator. */
  field = ape_par->field;

  do {
    /* Find next unquoted comma, which ends the field. */
    end = find_unquoted(begin, ",");

    /* Create a new list of tokens for this field, including terminating 0, and put it in the list of fields. */
    token_list = (char **) calloc(eNumTokens, sizeof(char *));

    if (0 == token_list) {
      ape_msg_debug("parse_par: %s.\n", "allocation failed for list of tokens in parameter in line at");
      ape_msg_debug("parse_par: %s.\n", begin);
      status = eDynAllocFailed;
      break;
    }
    ape_list_append(field, token_list);

    /* Parse the line until the end of a field, that is, the end of this list of tokens. */
    status = parse_field(begin, end, token_list);

    /* Continue parsing after this comma. */
    begin = end;
    if ('\0' != *begin) ++begin;
  } while ('\0' != *begin || ',' == *end);

  return status;
}

static int parse_field(const char * begin, const char * end, char ** token_list) {
  int status = eOK;
  const char * v_begin = begin;
  const char * v_end = end;
  const char * prefix_end = 0;
  const char * suffix_begin = 0;
  char * token = 0;

  /* Make sure range is valid. */
  if (begin > end) status = eInvalidArgument;

  if (eOK == status) {
    /* Locate beginning of significant characters in the field. */
    while (v_begin < v_end && 0 != isspace(*v_begin)) ++v_begin;

    /* Beginning of field signifies end of prefix (leading space). */
    prefix_end = v_begin;

    /* Locate end of significant characters in the field. */
    while (v_end > v_begin && 0 != isspace(*(v_end - 1))) --v_end;

    /* End of field signifies beginning of suffix (trailing space). */
    suffix_begin = v_end;

    /* See if field begins and ends with distinct matching quotes. If so, they will be handled separately. */
    if (1 < (v_end - v_begin) && *v_begin == *(v_end - 1) && ('\'' == *v_begin || '"' == *v_begin)) {
      ++v_begin; --v_end;
    }
  }

  /* Create string to hold the prefix and add it to the token list. */
  if (eOK == status) status = ape_util_copy_range(begin, prefix_end, &token);
  if (eOK == status) token_list[eTokenPrefix] = token;

  /* Create string to hold leading quote, if any. */
  if (eOK == status) status = ape_util_copy_range(prefix_end, v_begin, &token);
  if (eOK == status) token_list[eTokenLeadingQuote] = token;

  /* Create string to hold the value and add it to the token list. */
  if (eOK == status) status = ape_util_copy_range(v_begin, v_end, &token);
  if (eOK == status) token_list[eTokenValue] = token;

  /* Create string to hold trailing quote, if any. */
  if (eOK == status) status = ape_util_copy_range(v_end, suffix_begin, &token);
  if (eOK == status) token_list[eTokenTrailingQuote] = token;

  /* Create string to hold the suffix and add it to the token list. */
  if (eOK == status) status = ape_util_copy_range(suffix_begin, end, &token);
  if (eOK == status) token_list[eTokenSuffix] = token;

  return status;
}

/* Sets iterators to the beginning and ending of significant content of the given field. */
static int find_field(const ApePar * ape_par, ParFieldId field_id, char *** token_list, int * quote) {
  int status = eOK;
  ApeList * field = ape_par->field;
  unsigned int field_index = 0;
  ApeListIterator field_itor = ape_list_begin(field);
  ApeListIterator end = ape_list_end(field);

  /* Note: This is not an API call; doesn't check arguments! */
  *token_list = 0;

  /* Iterate over the fields until the correct field is found, or the end of the field list is reached. */
  for (; field_itor != end && field_index != field_id; field_itor = ape_list_next(field_itor), ++field_index) {}

  /* If the field was found, set the output iterators accordingly. */
  if (field_itor != end) {
    char ** token = (char **) ape_list_get(field_itor);

    if (0 != token) {
      /* Flag quotes. */
      const char * cp = token[eTokenLeadingQuote];
      if ('\'' == *cp) *quote = eSingleQuote;
      else if ('"' == *cp) *quote = eDoubleQuote;
      else *quote = eNoQuote;

      /* Extract value string. */
      *token_list = token;
    }

  } else {
    status = eFieldNotFound;
  }

  return status;
}

static int get_field(const ApePar * ape_par, ParFieldId field_id, char ** field_string) {
  int status = eOK;
  char ** token_list = 0;
  const char * token = 0;
  char * tmp_field = 0;
  int quote = eNoQuote;

  /* Check arguments. */
  if (0 != field_string) {
    *field_string = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == ape_par) {
    status = eNullPointer;
  }

  if (eOK == status && (field_id < eName || ePrompt < field_id)) {
    status = eInvalidArgument;
  }

  if (eOK == status) {
    /* Extract the significant part of the field. */
    status = find_field(ape_par, field_id, &token_list, &quote);
    if (eOK == status) token = token_list[eTokenValue];
  }

  if (eOK == status && eSingleQuote != quote) {
    /* Expand environment variables. */
    status = ape_util_expand_env_var(token, &tmp_field);
    if (eOK == status) {
      token = tmp_field;
    }
  }

  if (eOK == status) {
    /* Handle escaped characters. (Remove one level of backslashes.) */
    status = collapse_escape(token, "'\",#$", field_string);
  }

  /* Clean up. */
  free(tmp_field); tmp_field = 0;

  return status;
}

static int is_blank(const char * s) {
  int blank = 1;
  if (0 != s) {
    for (; '\0' != *s && 0 != isspace(*s); ++s) {}
    /* If loop terminated by reaching the end without finding non-space, string is blank. */
    blank = ('\0' == *s) ? 1 : 0;
  }
  return blank;
}

#if 0
/* TODO: reimplement this to support vectors. */
static int get_field_array(const ApePar * ape_par, ParFieldId field_id, char *** field_array) {
  int status = eOK;
  ApeListIterator begin = 0;
  ApeListIterator end = 0;
  int quote = eNoQuote;

  /* Extract the significant part of the field. */
  status = find_field(ape_par, field_id, &begin, &end, &quote);
  if (eOK == status) {
    /* The field was found, so concatenate the non-white space tokens to form the returned value string. */
    ApeListIterator itor = 0;
    size_t num_tokens = 1; /* At least one token needed, for terminating null. */
    size_t idx = 0;

    /* The range spanned by the iterators contain tokens separated by white space. Make one pass to
       determine the size of the output. */
    for (itor = begin; itor != end; itor = ape_list_next(itor)) {
      /* In array context, white space is not significant. */
      const char * token = (const char *) ape_list_get(itor);
      if (0 == is_blank(token)) ++num_tokens;
    }

    /* Allocate output string. */
    *field_array = (char **) calloc(num_tokens, sizeof(char *));

    if (0 != *field_array) {
      /* Make second pass to copy the tokens to the output array. */
      for (itor = begin; eOK == status && itor != end; itor = ape_list_next(itor)) {
        const char * token = (const char *) ape_list_get(itor);
        if (0 == is_blank(token)) {
          status = ape_util_copy_string(token, (*field_array) + idx);
          ++idx;
        }
      }
    } else {
      status = eDynAllocFailed;
    }

    /* Clean up in the event of an error. */
    if (eOK != status) {
      ape_util_free_string_array(*field_array);
      *field_array = 0;
    }
  }
  return status;
}
#endif

static int set_field(const ApePar * ape_par, ParFieldId field_id, const char * field_text) {
  int status = eOK;
  ApeList * field = ape_par->field;
  unsigned int field_index = 0;
  ApeListIterator field_itor = 0;
  char ** token = 0;
  char ** new_token = 0;
  char * prefix = 0;
  char * lead_quote = 0;
  char * current_value = 0;
  char * trail_quote = 0;
  char * suffix = 0;
  const char * begin = field_text;
  const char * end = field_text + strlen(field_text);
  const char * delim = "'\",#\\"; /* If no unescaped quotes surround an expression, all these must be escaped. */
  char * value = 0;
  ApeListIterator field_end = ape_list_end(field);

  /* Iterate over the fields until the correct field is found, or the end of the field list is reached. */
  for (field_itor = ape_list_begin(field); field_itor != field_end && field_index != field_id;
    field_itor = ape_list_next(field_itor), ++field_index) {}

  if (field_itor == field_end) {
    status = eFieldNotFound;
  }

  if (eOK == status) {
    /* Get or create content of this field. */
    token = (char **) ape_list_get(field_itor);
    if (0 == token) {
      token = (char **) calloc(eNumTokens, sizeof(char *));
      if (0 == token) status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* Get 0th token, the leading space. */
    prefix = token[eTokenPrefix];
    if (0 == prefix) {
      /* No prefix, so create blank one. */
      prefix = (char *) calloc(1, sizeof(char));
      if (0 == prefix) status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* Get next token, optional leading quote. */
    lead_quote = token[eTokenLeadingQuote];
    if (0 == lead_quote) {
      /* No leading quote, so create one. */
      lead_quote = (char *) calloc(2, sizeof(char));
      if (0 == lead_quote) status = eDynAllocFailed;
      *lead_quote = '"';
    }
  }

  if (eOK == status) {
    /* Go to next token in field, which is the current value (if any). */
    current_value = token[eTokenValue];
  }

  if (eOK == status) {
    /* Get next token, optional trailing quote. */
    trail_quote = token[eTokenTrailingQuote];
    if (0 == trail_quote) {
      /* No trailing quote, so create one. */
      trail_quote = (char *) calloc(2, sizeof(char));
      if (0 == trail_quote) status = eDynAllocFailed;
      *trail_quote = '"';
    }
  }

  if (eOK == status) {
    /* Get next token, optional quote followed by trailing space. */
    suffix = token[eTokenSuffix];

    if (0 == suffix) {
      /* No suffix, so create blank one. */
      suffix = (char *) calloc(1, sizeof(char));
      if (0 == suffix) status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* If value is between unescaped matching quotes, it's only necessary to escape those. */
    if (('\'' == *lead_quote || '"' == *lead_quote) && *lead_quote == *trail_quote) delim = lead_quote;
  }

  if (eOK == status) {

    /* Special case: if backslash is the last non-white space character, keep one more character
       because the backslash escaped it. */
    if ('\0' != *end && '\\' == *(end - 1)) ++end;

    /* Escape sensitive characters. */
    status = expand_escape(begin, end, delim, &value);
  }

  if (eOK == status) {
    /* Create a new list, including room for terminating 0. */
    new_token = (char **) calloc(eNumTokens, sizeof(char *));
    if (0 == new_token) status = eDynAllocFailed;
  }

  if (eOK == status) {
    new_token[eTokenPrefix] = prefix;
    new_token[eTokenLeadingQuote] = lead_quote;
    new_token[eTokenValue] = value;
    new_token[eTokenTrailingQuote] = trail_quote;
    new_token[eTokenSuffix] = suffix;
  }

  if (eOK == status) {
    /* Successful assignment; put new token into list at position given by field_iterator. */
    ape_list_set(field_itor, new_token);
    free(token); token = 0;
    free(current_value); current_value = 0;
  } else {
    /* Problem with assignment; clean up from the attempt. */
    /* TODO: prefix, suffix, lead_quote, trail_quote are potential memory leaks. Are they even necessary? */
    free(new_token); new_token = 0;
    free(value); current_value = 0;
  }

  return status;
}

/* Look for the first occurance in the string of any of a set of delimiters, but ignore characters within a
   quoted string (single or double quoted). */
/* TODO: Test this exhaustively. */
static const char * find_unquoted(const char * s, const char * delim) {
  if (0 != s && 0 != delim) {
    char double_quote = 0;
    char single_quote = 0;
    char quote = 0;

    for (; '\0' != *s; ++s) {
      /* Backslash escapes the next character no matter what, whether inside or outside a quote. */
      if ('\\' == *s && '\0' != *(s + 1)) {
        ++s;
      } else if (0 == quote) {
        /* Not currently inside a quoted string, so see if this is an opening quote. */
        if ('\'' == *s) {
          quote = 1;
          single_quote = 1;
        } else if ('"' == *s) {
          quote = 1;
          double_quote = 1;
        } else {
          /* Not a quote: check if this is one of the delimiters. */
          const char * d = delim;
          for (; '\0' != *d && *d != *s; ++d) {}
          if ('\0' != *d) break; /* This is the first found delimiter. */
        }
      } else {
        /* Currently inside a quoted string, so see if this is the closing quote. */
        if ('\'' == *s && 0 != single_quote) {
          quote = 0;
          single_quote = 0;
        } else if ('"' == *s && 0 != double_quote) {
          quote = 0;
          double_quote = 0;
        }
      }
    }
  }
  return s;
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_par.c,v $
 * Revision 1.101  2011/02/18 19:34:52  irby
 * Clarify that eInvalidName refers to a bad entry in the par file.
 *
 * Revision 1.100  2011/02/01 18:01:47  jpeachey
 * Add support and tests for functions to convert string values into int
 * and short types, and to convert a parameter value into short type.
 *
 * Revision 1.99  2011/01/21 21:33:13  irby
 * Strip leading/trailing white space only in ape_par_set_value_string, not in
 * other related functions (ape_par_set_field).  This fixes a bug in which a
 * default parameter value of " " (single space) was truncated when the file
 * was saved (ftlist).
 *
 * Revision 1.98  2010/11/23 22:05:55  jpeachey
 * Replace tabs with spaces to make code line up correctly.
 *
 * Revision 1.97  2010/11/23 19:22:21  jpeachey
 * Move read_from_stdin static function from ape_par to ape_util. This is
 * only used in the readline-free version.
 *
 * Revision 1.96  2010/11/12 20:53:03  irby
 * Moved some internal static functions from ape_par to ape_util, and allowed
 * new custom_get_text routine to be used instead of ape's standard getter in
 * order to handle xpi call-back.
 *
 * Revision 1.95  2010/01/19 17:53:37  peachey
 * Be more careful about handling and reporting values that are interpreted
 * as undefined or infinite.
 *
 * Revision 1.94  2009/07/09 21:19:05  peachey
 * Change conditional compilation so that readline clean-up
 * function is not compiled at all when readline is not used.
 *
 * Revision 1.93  2009/07/08 19:13:30  peachey
 * Work in progress to clean up after readline allocations.
 *
 * Revision 1.92  2007/11/16 18:20:43  peachey
 * Refactor ape_par_check and associated helper functions. Check more error
 * conditions and make the messages more specific. Streamline and simplify
 * the code for displaying information about the code. Decouple the
 * reporting of errors from code which determines whether an error is
 * related to a user input, and hence may be recoverable.
 *
 * Revision 1.91  2007/11/12 20:02:35  peachey
 * Remove unneeded extra call to ape_par_check from ape_par_get_string_case,
 * and uncomment some tests that didn't use to pass but now do.
 *
 * Revision 1.90  2007/11/12 19:27:18  peachey
 * Be more careful about calling ape_par_check from within the ape_par*
 * functions. Always call ape_par_check regardless of status, but only
 * report the status from ape_par_check if the status was not already non-0.
 *
 * Revision 1.89  2007/10/09 18:27:58  peachey
 * Add checks for various file modes (enrw).
 *
 * Revision 1.88  2007/10/09 16:43:45  peachey
 * Add function ape_par_redirect_prompt_stream to allow prompts to be
 * redirected. This also encapsulates all readline usage to within the ape_par module.
 *
 * Revision 1.87  2007/08/21 19:33:55  peachey
 * Add ape_par_get_name function.
 *
 * Revision 1.86  2007/07/26 16:22:32  peachey
 * Remove block that defines USE_READLINE, since that #define is now
 * provided by the Makefile.
 *
 * Revision 1.85  2007/05/15 20:52:48  peachey
 * Additional speed optimizations:
 * o Reduce number of calls to ape_list_end by calling it once and storing
 *   the returned iterator for use in loop logic when it is safe to do so.
 * o Make ape_list_destroy more efficient by iterating through the list
 *   back to front once and freeing the elements instead of repeatedly
 *   calling ape_list_remove_entry.
 *
 * Revision 1.84  2007/02/20 17:14:18  peachey
 * Add multiple entries to enumerated test.
 *
 * Revision 1.83  2007/02/15 21:00:17  peachey
 * Add ape_par_get/set_comment for getting and setting the
 * comment field from a parameter.
 *
 * Revision 1.82  2007/02/01 19:02:26  peachey
 * Streamline and refactor slightly the parameter checking code.
 *
 * Revision 1.81  2007/01/25 17:49:57  peachey
 * Fixed bug in ape_par_get_string_case that caused it not to behave
 * correctly when called for a non-enumerated parameter. Added functions
 * ape_trad_prompt_string_case and ape_trad_get_string_case for getting
 * parameters with control over case of returned string.
 *
 * Revision 1.80  2006/12/21 21:43:43  peachey
 * Add ape_par_get_string_case, similar to ape_par_get_string but more
 * flexible with regard to converting case of returned string values. Change
 * ape_par_get_string to call ape_par_get_string_case internally.
 *
 * Revision 1.79  2006/11/30 20:38:35  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.78  2006/11/24 19:51:27  peachey
 * Add ape_par_get_float.
 *
 * Revision 1.77  2006/11/24 18:10:57  peachey
 * Add a cast to make Visual Studio happy.
 *
 * Revision 1.76  2006/11/24 15:55:26  peachey
 * Add command line history within the session.
 *
 * Revision 1.75  2006/11/07 15:38:27  peachey
 * Enable re-prompting if a prompt is issued and user enters bad input.
 *
 * Revision 1.74  2006/11/03 14:29:15  peachey
 * Enable reprompting if errors occur in the first prompt.
 *
 * Revision 1.73  2006/11/01 18:57:09  peachey
 * Add support for fe and fn file types, as well as all blends of the
 * file modifiers, e, n, r, w, en, er etc.
 *
 * Revision 1.72  2006/10/30 21:55:31  peachey
 * Add tolerance for the g parameter type, as needed by xselect.
 * Treat it as a string, but do not consider it an error.
 *
 * Revision 1.71  2006/08/25 19:45:20  peachey
 * Change check_as helper function for check_range so that it
 * does not consider undefined or infinite numerical parameters as
 * an error.
 *
 * Revision 1.70  2006/06/23 20:37:51  peachey
 * Interpret output from ape_util_check_file_access according to new standard.
 *
 * Revision 1.69  2006/06/23 20:19:07  peachey
 * Include file types in those quoted by ape_par_set_value_string.
 *
 * Revision 1.68  2006/06/23 20:09:16  peachey
 * For string parameters, use quote_field to rationalize values which are
 * assigned using ape_par_set_value_string.
 *
 * Revision 1.67  2006/06/22 21:00:03  peachey
 * Uppercase enumerated string parameters.
 *
 * Revision 1.66  2006/06/22 20:25:41  peachey
 * Use $ENV{VAR_NAME} instead of simply $VAR_NAME for environment variables.
 *
 * Revision 1.65  2006/06/22 20:04:08  peachey
 * Add and test quote_field, for adding outside double quotes
 * correctly to a string, regardless of how it was created.
 *
 * Revision 1.64  2006/06/22 15:39:00  peachey
 * Change extract_field to find_field, which returns whole field as a
 * vector of tokens.
 *
 * Revision 1.63  2006/06/22 15:08:07  peachey
 * Save quotes separately from white space.
 *
 * Revision 1.62  2006/06/22 14:32:59  peachey
 * Add more torture tests for escape expansion.
 *
 * Revision 1.61  2006/06/22 13:51:38  peachey
 * Make collapse_escape more selective about which characters it escapes.
 *
 * Revision 1.60  2006/06/21 15:57:13  peachey
 * Add enum for terminating null in token lists.
 *
 * Revision 1.59  2006/06/21 14:54:10  peachey
 * Change internal storage of parameter fields to use a simple array of
 * three strings, prefix (white space + quote), value, suffix (quote + white space).
 *
 * Revision 1.58  2006/06/20 19:26:38  peachey
 * Test case of parsing empty double-quoted string, fix bug which was
 * causing the quotes not to be recognized.
 *
 * Revision 1.57  2006/06/20 17:22:28  peachey
 * Expand environment variables when getting values from parameters.
 *
 * Revision 1.56  2006/06/20 02:53:24  peachey
 * Rework get_field/set_field to unescape/escape characters when interpreting fields.
 *
 * Revision 1.55  2006/06/19 20:55:54  peachey
 * Change arguments of expand_escape to accept an input range instead of a
 * assumed 0-terminated string.
 *
 * Revision 1.54  2006/06/18 01:11:43  peachey
 * Complete and test expand_escape, for adding escape (backslash)
 * where needed to prevent quotes from being interpreted as quotes.
 *
 * Revision 1.53  2006/06/17 18:13:09  peachey
 * Add collapse_escape, for removing excess backslashes from a string.
 *
 * Revision 1.52  2006/06/17 03:44:51  peachey
 * Correct strings describing some tests of ape_par_set_field.
 *
 * Revision 1.51  2006/06/17 03:32:17  peachey
 * Clean up unused static functions, correct some comments.
 *
 * Revision 1.50  2006/06/17 03:19:17  peachey
 * Do not add any fields for a completely blank line.
 *
 * Revision 1.49  2006/06/17 03:01:30  peachey
 * Refactor parameter parsing code to be simpler and more robust with respect
 * to quotes.
 *
 * Revision 1.48  2006/06/15 19:45:00  peachey
 * Handle # embedded within fields, even without quotes.
 *
 * Revision 1.47  2006/06/15 13:40:56  peachey
 * Change signature of find_unquoted to work with constant strings.
 *
 * Revision 1.46  2006/06/14 16:03:51  peachey
 * Make prompts more informative in prompt tests.
 *
 * Revision 1.45  2006/06/13 17:16:48  peachey
 * Use int instead of char for comparison, to avoid portability problems.
 *
 * Revision 1.44  2006/06/13 14:51:44  peachey
 * If string includes unquoted comma, add quotes.
 *
 * Revision 1.43  2006/06/09 04:09:45  peachey
 * Check parameter for validity right after prompting.
 *
 * Revision 1.42  2006/06/07 06:06:47  peachey
 * Reverse the output of ape_util_file_check_access: 1 means access is OK,
 * 0 means it isn't.
 *
 * Revision 1.41  2006/06/07 02:44:03  peachey
 * Check file accessibility as part of a value check.
 *
 * Revision 1.40  2006/06/06 13:28:07  peachey
 * Replace ape_par_get/set_prompt_mode with ape_par_get/set_prompt_style
 * and ape_par_get/set_default_prompt_style. Support no-prompt, query-for-hidden
 * and multi-query styles.
 *
 * Revision 1.39  2006/06/03 01:52:22  peachey
 * Reenable readline, which was switched off by accident.
 *
 * Revision 1.38  2006/05/31 18:10:09  peachey
 * Add ape_par_get_prompt_mode/ape_par_set_prompt_mode.
 *
 * Revision 1.37  2006/05/31 03:21:44  peachey
 * Rationalize ape_par_get_* family of functions, and use them to implement
 * ape_trad_get_* and ape_trad_query_*.
 *
 * Revision 1.36  2006/05/31 01:41:36  peachey
 * Rename ape_par_set_string to ape_par_set_field.
 *
 * Revision 1.35  2006/05/31 01:36:29  peachey
 * Rename ape_par_get_string to ape_par_get_field and ape_par_get_string_array to
 * ape_par_get_field_array.
 *
 * Revision 1.34  2006/05/19 17:30:49  peachey
 * Update TODO which was done.
 *
 * Revision 1.33  2006/05/18 15:27:25  peachey
 * Test using single quotes for prompt to allow double quotes inside prompt.
 *
 * Revision 1.32  2006/05/18 15:01:11  peachey
 * Add support for single quotes instead of double quotes.
 *
 * Revision 1.31  2006/05/18 03:17:59  peachey
 * Improve error reporting in debug statements in
 * ape_par_check.
 *
 * Revision 1.30  2006/05/12 03:31:07  peachey
 * Correct bugs in tests for ape_par_query, and a bug in ape_par_query.
 *
 * Revision 1.29  2006/05/12 00:20:25  peachey
 * Add and test ape_par_query.
 *
 * Revision 1.28  2006/05/11 20:14:52  peachey
 * Add ape_set_value_string, for setting value and flagging it as modified.
 * Change ape_par_get_eff_mode to use this flag to make the effective mode hidden
 * after a parameter was set either on the command line or by a prompt.
 *
 * Revision 1.27  2006/05/11 19:45:48  peachey
 * Add ape_par_get_eff_mode, for getting effective mode of parameter,
 * taking into account mode parameter and command line arguments.
 *
 * Revision 1.26  2006/05/10 16:48:50  peachey
 * Generalize format tests to allow them to be run with or without
 * value checking, to let them be used both in contexts where the whole parameter
 * needs to be checked (checking user input) and contexts in which a bad value
 * is not a show-stopper, because the user or client may yet correct the problem.
 *
 * Revision 1.25  2006/05/10 13:56:51  peachey
 * Patch memory leak.
 *
 * Revision 1.24  2006/05/10 04:44:56  peachey
 * Include range, if any, and current value in the prompt.
 *
 * Revision 1.23  2006/05/10 03:57:51  peachey
 * Add some extra comments on the whole parameter checking mechanism.
 *
 * Revision 1.22  2006/05/10 03:37:46  peachey
 * Cosmetic rewording of error message.
 *
 * Revision 1.21  2006/05/10 03:17:45  peachey
 * Generalize integer-specific range checker.
 *
 * Revision 1.20  2006/05/10 00:30:31  peachey
 * Rename eMin/MaxTypeMismatch to the more generic eMin/MaxConversionError.
 *
 * Revision 1.19  2006/05/10 00:25:50  peachey
 * If an error occurs when converting value string, just return the error as is.
 *
 * Revision 1.18  2006/05/09 19:24:26  peachey
 * Add ape_par_check function, for determining whether parameter format is valid.
 *
 * Revision 1.17  2006/05/05 01:07:44  peachey
 * Make get_field more safe, by returning int code and using a passed
 * pointer for the returned string.
 *
 * Revision 1.16  2006/05/03 02:06:28  peachey
 * Remove echoed newline, which it turns out only appears to be needed.
 *
 * Revision 1.15  2006/05/03 01:58:49  peachey
 * Echo newline outside the ifdef so that it gets echoed either way.
 *
 * Revision 1.14  2006/05/03 01:29:13  peachey
 * Add ape_par_prompt, which uses readline on unix and homemade code on windows
 * to get a string from the user.
 *
 * Revision 1.13  2006/04/28 01:24:04  peachey
 * Add a couple casts to get the constness right when calling free.
 * Enable al mode, and remove eFileOpenError, using eFileReadError in lieu thereof.
 * The latter change was to make error messages consistent between Windows and
 * Unix, because although the two both fail to open a directory as if it were
 * a file, they fail at different points in the process of opening the file.
 *
 * Revision 1.12  2006/04/26 14:46:10  peachey
 * Plug memory leak.
 *
 * Revision 1.11  2006/04/26 14:40:55  peachey
 * Add ape_par_get_enum_string and ape_par_get_min_string for interpreting
 * the min field of the parameter.
 *
 * Revision 1.10  2006/04/26 01:34:25  peachey
 * Add ape_par_get_string_array, for getting arrays of strings from
 * fields in parameters. Add helper functions for this: extract_field
 * (derived from get_field, which was refactored to use the new function
 * as well), and get_field_array, which parallels get_field.
 *
 * Revision 1.9  2006/04/24 21:27:27  peachey
 * Add ape_par_get_mode, reusing most of ape_par_get_type in a static function,
 * encode_msg.
 *
 * Revision 1.8  2006/04/24 17:43:40  peachey
 * Add implementation and tests for ape_par_get_type.
 *
 * Revision 1.7  2006/04/22 01:36:56  peachey
 * Add handling of comments.
 *
 * Revision 1.6  2006/04/20 04:46:47  peachey
 * Add ape_par_set_string function, for setting arbitrary single fields in
 * a parameter.
 *
 * Revision 1.5  2006/04/20 03:43:59  peachey
 * Add and test ape_par_clone.
 *
 * Revision 1.4  2006/04/12 03:29:13  peachey
 * Split parse_line into two parts: parse_line + parse_field which parses
 * one individual field (list of tokens).
 *
 * Revision 1.3  2006/04/11 03:27:00  peachey
 * Add ape_par_get_double and its unit test.
 *
 * Revision 1.2  2006/04/11 01:54:56  peachey
 * Publish field identifier enum as part of the interface. Add and test
 * ape_par_get_string, for getting any field as a string.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
