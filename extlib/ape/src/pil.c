/** \file pil.c
    \brief Implementation of Ape's PIL emulating interface.

    Ape is backward compatible with a subset of the PIL library developed by
    ISDC. Replacements are provided for the most commonly used API functions.
    Less common functions and functions which rely on details of pil's internal
    structures are not provided. Ape's implementation of these functions does
    not always behave identically to PIL's implementation, but the behavior should be
    close enough to provide compatibility.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_test.h"
#include "ape/ape_trad.h"
#include "ape/ape_util.h"

#include "pil.h"
#include "pil_error.h"

#ifndef APE_COMPILE_FORTRAN_BINDINGS
#define APE_COMPILE_FORTRAN_BINDINGS
#endif

#ifdef APE_COMPILE_FORTRAN_BINDINGS
#include "cfortran.h"
#endif

#include <assert.h>
#include <float.h>
#include <stddef.h>
#include <string.h>

#ifdef WIN32
#define QC ";"
#define QS "|"
#define QD "\\"
#define QD2 "/"
#else
#define QC ":"
#define QS ";"
#define QD "/"
#define QD2 "/"
#endif

static int convert_ape_status_to_pil(int ape_status) {
  typedef struct lookup_struct {
    int ape_status;
    int pil_status;
  } lookup_struct;

  static const lookup_struct s_lookup[] = {
    { eOK, PIL_OK },
    { eNullPointer, PIL_NUL_PTR },
    { eDynAllocFailed, PIL_NO_MEM },
    { eTooManyAtExitFunc, PIL_UNSPECIFIED_ERROR },
    { eAtExitError, PIL_UNSPECIFIED_ERROR },
    /* 5 */
    { eInvalidArgument, PIL_BAD_ARG },
    { eFileNotFound, PIL_FILE_NOT_EXIST },
    { eFileReadError, PIL_ERR_FREAD },
    { eFileWriteError, PIL_ERR_FWRITE },
    { eFileRenameError, PIL_ERR_FWRITE },
    /* 10 */
    { eListError, PIL_UNSPECIFIED_ERROR },
    { eStringRemainder, PIL_BAD_VALUE },
    { eTypeMismatch, PIL_BAD_VALUE },
    { eOverflow, PIL_BAD_VALUE },
    { eUnderflow, PIL_BAD_VALUE },
    /* 15 */
    { eConversionError, PIL_BAD_VALUE },
    { eVarNotSet, PIL_UNSPECIFIED_ERROR },
    { eLineTooLong, PIL_LINE_ERROR },
    { eTooManyArguments, PIL_LINE_TOO_MANY },
    { eParameterDuplicated, PIL_BOGUS_CMDLINE },
    /* 20 */
    { eParNotFound, PIL_NOT_FOUND },
    { eUnnamedPar, PIL_BAD_NAME },
    { eInvalidParName, PIL_BAD_NAME },
    { eFieldNotFound, PIL_BAD_LINE },
    { eNewParameter, PIL_UNSPECIFIED_ERROR },
    /* 25 */
    { eUnknownMode, PIL_BAD_MODE },
    { eUnknownType, PIL_BAD_TYPE },
    { eFormatError, PIL_BAD_LINE },
    { eRangeEnum, PIL_UNSPECIFIED_ERROR },
    { eRangeNoEnum, PIL_UNSPECIFIED_ERROR },
    /* 30 */
    { eInputFailure, PIL_ERR_FREAD },
    { eTooManyFields, PIL_BAD_LINE },
    { eMinConversionError, PIL_BAD_VALUE },
    { eMaxConversionError, PIL_BAD_VALUE },
    { eValueBelowMin, PIL_OFF_RANGE },
    /* 35 */
    { eValueAboveMax, PIL_OFF_RANGE },
    { eInvalidRange, PIL_OFF_RANGE },
    { eInvalidChoice, PIL_BAD_ENUM_VALUE },
    { eInvalidName, PIL_BAD_NAME },
    { eNoParLoaded, PIL_UNSPECIFIED_ERROR },
    /* 40 */
    { eUninitialized, PIL_UNSPECIFIED_ERROR },
    { eUndefinedValue, PIL_VALUE_UNDEFINED },
    { eNan, PIL_BAD_VAL_REAL },
    { eFileNotAccessible, PIL_BAD_FILE_ACCESS },
    { eNoModePar, PIL_UNSPECIFIED_ERROR },
    /* 45 */
    { eAmbiguousParName, PIL_BAD_NAME }
  };
  static const size_t s_num_code = sizeof(s_lookup) / sizeof(s_lookup[0]);

  size_t idx = 0;
  int pil_status = PIL_UNSPECIFIED_ERROR;
  for (idx = 0; idx != s_num_code; ++idx) {
    if (s_lookup[idx].ape_status == ape_status) {
      pil_status = s_lookup[idx].pil_status;
      break;
    }
  }
  return pil_status;
}

#ifdef __cplusplus
extern "C" {
#endif

int PILInit(int argc, char ** argv) {
  return convert_ape_status_to_pil(ape_trad_init(argc, argv));
}

int PILClose(int status) {
  return convert_ape_status_to_pil(ape_trad_close(0 == status ? 1 : 0));
}

int PILGetBool(const char * name, int * result) {
  char tmp_result = 0;
  int status = 0;

  if (0 != result) {
    *result = 0;
    status = ape_trad_query_bool(name, &tmp_result);
    *result = tmp_result;
  }

  return convert_ape_status_to_pil(status);
}

int PILGetInt(const char * name, int * result) {
  long tmp_result = 0;
  int status = 0;

  if (0 != result) {
    *result = 0;
    status = ape_trad_query_long(name, &tmp_result);
    *result = tmp_result;
  }

  return convert_ape_status_to_pil(status);
}

int PILGetLong(const char * name, long * result) {
  return convert_ape_status_to_pil(ape_trad_query_long(name, result));
}

int PILGetReal(const char * name, double * result) {
  return convert_ape_status_to_pil(ape_trad_query_double(name, result));
}

int PILGetReal4(const char * name, float * result) {
  return convert_ape_status_to_pil(ape_trad_query_float(name, result));
}

int PILGetString(const char * name, char * result) {
  char * tmp_result = 0;
  int status = 0;

  if (0 != result) {
    *result = '\0';
    status = ape_trad_query_string(name, &tmp_result);
    if (0 != tmp_result) {
      size_t tmp_len = PIL_LINESIZE < strlen(tmp_result) + 1 ? PIL_LINESIZE : strlen(tmp_result) + 1;
      strncpy(result, tmp_result, tmp_len);
    }
    free(tmp_result); tmp_result = 0;
  }

  return convert_ape_status_to_pil(status);
}

int PILGetFname(const char * name, char * result) {
  char * tmp_result = 0;
  int status = 0;

  if (0 != result) {
    *result = '\0';
    status = ape_trad_query_file_name(name, &tmp_result);
    if (0 != tmp_result) {
      size_t tmp_len = PIL_LINESIZE < strlen(tmp_result) + 1 ? PIL_LINESIZE : strlen(tmp_result) + 1;
      strncpy(result, tmp_result, tmp_len);
    }
    free(tmp_result); tmp_result = 0;
  }

  return convert_ape_status_to_pil(status);
}

int PILGetAsString(const char * name, char * result) {
  return convert_ape_status_to_pil(PILGetString(name, result));
}

int PILPutBool(const char * name, int value) {
  return convert_ape_status_to_pil(ape_trad_set_bool(name, 0 != value ? 1 : 0));
}

int PILPutInt(const char * name, int value) {
  return convert_ape_status_to_pil(ape_trad_set_long(name, value));
}

int PILPutReal(const char * name, double value) {
  return convert_ape_status_to_pil(ape_trad_set_double(name, value));
}

int PILPutString(const char * name, const char * value) {
  return convert_ape_status_to_pil(ape_trad_set_string(name, value));
}

int PILPutFname(const char * name, const char * value) {
  return convert_ape_status_to_pil(ape_trad_set_file_name(name, value));
}

int PILGetQueryMode(void) {
  int prompt_style = eDefaultPrompt;
  ape_par_get_default_prompt_style(&prompt_style);
  return 0 != (eNoPrompt & prompt_style) ? PIL_QUERY_OVERRIDE : PIL_QUERY_DEFAULT;
}

int PILOverrideQueryMode(int new_mode) {
  int status = PIL_OK;
  int prompt_style = eDefaultPrompt;
  int ape_query_override = eNoPrompt;

  if (eOK == status) {
    /* Get the current default prompt style. */
    ape_par_get_default_prompt_style(&prompt_style);
  }

  if (eOK == status) {
    /* Map PIL override query mode to ape's version. */
    switch (new_mode) {
      case PIL_QUERY_DEFAULT:
        prompt_style &= ~ape_query_override;
        break;
      case PIL_QUERY_OVERRIDE:
        prompt_style |= ape_query_override;
        break;
      default:
        status = eInvalidArgument;
        break;
    }
  }

  if (eOK == status) {
    /* Change ape's default prompt style. */
    status = ape_par_set_default_prompt_style(prompt_style);
  }

  return convert_ape_status_to_pil(status);
}

static int (*s_out_func)(char *) = 0;

static void out_writer(const char * msg) {
  if (0 != s_out_func) (*s_out_func)((char *) msg);
}

int PILSetLoggerFunction(int (*func)(char *)) {
  s_out_func = func;
  ape_msg_set_out_handler(&out_writer);
  return PIL_OK;
}

int PILSetFileAccessFunction(int (*func)(const char *, const char *)) {
  return convert_ape_status_to_pil(ape_util_set_file_check_func(func));
}

int PILSetReprompt(const char * par_name, int reprompt) {
  int status = eOK;
  ApePar * par = 0;
  int prompt_style = eDefaultPrompt;
  int ape_reprompt = eMultiQuery | eQueryHidden; /* PIL reprompt is the same as ape's query multiple times + query for hidden. */

  /* Find the parameter. */
  status = ape_trad_find_par(par_name, &par);

  if (eOK == status) {
    /* Get the current prompt style. */
    ape_par_get_prompt_style(par, &prompt_style);
  }

  if (eOK == status) {
    /* Map PIL reprompt to ape_reprompt flags. */
    if (0 != reprompt) prompt_style |= ape_reprompt;
    else prompt_style &= ~ape_reprompt;

    /* Change ape's prompt style. */
    status = ape_par_set_prompt_style(par, prompt_style);
  }

  return convert_ape_status_to_pil(status);
}

int PILFlushParameters(void) {
  return convert_ape_status_to_pil(ape_trad_save());
}

#ifdef APE_COMPILE_FORTRAN_BINDINGS

#define pilinit_STRV_A2 NUM_ELEM_ARG(1)
FCALLSCFUN2(INT, PILInit, PILINIT, pilinit, INT, PSTRINGV)

FCALLSCFUN1(INT, PILClose, PILCLOSE, pilclose, INT)

FCALLSCFUN2(INT, PILGetBool, PILGETBOOL, pilgetbool, STRING, INTV)

FCALLSCFUN2(INT, PILGetInt, PILGETINT, pilgetint, STRING, INTV)

FCALLSCFUN2(INT, PILGetReal, PILGETREAL, pilgetreal, STRING, DOUBLEV)

FCALLSCFUN2(INT, PILGetReal4, PILGETREAL4, pilgetreal4, STRING, FLOATV)

FCALLSCFUN2(INT, PILGetString, PILGETSTRING, pilgetstring, STRING, PSTRING)

FCALLSCFUN2(INT, PILGetFname, PILGETFNAME, pilgetfname, STRING, PSTRING)

FCALLSCFUN2(INT, PILGetAsString, PILGETASSTRING, pilgetasstring, STRING, PSTRING)

FCALLSCFUN2(INT, PILPutBool, PILPUTBOOL, pilputbool, STRING, INT)

FCALLSCFUN2(INT, PILPutInt, PILPUTINT, pilputint, STRING, INT)

FCALLSCFUN2(INT, PILPutReal, PILPUTREAL, pilputreal, STRING, DOUBLE)

FCALLSCFUN2(INT, PILPutString, PILPUTSTRING, pilputstring, STRING, STRING)

FCALLSCFUN2(INT, PILPutFname, PILPUTFNAME, pilputfname, STRING, STRING)

#endif

void ape_pil_test(void) {
  int status = eOK;
  char * pfiles = 0;
  status = ape_io_get_pfiles(&pfiles);
  if (eOK == status) {
    status = ape_io_set_pfiles(QS".");
  }
  if (eOK == status) {
    char * argv[] = { "ape_test", "sh=sh" };
    int argc = sizeof(argv) / sizeof(argv[0]);
    char value[PIL_LINESIZE] = "";
    const char * init_value = 0;

    status = PILInit(argc, argv);
    ape_test_cmp_long("PILInit(ape_test sh=sh)", 0, 0, status, PIL_OK);

    /* Should prompt. */
    init_value = "0: You should see this; enter \"sq\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sq", value);
    /* For cosmetic reasons, add a newline here. */
    ape_msg_debug("\n");
    ape_test_cmp_string("After PILInit(ape_test sh=sh), PILGetString(\"sq\")", value, "sq", status, PIL_OK);

    /* Should not prompt. */
    init_value = "1: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sh", init_value);
    ape_test_cmp_long("Setting sh=sh", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sh", value);
    ape_test_cmp_string("After PILInit(ape_test sh=sh), PILGetString(\"sh\")", value, init_value, status, PIL_OK);

    /* Suppress prompts with override query mode. */
    status = PILOverrideQueryMode(PIL_QUERY_OVERRIDE);
    ape_test_cmp_long("PILOverrideQueryMode(PIL_QUERY_OVERRIDE)", 0, 0, status, PIL_OK);

    /* Should not prompt. */
    init_value = "2: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sq", value);
    ape_test_cmp_string("After PILOverrideQueryMode, PILGetString(\"sq\")", value, init_value, status, PIL_OK);

    init_value = "2.5: You should not see this$ENV{NON_EXISTENT_ENV_VAR}; enter \"wrong\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    init_value = "2.5: You should not see this; enter \"wrong\" if you do";
    status = PILGetString("sq", value);
    ape_test_cmp_string("After PILOverrideQueryMode, PILGetString(\"sq\")", value, init_value, status, PIL_OK);

    /* Should not prompt. */
    init_value = "3: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sh", init_value);
    ape_test_cmp_long("Setting sh=sh", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sh", value);
    ape_test_cmp_string("After PILOverrideQueryMode, PILGetString(\"sh\")", value, init_value, status, PIL_OK);

    /* Set reprompt mode to allow infinite reprompting. This should not have an immediate effect because all prompts disabled. */
    status = PILSetReprompt("sq", 1);
    ape_test_cmp_long("PILSetReprompt(sq, 1)", 0, 0, status, PIL_OK);

    /* Should not prompt. */
    init_value = "4: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sq", value);
    ape_test_cmp_string("After PILSetReprompt, PILGetString(\"sq\")", value, init_value, status, PIL_OK);

    /* Should not prompt. */
    init_value = "5: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sh", init_value);
    ape_test_cmp_long("Setting sh=sh", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sh", value);
    ape_test_cmp_string("After PILSetReprompt, PILGetString(\"sh\")", value, init_value, status, PIL_OK);

    /* Enable prompts by switching off override query mode. */
    status = PILOverrideQueryMode(PIL_QUERY_DEFAULT);
    ape_test_cmp_long("PILOverrideQueryMode(PIL_QUERY_DEFAULT)", 0, 0, status, PIL_OK);

    /* Should prompt. */
    init_value = "6: You should see this; enter \"sq\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sq", value);
    /* For cosmetic reasons, add a newline here. */
    ape_msg_debug("\n");
    ape_test_cmp_string("After PILOverrideQueryMode(PIL_QUERY_DEFAULT), PILGetString(\"sq\")", value, "sq", status, PIL_OK);

    /* Should not prompt. */
    init_value = "7: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sh", init_value);
    ape_test_cmp_long("Setting sh=sh", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sh", value);
    ape_test_cmp_string("After PILOverrideQueryMode(PIL_QUERY_DEFAULT), PILGetString(\"sh\")", value, init_value, status, PIL_OK);

    /* Don't save anything, but do clean up. */
    PILClose(-1);
    status = PIL_OK;
  }
  if (PIL_OK == status) {
    char * argv[] = { "ape_test", "sq=sq", "sh=sh" };
    int argc = sizeof(argv) / sizeof(argv[0]);
    char value[PIL_LINESIZE] = "";
    const char * init_value = 0;

    status = PILInit(argc, argv);
    ape_test_cmp_long("PILInit(ape_test sq=sq sh=sh)", 0, 0, status, PIL_OK);

    /* Should not prompt, since sq set on command line. */
    init_value = "8: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sq", value);
    ape_test_cmp_string("After PILInit(ape_test sq=sq sh=sh), PILGetString(\"sq\")", value, init_value, status, PIL_OK);

    /* Should not prompt. */
    init_value = "9: You should not see this; enter \"wrong\" if you do";
    status = PILPutString("sh", init_value);
    ape_test_cmp_long("Setting sh=sh", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sh", value);
    ape_test_cmp_string("After PILInit(ape_test sq=sq sh=sh), PILGetString(\"sh\")", value, init_value, status, PIL_OK);

    /* Set reprompt mode to allow infinite reprompting. This should not have an immediate effect because all prompts disabled. */
    status = PILSetReprompt("sq", 1);
    ape_test_cmp_long("PILSetReprompt(sq, 1)", 0, 0, status, PIL_OK);

    status = PILSetReprompt("sh", 1);
    ape_test_cmp_long("PILSetReprompt(sh, 1)", 0, 0, status, PIL_OK);

    /* Should prompt, since reprompt mode set. */
    init_value = "10: You should see this; enter \"sq\" if you do";
    status = PILPutString("sq", init_value);
    ape_test_cmp_long("Setting sq=sq", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sq", value);
    /* For cosmetic reasons, add a newline here. */
    ape_msg_debug("\n");
    ape_test_cmp_string("After PILInit(ape_test sq=sq sh=sh), PILGetString(\"sq\")", value, "sq", status, PIL_OK);

    /* Should prompt, since reprompt mode set. */
    init_value = "11 You should see this; enter \"sh\" if you do";
    status = PILPutString("sh", init_value);
    ape_test_cmp_long("Setting sh=sh", 0, 0, status, PIL_OK);

    *value = '\0';
    status = PILGetString("sh", value);
    /* For cosmetic reasons, add a newline here. */
    ape_msg_debug("\n");
    ape_test_cmp_string("After PILInit(ape_test sq=sq sh=sh), PILGetString(\"sh\")", value, "sh", status, PIL_OK);

    /* Test an error case. */
    *value = '\0';
    status = PILGetString("non-existent", value);
    ape_test_cmp_string("PILGetString(\"non-existent\")", value, "", status, PIL_NOT_FOUND);

    /* Don't save anything, but do clean up. */
    PILClose(-1);
  }

  if (PIL_OK == status) {
    char * argv[] = { "non-existent" };
    int argc = sizeof(argv) / sizeof(argv[0]);

    status = PILInit(argc, argv);
    ape_test_cmp_long("PILInit(ape_test sq=sq sh=sh)", 0, 0, status, PIL_FILE_NOT_EXIST);
    PILClose(-1);
  }

  /* Final clean up. */
  ape_par_set_default_prompt_style(eDefaultPrompt);
  status = ape_io_set_pfiles(pfiles);
  free(pfiles); pfiles = 0;

  /* Test convert_ape_status_to_pil. */
  status = convert_ape_status_to_pil(eTooManyArguments);
  ape_test_cmp_long("convert_ape_status_to_pil(eTooManyArguments)", status, PIL_LINE_TOO_MANY, 0, 0);
}

int PILGetRealVector(const char * parameter, int length, double * output) {
  /* TODO: either really implement these vector functions are remove them completely. */
  assert(0);
  return 0;
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: pil.c,v $
 * Revision 1.18  2011/02/18 19:39:16  irby
 * Update the lookup table for converting Ape error codes to PIL codes.
 *
 * Revision 1.17  2009/12/04 19:15:26  peachey
 * Add PILGetLong function.
 *
 * Revision 1.16  2007/10/09 16:46:55  peachey
 * Cosmetically fix output of unit test to look OK when running
 * the unit test with a log file.
 *
 * Revision 1.15  2006/11/30 16:40:27  peachey
 * Fix bug in test code: pfiles needs to be set correctly for Windows test.
 *
 * Revision 1.14  2006/11/24 19:48:52  peachey
 * Use ape_trad_query_float, not ape_trad_query_double, in PILGetReal4.
 *
 * Revision 1.13  2006/06/22 20:25:41  peachey
 * Use $ENV{VAR_NAME} instead of simply $VAR_NAME for environment variables.
 *
 * Revision 1.12  2006/06/20 17:28:32  peachey
 * Add a test for environment variables.
 *
 * Revision 1.11  2006/06/07 05:54:39  peachey
 * Be more careful about strncpy, which pads output string.
 *
 * Revision 1.10  2006/06/07 04:42:10  peachey
 * Correct signatures for PILPut* fortran wrappers.
 *
 * Revision 1.9  2006/06/07 04:39:12  peachey
 * Add stub PILGetRealVector for backward compatibility.
 *
 * Revision 1.8  2006/06/07 02:19:50  peachey
 * Cast ape errors into pil's error range for pil wrappers.
 *
 */
