/** \file ape_trad.c
    \brief Implementation of traditional XPI/PIL-compliant interface.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_io.h"
#include "ape/ape_list.h"
#include "ape/ape_par.h"
#include "ape/ape_session.h"
#include "ape/ape_test.h"
#include "ape/ape_trad.h"
#include "ape/ape_util.h"

#include <stdio.h>
#include <stdlib.h>

#ifndef APE_COMPILE_FORTRAN_BINDINGS
#define APE_COMPILE_FORTRAN_BINDINGS
#endif

#ifdef APE_COMPILE_FORTRAN_BINDINGS
#include "cfortran.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
static const char s_dir_separator[] = { '/', '\\' };
#define QC ";"
#define QS "|"
#define QD "\\"
#define QD2 "/"
#else
static const char s_dir_separator[] = { '/' };
#define QC ":"
#define QS ";"
#define QD "/"
#define QD2 "/"
#endif

static ApeSession * s_default_session = 0;

static int ape_trad_check_init(void) {
  int status = eOK;
  if (0 == s_default_session) status = eUninitialized;
  return status;
}

static void ape_trad_atexit(void) {
  ape_trad_close(0);
}

int ape_trad_init(int argc, char ** argv) {
  ape_util_atexit(&ape_trad_atexit);
  return ape_session_init(&s_default_session, argc, argv);
}

int ape_trad_close(int save_flag) {
  int status = ape_session_close(s_default_session, save_flag);
  s_default_session = 0;
  return status;
}

int ape_trad_get_current(ApeParFile ** par_file) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_current(s_default_session, par_file);
  }
  return status;
}

int ape_trad_get_sys_pfile(ApeParFile ** par_file) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_sys_pfile(s_default_session, par_file);
  }
  return status;
}

int ape_trad_get_par_names(char *** par_name) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_par_names(s_default_session, par_name);
  }
  return status;
}

int ape_trad_get_type(const char * par_name, char * par_type) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_type(s_default_session, par_name, par_type);
  }
  return status;
}

int ape_trad_query_bool(const char * par_name, char * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_bool(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_double(const char * par_name, double * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_double(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_float(const char * par_name, float * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_float(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_file_name(const char * par_name, char ** value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_file_name(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_int(const char * par_name, int * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_int(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_long(const char * par_name, long * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_long(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_string(const char * par_name, char ** value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_string(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_query_string_case(const char * par_name, char ** value, char case_code) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_query_string_case(s_default_session, par_name, value, case_code);
  }
  return status;
}

int ape_trad_get_bool(const char * par_name, char * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_bool(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_double(const char * par_name, double * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_double(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_float(const char * par_name, float * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_float(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_file_name(const char * par_name, char ** value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_file_name(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_int(const char * par_name, int * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_int(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_long(const char * par_name, long * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_long(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_string(const char * par_name, char ** value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_string(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_get_string_case(const char * par_name, char ** value, char case_code) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_get_string_case(s_default_session, par_name, value, case_code);
  }
  return status;
}

int ape_trad_set_bool(const char * par_name, char value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_bool(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_double(const char * par_name, double value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_double(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_float(const char * par_name, float value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_float(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_file_name(const char * par_name, const char * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_file_name(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_int(const char * par_name, int value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_int(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_long(const char * par_name, long value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_long(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_short(const char * par_name, short value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_short(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_set_string(const char * par_name, const char * value) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_set_string(s_default_session, par_name, value);
  }
  return status;
}

int ape_trad_save(void) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_save(s_default_session);
  }
  return status;
}

int ape_trad_find_par(const char * par_name, ApePar ** par) {
  int status = ape_trad_check_init();
  if (eOK == status) {
    status = ape_session_find_par(s_default_session, par_name, par);
  }
  return status;
}

#ifdef APE_COMPILE_FORTRAN_BINDINGS

#define ape_trad_init_STRV_A2 NUM_ELEM_ARG(1)
FCALLSCFUN2(INT, ape_trad_init, APE_TRAD_INIT, ape_trad_init, INT, PSTRINGV)

FCALLSCFUN1(INT, ape_trad_close, APE_TRAD_CLOSE, ape_trad_close, INT)

FCALLSCFUN0(INT, ape_trad_save, APE_TRAD_SAVE, ape_trad_save)

/* ape_trad_get_bool par_value is char, but cfortran has no char, only byte: */
static int ape_trad_get_byte(const char * par_name, signed char * par_value) {
  char par_value_byte = 0;
  int status = eOK;
  if (0 == par_name || 0 == par_value) {
    status = eNullPointer;
  } else {
    status = ape_trad_get_bool(par_name, &par_value_byte);
    *par_value = par_value_byte;
  }
  return status;
}
FCALLSCFUN2(INT, ape_trad_get_byte, APE_TRAD_GET_BYTE, ape_trad_get_byte, STRING, BYTEV)

FCALLSCFUN2(INT, ape_trad_get_int, APE_TRAD_GET_INT, ape_trad_get_int, STRING, INTV)

FCALLSCFUN2(INT, ape_trad_get_double, APE_TRAD_GET_DOUBLE, ape_trad_get_double, STRING, DOUBLEV)

FCALLSCFUN2(INT, ape_trad_get_float, APE_TRAD_GET_FLOAT, ape_trad_get_float, STRING, FLOATV)

static int ape_trad_query_string_cf(const char * par_name, char * value) {
  int status = eOK;
  char * cf_value = 0;
  status = ape_trad_query_string(par_name, &cf_value);
  if (eOK == status) {
    strcpy(value, cf_value);
  }
  free(cf_value);
  return status;
}
FCALLSCFUN2(INT, ape_trad_query_string_cf, APE_TRAD_QUERY_STRING, ape_trad_query_string, STRING, PSTRING)

static int ape_trad_get_string_cf(const char * par_name, char * value) {
  int status = eOK;
  char * cf_value = 0;
  status = ape_trad_get_string(par_name, &cf_value);
  if (eOK == status) {
    strcpy(value, cf_value);
  }
  free(cf_value);
  return status;
}
FCALLSCFUN2(INT, ape_trad_get_string_cf, APE_TRAD_GET_STRING, ape_trad_get_string, STRING, STRING)

static int ape_trad_query_file_name_cf(const char * par_name, char * value) {
  int status = eOK;
  char * cf_value = 0;
  status = ape_trad_query_file_name(par_name, &cf_value);
  if (eOK == status) {
    strcpy(value, cf_value);
  }
  free(cf_value);
  return status;
}
FCALLSCFUN2(INT, ape_trad_query_file_name_cf, APE_TRAD_QUERY_FILE_NAME, ape_trad_query_file_name, STRING, PSTRING)

static int ape_trad_get_file_name_cf(const char * par_name, char * value) {
  int status = eOK;
  char * cf_value = 0;
  status = ape_trad_get_file_name(par_name, &cf_value);
  if (eOK == status) {
    strcpy(value, cf_value);
  }
  free(cf_value);
  return status;
}
FCALLSCFUN2(INT, ape_trad_get_file_name_cf, APE_TRAD_GET_FILE_NAME, ape_trad_get_file_name, STRING, STRING)


FCALLSCFUN2(INT, ape_trad_set_bool, APE_TRAD_SET_BOOL, ape_trad_set_bool, STRING, INT)

FCALLSCFUN2(INT, ape_trad_set_int, APE_TRAD_SET_INT, ape_trad_set_int, STRING, INT)

FCALLSCFUN2(INT, ape_trad_set_double, APE_TRAD_SET_DOUBLE, ape_trad_set_double, STRING, DOUBLE)

FCALLSCFUN2(INT, ape_trad_set_float, APE_TRAD_SET_FLOAT, ape_trad_set_float, STRING, FLOAT)

FCALLSCFUN2(INT, ape_trad_set_long, APE_TRAD_SET_LONG, ape_trad_set_long, STRING, LONG)

FCALLSCFUN2(INT, ape_trad_set_short, APE_TRAD_SET_SHORT, ape_trad_set_short, STRING, SHORT)

FCALLSCFUN2(INT, ape_trad_set_string, APE_TRAD_SET_STRING, ape_trad_set_string, STRING, STRING)

FCALLSCFUN2(INT, ape_trad_set_file_name, APE_TRAD_SET_FILE_NAME, ape_trad_set_file_name, STRING, STRING)

#endif

void ape_trad_test(void) {
  int status = eOK;

  /* Attempt to initialize a task using invalid inputs. */
  /* Null argc and argv. */
  { int argc = 0;
    char ** argv = 0;
    status = ape_trad_init(argc, argv);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_trad_init(0, 0) returned status %d, not %d as expected.\n", status, eInvalidArgument);
    }
  }
  /* Null argv. */
  { int argc = 1;
    char ** argv = 0;
    status = ape_trad_init(argc, argv);
    if (eNullPointer != status) {
      ape_test_failed("ape_trad_init(1, 0) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  /* Null argv[0]. */
  { int argc = 1;
    char * argv[] = { 0 };
    status = ape_trad_init(argc, argv);
    if (eNullPointer != status) {
      ape_test_failed("ape_trad_init(1, argv == { 0 }) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  /* Empty argv[0]. */
  { int argc = 1;
    char * argv[] = { "" };
    status = ape_trad_init(argc, argv);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_trad_init(1, argv == \"\") returned status %d, not %d as expected.\n", status, eInvalidArgument);
    }
  }

  /* Non-existent parameter file. */
  { int argc = 1;
    char * argv[] = { "non-existent-task" };
    status = ape_trad_init(argc, argv);
    if (eFileReadError != status) {
      ape_test_failed("ape_trad_init(1, argv == \"%s\") returned status %d, not %d as expected.\n",
        argv[0], status, eInvalidArgument);
    }
  }

  /* Existing parameter file. */
  { char * argv[] = {
/* Note this #ifdef appears backwards, but this is deliberate, to cross-up the parameter files to ensure
   that parameter files written on Windows can be read on Unix and vice versa. */
#ifdef WIN32
      "ape_test_unix",
#else
      "ape_test_win32",
#endif
      "sa",
      "sal",
      "sq",
      "sql",
      "bh",
      "=",
      "yes",
      "fh=.",
      0
    };
    const char * name[] = { "sa", "sal", "sq", "sql", "bh", "fh" };
    const char * expected[] = { "sa", "sal", "sq", "sql", "yes", "." };
    int argc = sizeof(argv) / sizeof(char *) - 1;
    int idx = 0;
    char * pfiles_orig = 0;
    status = ape_io_get_pfiles(&pfiles_orig);
    if (eOK == status) {
      status = ape_io_set_pfiles(QS".");
    }
    if (eOK == status) {
      status = ape_trad_init(argc, argv);
      if (eOK == status) {
        ApeParFile * current = 0;
        status = ape_session_get_current(s_default_session, &current);
        if (eOK == status) {
          for (idx = 0; idx != sizeof(name) / sizeof(char *); ++idx) {
            ApeListIterator itor = 0;
            status = ape_io_find_par(name[idx], current, &itor);
            if (eOK == status) {
              ApePar * par = (ApePar *) ape_list_get(itor);
              char * value = 0;
              char buf[128];
              status = ape_par_get_field(par, eValue, &value);
              sprintf(buf, "after ape_trad_init(8 args) parameter %s", name[idx]);
              ape_test_cmp_string(buf, value, expected[idx], status, eOK);
              free(value); value = 0;
            }
          }
        } else {
          ape_test_failed("after ape_trad_init(8 args) unable to verify parameter values.");
        }
        /* Tests for ape_trad_query_bool. */
        { char value = '\0';
          status = ape_trad_query_bool("bh", &value);
          ape_test_cmp_long("ape_trad_query_bool(\"bh\")", value, '\1', status, eOK);
        }
        { char value = '\1';
          status = ape_trad_query_bool("nobool", &value);
          ape_test_cmp_long("ape_trad_query_bool(\"nobool\")", value, '\0', status, eParNotFound);
        }
        { char value = '\1';
          status = ape_trad_query_bool("sa", &value);
          ape_test_cmp_long("ape_trad_query_bool(\"sa\")", value, '\0', status, eConversionError);
        }
        /* Tests for ape_trad_query_double and ape_trad_set_double. */
        { double value = 1.;
          char msg[128] = "ape_trad_query_double(\"dh\", &value)";
          status = ape_trad_query_double("dh", &value);
          ape_test_cmp_double(msg, value, 0., status, eUndefinedValue);

          value = 1.23456789012345e200;
          sprintf(msg, "ape_trad_set_double(\"dh\", %1.15g)", value);
          status = ape_trad_set_double("dh", value);
          if (eOK == status) {
            char * string_value = 0;
            char expected_string_value[64] = "";
            double expected_value = value;
            sprintf(expected_string_value, "%1.15g", expected_value);
            status = ape_trad_get_string("dh", &string_value);
            ape_test_cmp_string("ape_trad_get_string(\"dh\")", string_value, expected_string_value, status, eOK);
            free(string_value); string_value = 0;
            value = 1.;
            status = ape_trad_get_double("dh", &value);
            ape_test_cmp_double("ape_trad_get_double(\"dh\")", value, expected_value, status, eOK);
          } else {
            ape_test_failed("%s returned status %d, not %d, as expected.\n", msg, status, eOK);
          }
        }

        /* Test error case with ape_trad_query_double. */
        { double value = 0.;
          status = ape_trad_query_double("rhinvalid", &value);
          ape_test_cmp_double("ape_trad_query_double(\"rhinvalid\", ...)", value, 5., status, eStringRemainder);
        }

        /* Tests for ape_trad_query_int and ape_trad_set_int. */
        { int value = 0;
          status = ape_trad_query_int("ih", &value);
          ape_test_cmp_long("ape_trad_query_int(\"ih\")", value, 0, status, eUndefinedValue);
        }
        { int value = 10000;
          status = ape_trad_set_int("ih", value);
          ape_test_cmp_long("ape_trad_set_int(\"ih\", 10000)", value, value, status, eOK);
        }

        /* Tests for ape_trad_query_long and ape_trad_set_long. */
        { long value = 0;
          status = ape_trad_query_long("lh", &value);
          ape_test_cmp_long("ape_trad_query_long(\"lh\")", value, 0, status, eUndefinedValue);
        }
        { long value = 1000000001l;
          status = ape_trad_set_long("lh", value);
          ape_test_cmp_long("ape_trad_set_long(\"lh\", 1000000001l)", value, value, status, eOK);
        }
        if (eOK == status) {
          char * value = 0;
          status = ape_trad_get_string("lh", &value);
          ape_test_cmp_string("ape_trad_get_string(\"lh\")", value, "1000000001", status, eOK);
          free(value); value = 0;
        }
        /* Test ape_trad_get_par_names. */
        { const char * expected[] = {
            "sa", "sal", "sh", "shl", "sq",
            "sql", "sh_comment", "bh", "dh", "fh",
            "frh", "fwh", "ih", "lh", "rh", "bhmin",
            "dhmin", "fhmin", "frhmin", "fwhmin", "ihmin",
            "rhmin", "shmin", "bhenum", "dhenum", "fhenum",
            "frhenum", "fwhenum", "ihenum", "rhenum", "shenum",
            "bhmax", "dhmax", "fhmax", "frhmax", "fwhmax",
            "ihmax", "rhmax", "shmax", "bhrange", "dhrange",
            "fhrange", "frhrange", "fwhrange", "ihrange", "rhrange",
            "shrange", "bhvalid", "dhvalid", "fhvalid", "frhvalid",
            "fwhvalid", "ihvalid", "rhvalid", "shvalid", "bhlow",
            "dhlow", "fhlow", "frhlow", "fwhlow", "ihlow",
            "rhlow", "shlow", "bhhigh", "dhhigh", "fhhigh",
            "frhhigh", "fwhhigh", "ihhigh", "rhhigh", "shhigh",
            "dhinvalid", "ihinvalid", "rhinvalid", "ihindef", "ihinf",
            "ihinfinity", "ihnan", "ihnone", "ihundef", "ihundefined",
            "rhindef", "rhinf", "rhinfinity", "rhnan", "rhnone",
            "rhundef", "rhundefined", "shspace", "mode",
            0
          };
          char ** par_name = 0;
          status = ape_trad_get_par_names(&par_name);
          ape_test_cmp_string_array("ape_trad_get_par_names", par_name, expected, status, eOK);
          ape_util_free_string_array(par_name); par_name = 0;
        }
        /* Test ape_trad_get_type. */
        { const char * par_name = "frh";
          char par_type[APE_PAR_TYPE_CODE_LEN] = "qZX";
          const char * expected_type = "fr";
          int expected_status = eOK;
          status = ape_trad_get_type(par_name, par_type);
          ape_test_cmp_string("ape_trad_get_type(\"frh\")", par_type, expected_type, status, expected_status);

          par_name = "non-existent-par";
          expected_type = "";
          expected_status = eParNotFound;
          status = ape_trad_get_type(par_name, par_type);
          ape_test_cmp_string("ape_trad_get_type(\"non-existent-par\")", par_type, expected_type, status, expected_status);
        }
      } else {
        ape_test_failed("ape_trad_init(argv == \"%s\", ...) returned status %d, not %d as expected.\n",
          argv[0], status, eOK);
      }
      ape_trad_close(0);
      ape_io_set_pfiles(pfiles_orig);
    } else {
      ape_test_failed("Could not set up to test ape_trad_init (status is %d.)\n", status);
    }
    free(pfiles_orig); pfiles_orig = 0;
  }
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_trad.c,v $
 * Revision 1.50  2011/01/21 21:29:54  irby
 * Add test for parameter whose default value is a single space.  Test by
 * making sure a parameter file open/close with no other action leaves the
 * file completely unchanged even if saved.
 *
 * Revision 1.49  2010/11/23 22:05:56  jpeachey
 * Replace tabs with spaces to make code line up correctly.
 *
 * Revision 1.48  2010/11/12 20:47:45  irby
 * Add ape_trad_set_short() and cfortran macros for ape_trad_set_double,
 * ape_trad_set_long, and ape_trad_set_short.  Also, add string & filename
 * cfortran wrapper functions (including string handling code).
 *
 * Revision 1.47  2010/09/10 19:08:34  irby
 * Add ape_trad_get_byte() wrapper for ape_trad_get_bool() as cfortran has no
 * char, only byte (signed char).  Add ape_trad_get_int() & ape_trad_set_int()
 * (also ape_session_get_int() & ape_session_set_int() and related tests) for
 * Fortran wrappers.
 *
 * Revision 1.46  2010/08/27 20:40:56  irby
 * Add cfortran wrappers for some ape_trad_* functions (incomplete).
 *
 * Revision 1.45  2010/06/02 20:21:33  peachey
 * Replace original code with more lightweight calls to ape_session
 * equivalents, which now replace all the ape_trad functionality.
 *
 * Revision 1.44  2010/06/02 19:45:51  peachey
 * Correct typo in name of function being tested.
 *
 * Revision 1.43  2010/06/02 19:07:55  peachey
 * Remove unneeded enum and add explicit initialization for default session.
 *
 * Revision 1.42  2010/06/02 18:54:27  peachey
 * Use internal ApeSession structure to pave the way for multi-session support.
 *
 * Revision 1.41  2007/11/16 18:21:08  peachey
 * Don't call ape_io_check_file_format from ape_trad_init.
 *
 * Revision 1.40  2007/11/12 19:51:33  peachey
 * Call ape_util_interpret_env at beginning of ape_trad_init.
 *
 * Revision 1.39  2007/11/12 19:22:14  peachey
 * Remove extraneous calls to ape_par_check in the ape_trad_get* functions,
 * and test getting a double-valued parameter that has an error to make
 * sure the message appears correctly.
 *
 * Revision 1.38  2007/11/12 16:54:52  peachey
 * Check status more specifically when calling ape_io_get_default_mode.
 * Ignore error only if the mode parameter was not found.
 *
 * Revision 1.37  2007/05/15 20:52:48  peachey
 * Additional speed optimizations:
 * o Reduce number of calls to ape_list_end by calling it once and storing
 *   the returned iterator for use in loop logic when it is safe to do so.
 * o Make ape_list_destroy more efficient by iterating through the list
 *   back to front once and freeing the elements instead of repeatedly
 *   calling ape_list_remove_entry.
 *
 * Revision 1.36  2007/01/25 17:49:57  peachey
 * Fixed bug in ape_par_get_string_case that caused it not to behave
 * correctly when called for a non-enumerated parameter. Added functions
 * ape_trad_prompt_string_case and ape_trad_get_string_case for getting
 * parameters with control over case of returned string.
 *
 * Revision 1.35  2006/11/30 20:38:35  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.34  2006/11/30 16:40:27  peachey
 * Fix bug in test code: pfiles needs to be set correctly for Windows test.
 *
 * Revision 1.33  2006/11/30 16:28:47  peachey
 * Fix bug in test code: pfiles needs to be set correctly for Windows test.
 *
 * Revision 1.32  2006/11/24 19:47:28  peachey
 * Add ape_trad_query_float, ape_trad_get_float and ape_trad_set_float.
 *
 * Revision 1.31  2006/08/25 19:45:58  peachey
 * Add tests for undefined and infinite numerical parameters.
 *
 * Revision 1.30  2006/06/06 13:30:21  peachey
 * Add ape_trad_find_par, for getting underlying ape parameter object.
 *
 * Revision 1.29  2006/06/05 18:59:34  peachey
 * Add apd_trad_get_type for getting the type code of a parameter.
 *
 * Revision 1.28  2006/05/31 03:53:00  peachey
 * Add body of ape_trad_query_file_name.
 *
 * Revision 1.27  2006/05/31 03:21:44  peachey
 * Rationalize ape_par_get_* family of functions, and use them to implement
 * ape_trad_get_* and ape_trad_query_*.
 *
 * Revision 1.26  2006/05/31 01:36:51  peachey
 * Rename ape_par_get_string to ape_par_get_field and ape_par_get_string_array to
 * ape_par_get_field_array.
 *
 * Revision 1.25  2006/05/26 15:22:11  peachey
 * Add ape_trad_set* functions complete ape_trad_query* functions, for completeness
 * and to support PIL compatibility.
 *
 * Revision 1.24  2006/05/23 19:31:40  peachey
 * Add ape_trad_get_sys_pfile, for getting the system parameter file.
 *
 * Revision 1.23  2006/05/23 16:27:25  peachey
 * Update call to ape_io_write to add new force_write argument.
 *
 * Revision 1.22  2006/05/22 17:35:40  peachey
 * Add ape_trad_query_string, for prompting for and getting string parameters.
 *
 * Revision 1.21  2006/05/18 14:00:50  peachey
 * Rename ape_trad_get_bool ape_trad_query_bool. Add ape_trad_get_string
 * which only gets string, without prompting.
 *
 * Revision 1.20  2006/05/18 03:18:54  peachey
 * Correct handling of situation in which local pfiles is not empty
 * but no file exists in that location yet.
 *
 * Revision 1.19  2006/05/17 03:15:30  peachey
 * Add ape_trad_get_current, for getting the current parameter file object.
 *
 * Revision 1.18  2006/05/17 02:20:24  peachey
 * Add ape_trad_get_par_names, for getting all named parameters in a file.
 *
 * Revision 1.17  2006/05/16 19:47:22  peachey
 * Set pfiles to ;. to prevent local file from being modified before
 * testing ape_trad_init. Consolidate ape_trad_init tests.
 *
 * Revision 1.16  2006/05/16 17:22:27  peachey
 * Update ape_trad tests to use latest, greatest ape_test.par.
 *
 * Revision 1.15  2006/05/16 15:13:21  peachey
 * Add note reminding of need to write better test for ape_trad_init.
 *
 * Revision 1.14  2006/05/16 15:09:46  peachey
 * Add test for ape_trad_init in which the long and complex ape_test.par is opened.
 *
 * Revision 1.13  2006/05/13 18:45:56  peachey
 * Add ape_trad_save and use it to save parameter file after each prompt.
 *
 * Revision 1.12  2006/05/12 17:33:21  peachey
 * Add fourth default parameter file and apply command line arguments to that;
 * preserve pre-command line version as 'merged' version.
 *
 * Revision 1.11  2006/05/12 03:32:09  peachey
 * Add ape_trad_get_bool, and infrastructure for all the ape_trad_get_... family.
 *
 * Revision 1.10  2006/05/10 16:48:50  peachey
 * Generalize format tests to allow them to be run with or without
 * value checking, to let them be used both in contexts where the whole parameter
 * needs to be checked (checking user input) and contexts in which a bad value
 * is not a show-stopper, because the user or client may yet correct the problem.
 *
 * Revision 1.9  2006/05/09 19:32:08  peachey
 * In ape_trad_init, check the merged parameter file using ape_io_check_file_format.
 *
 * Revision 1.8  2006/05/06 00:43:16  peachey
 * Add most of implemention and tests for ape_trad_init, which is
 * modeled after PILInit.
 *
 * Revision 1.7  2006/04/28 01:24:04  peachey
 * Add a couple casts to get the constness right when calling free.
 * Enable al mode, and remove eFileOpenError, using eFileReadError in lieu thereof.
 * The latter change was to make error messages consistent between Windows and
 * Unix, because although the two both fail to open a directory as if it were
 * a file, they fail at different points in the process of opening the file.
 *
 * Revision 1.6  2006/04/19 16:04:20  peachey
 * Add rational shutdown/cleanup/atexit stuff. Public function
 * ape_trad_close performs all shutdown/cleanup. Static function ape_trad_destroy
 * does all deallocations. Static function ape_trad_atexit is registered with
 * ape_util_atexit, and calls ape_trad_close. System and local par files are now
 * kept in static storage by ape_trad_init.
 *
 * Revision 1.5  2006/04/19 03:07:09  peachey
 * Uncomment final effort to open local parameter file, and correct
 * the expected status.
 *
 * Revision 1.4  2006/04/15 01:56:29  peachey
 * Revamp ape_trad_init using recent developments in ape_io.c to get local
 * and system paths, and from them to load local and system parameter
 * files.
 *
 * Revision 1.3  2006/04/13 18:44:44  peachey
 * Continue to chip away at ape_trad_init, adding get_pfile_name helper function.
 *
 * Revision 1.2  2006/04/12 14:29:29  peachey
 * Add forgotten header file.
 *
 * Revision 1.1  2006/04/12 14:18:14  peachey
 * Add beginning of traditional facility for duplicating PIL/XPI behavior.
 *
*/
