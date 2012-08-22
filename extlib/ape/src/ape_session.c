/** \file ape_session.c
    \brief Implementation of traditional XPI/PIL-compliant interface.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_io.h"
#include "ape/ape_list.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_test.h"
#include "ape/ape_session.h"
#include "ape/ape_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

struct ApeSession {
  ApeParFile * system;
  ApeParFile * local;
  ApeParFile * merged;
  ApeParFile * current;
};

static void ape_session_destroy(ApeSession * session) {
  if (0 != session) {
    ape_io_destroy_file(session->current); session->current = 0;
    ape_io_destroy_file(session->merged); session->merged = 0;
    ape_io_destroy_file(session->local); session->local = 0;
    ape_io_destroy_file(session->system); session->system = 0;
    free(session);
  }
}

static int find_par(ApeSession * session, const char * par_name, ApePar ** par) {
  int status = eOK;
  ApeParFile * par_file = 0;
  ApeListIterator itor = 0;

  if (0 != par) *par = 0;
  else status = eNullPointer;

  if (eOK == status && 0 == session) status = eNullPointer;
  if (eOK == status && 0 == par_name) status = eNullPointer;

  if (eOK == status) {
    par_file = session->current;
  }

  if (eOK == status && 0 != par_file) {
    /* Get the named parameter. */
    status = ape_io_find_par(par_name, par_file, &itor);
  } else {
    status = eNoParLoaded;
  }

  if (eOK == status) {
    *par = (ApePar *) ape_list_get(itor);
  }

  return status;
}

/* Internal utility which takes the full path to the executable as an argument, and looks for the
   last directory separator (think slash). Everything after the last separator is considered the
   name of the tool. */
static int get_pfile_name(const char * full_exec_name, char ** pfile_name) {
  size_t num_sep = sizeof(s_dir_separator) / sizeof(const char);
  size_t idx = 0;
  size_t last_idx = 0;
  const char * last_sep = 0;
  char * exec_name = 0;
  int status = eOK;

  /* Iterate over valid separators. */
  for (idx = 0; idx != num_sep; ++idx) {
    /* Look for right-most separator. */
    const char * sep = strrchr(full_exec_name, s_dir_separator[idx]);

    /* In case multiple separators are found, pick the right-most one. */
    if (sep > last_sep) {
      last_idx = idx;
      last_sep = sep;
    }
  }

  /* If a separator was found, start after it, otherwise just use the whole file name. */
  if (0 != last_sep) last_sep += sizeof(s_dir_separator[last_idx]) / sizeof(const char);
  else last_sep = full_exec_name;

  /* Make a copy of the file name part of the full name. */
  status = ape_util_copy_string(last_sep, &exec_name);
  if (0 == status) {
    /* Terminate the copy at the last period, in order to truncate any extension (e.g. .exe) before adding .par. */
    char * extension = strrchr(exec_name, '.');
    if (0 != extension) *extension = '\0';

    /* Add the .par in lieu of the extension. */
    status = ape_util_cat_string(exec_name, ".par", pfile_name);
  }
  /* Clean up. */
  free(exec_name);
  return status;
}

int ape_session_init(ApeSession ** session, int argc, char ** argv) {
  int status = eOK;
  char * pfile_name = 0;
  ApeSession * new_session = 0;

  /* Perform initialization of ape based on environment variables etc. */
  status = ape_util_interpret_env();

  if (eOK == status) {
    if (0 != session) {
      /* Destroy whatever parameter objects are currently open, saving them first and ignoring status. */
      ape_session_close(*session, 1);
      *session = 0;
    }

    /* Check command line arguments for problems. */
    if (0 == session) {
      status = eNullPointer;
      ape_msg_debug("ape_session_init was called with session == 0.\n");
    }
    if (0 == argv) {
      status = eNullPointer;
      ape_msg_debug("ape_session_init was called with argv == 0.\n");
    }
    if (eOK == status) {
      if (0 == argv[0]) {
        status = eNullPointer;
        ape_msg_debug("ape_session_init was called with argv[0] == 0.\n");
      } else if (0 == *argv[0]) {
        status = eInvalidArgument;
        ape_msg_debug("ape_session_init was called with argv[0] == \"\".\n");
      }
    }
    if (0 >= argc) {
      status = eInvalidArgument;
      ape_msg_debug("ape_session_init was called with argc <= 0.\n");
    }
  }

  if (eOK == status) {
    /* Create a new session object. */
    new_session = (ApeSession *) calloc(1, sizeof(ApeSession));
    if (0 == new_session) {
      status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* Get the name of the parameter file, based on this executable name. */
    status = get_pfile_name(argv[0], &pfile_name);
  }

  if (eOK == status) {
    ApeList * loc_path = 0;
    ApeList * sys_path = 0;

    /* Get the system par file path, if any. */
    status = ape_io_get_sys_path(&sys_path);
    if (eOK == status) {
      /* Ignore the status because it's OK if this fails. */
      ape_io_read_file_path(pfile_name, sys_path, &new_session->system);
    }

    /* Get the local par file path, if any. */
    status = ape_io_get_loc_path(&loc_path);
    if (eOK == status) {
      /* Try to open a parameter file on the local path. */
      status = ape_io_read_file_path(pfile_name, loc_path, &new_session->local);
      if (eOK != status) {
        ApeListIterator itor = ape_list_begin(loc_path);
        if (itor != ape_list_end(loc_path)) {
          char * full_file_name = 0;

          /* Local path contains at least one directory, so clone the system file. */
          status = ape_io_clone_file(new_session->system, &new_session->local);

          if (eOK == status) {
            const char * dir_name = (const char *) ape_list_get(itor);
            /* Get the name of the local parameter file. */
            status = ape_util_append_file_name(dir_name, pfile_name, &full_file_name);
          }

          if (eOK == status) {
            /* Ignore the status because it's OK if this fails. */
            ape_io_set_file_name(new_session->local, full_file_name);
          }
          free(full_file_name); full_file_name = 0;
        }
      }
    }
    /* Some amount of error above may be tolerated. If something unrecoverable happened it will be detected below. */
    status = eOK;

    if (0 == new_session->local && 0 == new_session->system) {
      /* No luck searching the paths: try the parameter file name by itself. */
      status = ape_io_read(pfile_name, &new_session->local);
    }
  }

  if (eOK == status) {
    /* Merge the system and local parameter files to create a file which has the best set of parameter values. */
    status = ape_io_merge_par_files(new_session->system, new_session->local, &new_session->merged);
  }

  if (eOK == status) {
    /* Clone the merged parameter file to create the "current" parameter file. */
    status = ape_io_clone_file(new_session->merged, &new_session->current);
  }

  if (eOK == status) {
    /* Modify the merged parameter file using the command line arguments. */
    status = ape_io_apply_command_line(new_session->current, argc - 1, argv + 1);
  }

  if (eOK == status) {
    *session = new_session;
  } else {
    ape_session_destroy(new_session);
  }

  /* Clean up. */
  free(pfile_name); pfile_name = 0;

  return status;
}

int ape_session_close(ApeSession * session, int save_flag) {
  if (0 != save_flag) ape_session_save(session);
  ape_session_destroy(session);
  return eOK;
}

int ape_session_get_current(ApeSession * session, ApeParFile ** par_file) {
  int status = eOK;

  if (0 == session) {
    status = eNullPointer;
  } else if (0 != par_file) {
    *par_file = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    if (0 != session->current) {
      *par_file = session->current;
    } else {
      status = eUninitialized;
    }
  }

  return status;
}

int ape_session_get_sys_pfile(ApeSession * session, ApeParFile ** par_file) {
  int status = eOK;

  if (0 == session) {
    status = eNullPointer;
  } else if (0 != par_file) {
    *par_file = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    if (0 != session->system) {
      *par_file = session->system;
    } else {
      status = eUninitialized;
    }
  }

  return status;
}

int ape_session_get_par_names(ApeSession * session, char *** par_name) {
  int status = eOK;
  ApeList * par_cont = 0;
  ApeListIterator itor = 0;
  ApeListIterator end = 0;
  size_t num_par = 0;
  size_t idx = 0;

  if (0 == session) {
    status = eNullPointer;
  } else if (0 != par_name) {
    *par_name = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == session->current) {
    status = eUninitialized;
  }

  if (eOK == status) {
    status = ape_io_get_par_cont(session->current, &par_cont);
  }

  if (eOK == status) {
    num_par = ape_list_get_size(par_cont);
  }

  if (eOK == status) {
    *par_name = (char **) calloc(num_par + 1, sizeof(char *));
  }

  if (eOK == status) {
    end = ape_list_end(par_cont);
  }

  for (itor = ape_list_begin(par_cont); itor != end; itor = ape_list_next(itor)) {
    char * name = 0;
    ApePar * par = (ApePar *) ape_list_get(itor);
    int local_status = ape_par_get_field(par, eName, &name);
    if (eOK == local_status) {
      if (0 == name || '\0' == *name) {
        free(name); name = 0;
      } else {
        (*par_name)[idx] = name;
        ++idx;
      }
    } else {
      free(name); name = 0;
      status = eOK == status ? local_status : status;
    }
  }

  return status;
}

int ape_session_get_type(ApeSession * session, const char * par_name, char * par_type) {
  int status = eOK;
  ApePar * par = 0;

  if (0 != par_type) *par_type = 0;
  else status = eNullPointer;

  if (eOK == status) status = find_par(session, par_name, &par);

  if (eOK == status) status = ape_par_get_type(par, par_type);

  return status;
}

static int query_value(ApeSession * session, const char * par_name, ApePar ** par) {
  int status = eOK;
  ApeListIterator itor = 0;
  ApeParFile * par_file = 0;
  char default_mode[APE_PAR_MODE_CODE_LEN] = "";

  /* Check arguments. */
  if (0 != par) {
    *par = '\0';
  } else if (eOK == status) {
    status = eNullPointer;
  }
  if (eOK == status && 0 == session) status = eNullPointer;
  if (eOK == status && 0 == par_name) status = eNullPointer;

  par_file = session->current;

  if (eOK == status && 0 != par_file) {
    /* Find the named parameter in the merged (working) parameter file. */
    status = ape_io_find_par(par_name, par_file, &itor);
  } else {
    status = eNoParLoaded;
  }

  if (eOK == status) {
    /* Get the default mode from the parameter file. */
    status = ape_io_get_default_mode(par_file, default_mode);

    if (eParNotFound == status) {
      /* No default mode defined, so use ql. */
      status = eOK;
      strcpy(default_mode, "ql");
    }
  }

  if (eOK == status) {
    *par = (ApePar *) ape_list_get(itor);

    /* Query for the parameter. */
    status = ape_par_query(*par, default_mode);

    /* Ignore status; try to save the current parameters here. */
    ape_session_save(session);
  }

  return status;
}

int ape_session_query_bool(ApeSession * session, const char * par_name, char * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_bool(par, value);
  }

  return status;
}

int ape_session_query_double(ApeSession * session, const char * par_name, double * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_double(par, value);
  }

  return status;
}

int ape_session_query_float(ApeSession * session, const char * par_name, float * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_float(par, value);
  }

  return status;
}

int ape_session_query_file_name(ApeSession * session, const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_file_name(par, value);
  }

  return status;
}

int ape_session_query_int(ApeSession * session, const char * par_name, int * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_int(par, value);
  }

  return status;
}

int ape_session_query_long(ApeSession * session, const char * par_name, long * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_long(par, value);
  }

  return status;
}

int ape_session_query_string(ApeSession * session, const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_string(par, value);
  }

  return status;
}

int ape_session_query_string_case(ApeSession * session, const char * par_name, char ** value, char case_code) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(session, par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_string_case(par, value, case_code);
  }

  return status;
}

int ape_session_get_bool(ApeSession * session, const char * par_name, char * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_bool(par, value);

  return status;
}

int ape_session_get_double(ApeSession * session, const char * par_name, double * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_double(par, value);

  return status;
}

int ape_session_get_float(ApeSession * session, const char * par_name, float * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_float(par, value);

  return status;
}

int ape_session_get_file_name(ApeSession * session, const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_file_name(par, value);

  return status;
}

int ape_session_get_int(ApeSession * session, const char * par_name, int * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_int(par, value);

  return status;
}

int ape_session_get_long(ApeSession * session, const char * par_name, long * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_long(par, value);

  return status;
}

int ape_session_get_string(ApeSession * session, const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_string(par, value);

  return status;
}

int ape_session_get_string_case(ApeSession * session, const char * par_name, char ** value, char case_code) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(session, par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_string_case(par, value, case_code);

  return status;
}

int ape_session_set_bool(ApeSession * session, const char * par_name, char value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    if (0 != value) status = ape_par_set_value_string(par, "yes");
    else status = ape_par_set_value_string(par, "no");
  }

  return status;
}

int ape_session_set_double(ApeSession * session, const char * par_name, double value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%1.15g", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_session_set_float(ApeSession * session, const char * par_name, float value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%1.7g", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_session_set_file_name(ApeSession * session, const char * par_name, const char * value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    status = ape_par_set_value_string(par, value);
  }

  return status;
}

int ape_session_set_int(ApeSession * session, const char * par_name, int value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%d", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_session_set_long(ApeSession * session, const char * par_name, long value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%ld", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_session_set_short(ApeSession * session, const char * par_name, short value) {
  return ape_session_set_int(session, par_name, value);
}

int ape_session_set_string(ApeSession * session, const char * par_name, const char * value) {
  ApePar * par = 0;
  int status = find_par(session, par_name, &par);

  if (eOK == status) {
    status = ape_par_set_value_string(par, value);
  }

  return status;
}

int ape_session_save(ApeSession * session) {
  int status = eOK;
  ApeParFile * tmp_file = 0;

  if (0 != session) {
    /* Make a temporary copy of current parameters. */
    status = ape_io_clone_file(session->current, &tmp_file);

    if (eOK == status) {
      /* Reset unlearned parameters in the copy to their pre-command-line state. */
      status = ape_io_revert_unlearned(tmp_file, session->merged);
    }

    if (eOK == status) {
      /* Write temporary copy. */
      status = ape_io_write(tmp_file, 0);
    }

    /* Clean up. */
    ape_io_destroy_file(tmp_file); tmp_file = 0;
  }

  return status;
}

int ape_session_find_par(ApeSession * session, const char * par_name, ApePar ** par) {
  return find_par(session, par_name, par);
}

void ape_session_test(void) {
  int status = eOK;

  /* Test get_pfile_name. */
  /* TODO: test more thoroughly. */
  { char * pfile_name = 0;
    const char * expected_pfile_name = "binary.par";
    status = get_pfile_name("c:/some/path/to\\a/binary.exe", &pfile_name);
    if (0 == pfile_name) {
      ape_test_failed("ape_session_test: get_pfile_name(long path) returned pfile_name == 0, not \"%s\" as expected.\n",
        expected_pfile_name);
    } else if (0 != strcmp(pfile_name, expected_pfile_name)) {
      ape_test_failed("ape_session_test: get_pfile_name(long path) returned \"%s\", not \"%s\" as expected.\n",
        pfile_name, expected_pfile_name);
    }
    free(pfile_name); pfile_name = 0;
  }
  { char * pfile_name = 0;
    const char * expected_pfile_name = "binary.par";
    status = get_pfile_name("binary.exe", &pfile_name);
    if (0 == pfile_name) {
      ape_test_failed("ape_session_test: get_pfile_name(short path) returned pfile_name == 0, not \"%s\" as expected.\n",
        expected_pfile_name);
    } else if (0 != strcmp(pfile_name, expected_pfile_name)) {
      ape_test_failed("ape_session_test: get_pfile_name(short path) returned \"%s\", not \"%s\" as expected.\n",
        pfile_name, expected_pfile_name);
    }
    free(pfile_name); pfile_name = 0;
  }

  /* Attempt to initialize a task using invalid inputs. */
  /* Null argc and argv. */
  { int argc = 0;
    char ** argv = 0;
    ApeSession * session = 0;
    status = ape_session_init(&session, argc, argv);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_session_init(0, 0, 0) returned status %d, not %d as expected.\n", status, eInvalidArgument);
    }
  }
  /* Null argv. */
  { int argc = 1;
    char ** argv = 0;
    ApeSession ** session = 0;
    status = ape_session_init(session, argc, argv);
    if (eNullPointer != status) {
      ape_test_failed("ape_session_init(0, 1, 0) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  /* Null argv[0]. */
  { int argc = 1;
    char * argv[] = { 0 };
    ApeSession * session = 0;
    status = ape_session_init(&session, argc, argv);
    if (eNullPointer != status) {
      ape_test_failed("ape_session_init(&session, 1, argv == { 0 }) returned status %d, not %d as expected.\n", status,
        eNullPointer);
    }
  }
  /* Empty argv[0]. */
  { int argc = 1;
    char * argv[] = { "" };
    ApeSession * session = 0;
    status = ape_session_init(&session, argc, argv);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_session_init(&session1, argv == \"\") returned status %d, not %d as expected.\n", status,
        eInvalidArgument);
    }
  }

  /* Non-existent parameter file. */
  { int argc = 1;
    char * argv[] = { "non-existent-task" };
    ApeSession * session = 0;
    status = ape_session_init(&session, argc, argv);
    if (eFileReadError != status) {
      ape_test_failed("ape_session_init(&session, 1, argv == \"%s\") returned status %d, not %d as expected.\n",
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
      ApeSession * session = 0;
      status = ape_session_init(&session, argc, argv);
      if (eOK == status) {
        for (idx = 0; idx != sizeof(name) / sizeof(char *); ++idx) {
          ApeListIterator itor = 0;
          status = ape_io_find_par(name[idx], session->current, &itor);
          if (eOK == status) {
            ApePar * par = (ApePar *) ape_list_get(itor);
            char * value = 0;
            char buf[128];
            status = ape_par_get_field(par, eValue, &value);
            sprintf(buf, "after ape_session_init(8 args) parameter %s", name[idx]);
            ape_test_cmp_string(buf, value, expected[idx], status, eOK);
            free(value); value = 0;
          }
        }
        /* Tests for ape_session_query_bool. */
        { char value = '\0';
          status = ape_session_query_bool(session, "bh", &value);
          ape_test_cmp_long("ape_session_query_bool(session, \"bh\")", value, '\1', status, eOK);
        }
        { char value = '\1';
          status = ape_session_query_bool(session, "nobool", &value);
          ape_test_cmp_long("ape_session_query_bool(session, \"nobool\")", value, '\0', status, eParNotFound);
        }
        { char value = '\1';
          status = ape_session_query_bool(session, "sa", &value);
          ape_test_cmp_long("ape_session_query_bool(session, \"sa\")", value, '\0', status, eConversionError);
        }
        /* Tests for ape_session_query_double and ape_session_set_double. */
        { double value = 1.;
          char msg[128] = "ape_session_query_double(session, \"dh\", &value)";
          status = ape_session_query_double(session, "dh", &value);
          ape_test_cmp_double(msg, value, 0., status, eUndefinedValue);

          value = 1.23456789012345e200;
          sprintf(msg, "ape_session_set_double(session, \"dh\", %1.15g)", value);
          status = ape_session_set_double(session, "dh", value);
          if (eOK == status) {
            char * string_value = 0;
            char expected_string_value[64] = "";
            double expected_value = value;
            sprintf(expected_string_value, "%1.15g", expected_value);
            status = ape_session_get_string(session, "dh", &string_value);
            ape_test_cmp_string("ape_session_get_string(session, \"dh\")", string_value, expected_string_value, status, eOK);
            free(string_value); string_value = 0;
            value = 1.;
            status = ape_session_get_double(session, "dh", &value);
            ape_test_cmp_double("ape_session_get_double(session, \"dh\")", value, expected_value, status, eOK);
          } else {
            ape_test_failed("%s returned status %d, not %d, as expected.\n", msg, status, eOK);
          }
        }

        /* Test error case with ape_session_query_double. */
        { double value = 0.;
          status = ape_session_query_double(session, "rhinvalid", &value);
          ape_test_cmp_double("ape_session_query_double(session, \"rhinvalid\", ...)", value, 5., status, eStringRemainder);
        }

        /* Tests for ape_session_query_int and ape_session_set_int. */
        { int value = 0;
          status = ape_session_query_int(session, "ih", &value);
          ape_test_cmp_long("ape_session_query_int(session, \"ih\")", value, 0, status, eUndefinedValue);
        }
        { int value = 10000;
          status = ape_session_set_int(session, "ih", value);
          ape_test_cmp_long("ape_session_set_int(session, \"ih\", 10000)", value, value, status, eOK);
        }

        /* Tests for ape_session_query_long and ape_session_set_long. */
        { long value = 0;
          status = ape_session_query_long(session, "lh", &value);
          ape_test_cmp_long("ape_session_query_long(session, \"lh\")", value, 0, status, eUndefinedValue);
        }
        { long value = 1000000001l;
          status = ape_session_set_long(session, "lh", value);
          ape_test_cmp_long("ape_session_set_long(session, \"lh\", 1000000001l)", value, value, status, eOK);
        }
        if (eOK == status) {
          char * value = 0;
          status = ape_session_get_string(session, "lh", &value);
          ape_test_cmp_string("ape_session_get_string(session, \"lh\")", value, "1000000001", status, eOK);
          free(value); value = 0;
        }
        /* Test ape_session_get_par_names. */
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
          status = ape_session_get_par_names(session, &par_name);
          ape_test_cmp_string_array("ape_session_get_par_names", par_name, expected, status, eOK);
          ape_util_free_string_array(par_name); par_name = 0;
        }
        /* Test ape_session_get_type. */
        { const char * par_name = "frh";
          char par_type[APE_PAR_TYPE_CODE_LEN] = "qZX";
          const char * expected_type = "fr";
          int expected_status = eOK;
          status = ape_session_get_type(session, par_name, par_type);
          ape_test_cmp_string("ape_session_get_type(session, \"frh\")", par_type, expected_type, status, expected_status);

          par_name = "non-existent-par";
          expected_type = "";
          expected_status = eParNotFound;
          status = ape_session_get_type(session, par_name, par_type);
          ape_test_cmp_string("ape_session_get_type(session, \"non-existent-par\")", par_type, expected_type, status,
            expected_status);
        }
      } else {
        ape_test_failed("ape_session_init(&session, argv == \"%s\", ...) returned status %d, not %d as expected.\n",
          argv[0], status, eOK);
      }
      ape_session_close(session, 0);
      ape_io_set_pfiles(pfiles_orig);
    } else {
      ape_test_failed("Could not set up to test ape_session_init (status is %d.)\n", status);
    }
    free(pfiles_orig); pfiles_orig = 0;
  }
  /* Test that an open/close with no other action leaves the parameter
   * file completely unchanged even if saved: */
  { char * argv[] = {
   /* For this test match OS since we're saving the file. */
#ifdef WIN32
      "ape_test_win32",
#else
      "ape_test_unix",
#endif
      0
    };
    int argc = sizeof(argv) / sizeof(char *) - 1;
    char * pfiles_orig = 0;
    status = ape_io_get_pfiles(&pfiles_orig);
    if (eOK == status) {
      status = ape_io_set_pfiles("."QS);
    }
    if (eOK == status) {
      ApeSession * session = 0;
      status = ape_session_init(&session, argc, argv);
      if (eOK != status) {
        ape_test_failed("ape_session_init(&session, argv == \"%s\", ...) returned status %d, not %d as expected.\n",
          argv[0], status, eOK);
      }
      ape_session_close(session, 1);
      ape_io_set_pfiles(pfiles_orig);
    } else {
      ape_test_failed("Could not set up to test ape_session_init (status is %d.)\n", status);
    }
    free(pfiles_orig); pfiles_orig = 0;
  }
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_session.c,v $
 * Revision 1.8  2011/01/21 21:29:54  irby
 * Add test for parameter whose default value is a single space.  Test by
 * making sure a parameter file open/close with no other action leaves the
 * file completely unchanged even if saved.
 *
 * Revision 1.7  2010/11/12 20:49:15  irby
 * Add ape_session_set_short() (wraps to ape_session_set_int).
 *
 * Revision 1.6  2010/09/10 19:08:34  irby
 * Add ape_trad_get_byte() wrapper for ape_trad_get_bool() as cfortran has no
 * char, only byte (signed char).  Add ape_trad_get_int() & ape_trad_set_int()
 * (also ape_session_get_int() & ape_session_set_int() and related tests) for
 * Fortran wrappers.
 *
 * Revision 1.5  2010/06/02 20:19:34  peachey
 * Remove spurious functions ape_session_get_value_string and ape_session_get_value_string_case.
 * Add correctly named ape_session_get_string_case in the correct spot.
 *
 * Revision 1.4  2010/06/02 20:03:57  peachey
 * Remove/replace usages of static global ApeSession. Instead,
 * use ape_session_init and ape_session_close to create sessions dynamically.
 *
 * Revision 1.3  2010/06/02 19:57:05  peachey
 * Convert remaining internals to use the new structures.
 *
 * Revision 1.2  2010/06/02 19:47:18  peachey
 * Add ApeSession argument as first parameter to all ape_session calls.
 * Use local pointer to ApeSession instead of the global (static) pointer
 * whenever possible.
 *
 * Revision 1.1  2010/06/02 19:22:10  peachey
 * Add initial (partial) implementation of ape_session module.
 *
 */
