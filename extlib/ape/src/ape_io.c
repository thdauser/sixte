/** \file ape_io.c
    \brief Implementation of input/output facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_io.h"
#include "ape/ape_error.h"
#include "ape/ape_list.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_test.h"
#include "ape/ape_util.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* OS-specific headers are for getting host name and process id, used to construct temporary file names. */
#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>

/* Macros with short and likely unused names, representing:
   o QC (mnemonic is "quote colon"): the token dividing directories in a path.
   o QD (mnemonic is "quote path"): the divider between sub-directories in a directory.
   o QS (mnemonic is "quote semicolon"): the divider between local and system parts of PFILES variable.
*/
/* TODO: Allow / and \ for windows. */
#ifdef WIN32
#define QC ";"
#define QS "|"
#define QD "\\"
#else
#define QC ":"
#define QS ";"
#define QD "/"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Structure representing a parameter file, with associated loaded parameter container structure.
*/
struct ApeParFile {
  char * file_name;
  ApeList * par_cont;
  char read_only;
  time_t mod_time;
};

static struct {
  char * pfiles;
  ApeList * loc_path;
  ApeList * sys_path;
  int init;
} s_par_path = { 0, 0, 0, 0 };

static void ape_io_atexit(void);
static void destroy_paths(void);
static int parse_path(ApeList * list, char * begin);
static void test_apply_command_line(void);
static void test_one_command_line(int argc, char ** argv, int num_par, const char ** par_name, const char ** par_value,
  int expected_status);
static void cat_cmd_line(int argc, char ** argv, char * line);

static int create_paths(const char * pfiles) {
  int status = eOK;

  /* Call clean-up code upon exit. Ignore errors if this fails. */
  ape_util_atexit(&ape_io_atexit);

  /* First, handle the pfiles string itself. */
  if (0 != pfiles) {
    /* Copy the value of pfiles supplied by the client. */
    status = ape_util_copy_string(pfiles, &s_par_path.pfiles);
  } else {
    /* Client supplied pfiles == 0 => get pfiles from environment. */
    status = ape_util_getenv("PFILES", &s_par_path.pfiles, "");

    /* It's OK in this case if PFILES was not set. */
    if (eVarNotSet == status) status = eOK;
  }

  /* If all goes well, create lists to hold the local path and the system path. */
  if (eOK == status) {
    s_par_path.loc_path = ape_list_create();
    if (0 != s_par_path.loc_path) {
      s_par_path.sys_path = ape_list_create();
      if (0 != s_par_path.sys_path) {
        /* Find division between local and system paths. */
        char * begin = s_par_path.pfiles;
        char * end = strstr(s_par_path.pfiles, QS);
        char * sys = 0;

        if (0 != end) {
          /* QS was found, so terminate the local part of the path at this spot. */
          *end = '\0';
          /* Keep track of the system part of the path, which starts after that. */
          sys = end + 1;
        }

        /* Chop up local path. */
        status = parse_path(s_par_path.loc_path, begin);
        if (eOK == status && 0 != sys) {
          /* Chop up system path if it exists. */
          status = parse_path(s_par_path.sys_path, sys);
        }
      }
    }
  }

  if (eOK == status) {
    s_par_path.init = 1;
  } else {
    destroy_paths();
  }
  return status;
}

static void ape_io_atexit(void) {
  destroy_paths();
}

static void destroy_paths(void) {
  /* Reset flag which indicates init is done. */
  s_par_path.init = 0;

  /* Free the path lists themselves, and reset the pointers. */
  ape_list_destroy(s_par_path.sys_path); s_par_path.sys_path = 0;
  ape_list_destroy(s_par_path.loc_path); s_par_path.loc_path = 0;

  /* Free the string which holds the value of pfiles. */
  free(s_par_path.pfiles); s_par_path.pfiles = 0;
}

static int parse_path(ApeList * list, char * begin) {
  int status = eOK;
  if (0 != list) {
    if (0 != begin && '\0' != *begin) {
      if (0 != ape_list_append(list, begin)) {
        char * end = 0;
        for (end = strstr(begin, QC); 0 != end; end = strstr(begin, QC)) {
          *end = '\0';
          begin = end + 1;
          if (0 == ape_list_append(list, begin)) {
            status = eDynAllocFailed;
            break;
          }
        }
      } else {
        status = eDynAllocFailed;
      }
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

static void destroy_par_list(ApeList * par_list) {
  ApeListIterator itor = 0;
  if (0 != par_list) {
    for (itor = ape_list_end(par_list); itor != ape_list_begin(par_list);) {
      itor = ape_list_prev(itor);
      ape_par_destroy((ApePar *) ape_list_get(itor));
    }
    ape_list_destroy(par_list); par_list = 0;
  }
}

/* Type of function which operates on parameters from two lists. */
typedef int (*par_operation_type)(ApeList *, ApeListIterator *, ApeList *, ApeListIterator *, char);

/* Utility which accepts two parameter files. It iterates over the parameters in the first "master" file,
   looks up each parameter in the second file, and then calls the supplied function on behalf of each pair
   of parameters. */
static int apply_operation(ApeParFile * master_file, ApeParFile * par_file, par_operation_type operation);

int ape_io_clone_file(ApeParFile * orig_file, ApeParFile ** clone_file) {
  int status = eOK;

  if (0 != clone_file) {
    *clone_file = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == orig_file) {
    status = eNullPointer;
  }

  if (eOK == status) {
    *clone_file = (ApeParFile *) calloc(1, sizeof(ApeParFile));
    if (0 == *clone_file) status = eDynAllocFailed;
  }

  if (eOK == status && 0 != orig_file->file_name) {
    /* Copy the file name. */
    status = ape_util_copy_string(orig_file->file_name, &(*clone_file)->file_name);
  }

  if (eOK == status && 0 != orig_file->par_cont) {
    (*clone_file)->par_cont = ape_list_create();
    if (0 == (*clone_file)->par_cont) status = eDynAllocFailed;
  }

  if (eOK == status && 0 != orig_file->par_cont) {
    ApeListIterator itor = ape_list_begin(orig_file->par_cont);
    ApeListIterator end = ape_list_end(orig_file->par_cont);
    for (; itor != end; itor = ape_list_next(itor)) {
      ApePar * clone_par = 0;
      status = ape_par_clone((ApePar *) ape_list_get(itor), &clone_par);
      if (eOK == status) {
        ape_list_append((*clone_file)->par_cont, clone_par);
      } else {
        break;
      }
    }
  }

  if (eOK == status) {
    /* Copy read-only flag. */
    (*clone_file)->read_only = orig_file->read_only;
  }

  if (eOK == status) {
    /* Copy the file modification time. */
    (*clone_file)->mod_time = orig_file->mod_time;
  }

  return status;
}

void ape_io_destroy_file(ApeParFile * par_file) {
  if (0 != par_file) {
    destroy_par_list(par_file->par_cont); par_file->par_cont = 0;
    free(par_file->file_name); par_file->file_name = 0;
    free(par_file);
  }
}

int ape_io_find_par(const char * par_name, ApeParFile * par_file, ApeListIterator * par_itor) {
  int status = eOK;
  ApeList * par_cont = 0;

  /* Check output argument. */
  if (0 != par_itor) {
    /* Initialize output. */
    *par_itor = 0;
  } else {
    status = eNullPointer;
  }

  /* Check input arguments. */
  if (eOK == status && (0 == par_name || 0 == par_file)) {
    status = eNullPointer;
  }

  /* Get the container of parameters from the file, and make sure it is valid. */
  if (eOK == status) {
    par_cont = par_file->par_cont;
    if (0 == par_cont) status = eNullPointer;
  }

  if (eOK == status) {
    ApePar * found = 0;
    /* Iterate over the parameter container. */
    ApeListIterator itor = ape_list_begin(par_cont);
    ApeListIterator end = ape_list_end(par_cont);
    for (; itor != end; itor = ape_list_next(itor)) {
      char * current_name = 0;

      /* Get the current parameter, and extract its name as a string. */
      found = (ApePar *) ape_list_get(itor);
      status = ape_par_get_field(found, eName, &current_name);

      /* See if the name of the current parameter matches the parameter being sought. */
      if (eOK == status && 0 == ape_util_strcmp(current_name, par_name, 1)) {
        /* Parameter was found; clean up and stop. */
        free(current_name); current_name = 0;
        break;
      }
      /* Parameter was not found; clean up and go on. */
      free(current_name); current_name = 0;
    }

    if (eOK == status) {
      /* No matter what, return the current iterator position. */
      *par_itor = itor;

      /* Indicate if parameter was not found. */
      if (itor == end) {
        status = eParNotFound;
      }
    }
  }

  return status;
}

int ape_io_get_par_cont(ApeParFile * par_file, ApeList ** par_cont) {
  int status = eOK;

  if (0 != par_cont) {
    *par_cont = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == par_file) {
    status = eNullPointer;
  }

  if (eOK == status) {
    *par_cont = par_file->par_cont;
  }

  return status;
}

int ape_io_get_file_name(const ApeParFile * par_file, const char ** file_name) {
  int status = eOK;
  if (0 != file_name) {
    *file_name = 0;
    if (0 != par_file) {
      *file_name = par_file->file_name;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

int ape_io_set_file_name(ApeParFile * par_file, const char * file_name) {
  int status = eOK;

  if (0 != par_file) {
    free(par_file->file_name); par_file->file_name = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == file_name) {
    status = eNullPointer;
  }

  if (eOK == status) {
    status = ape_util_copy_string(file_name, &par_file->file_name);
  }

  return status;
}

int ape_io_get_default_mode(ApeParFile * par_file, char * mode) {
  int status = eOK;
  ApeListIterator itor = 0;
  ApePar * par = 0;
  char * mode_string = 0;

  /* Check arguments. */
  if (0 != mode) {
    memset(mode, '\0', APE_PAR_MODE_CODE_LEN);
  } else {
    status = eNullPointer;
  }
  if (eOK == status && 0 == par_file) status = eNullPointer;

  if (eOK == status) {
    /* Get the mode parameter. */
    status = ape_io_find_par("mode", par_file, &itor);
  }

  if (eOK == status) {
    par = (ApePar *) ape_list_get(itor);

    /* Get the value of the mode parameter. */
    status = ape_par_get_field(par, eValue, &mode_string);
  }

  if (eOK == status) {
    strncpy(mode, mode_string, APE_PAR_MODE_CODE_LEN);
  }

  /* Clean up. */
  free(mode_string); mode_string = 0;

  return status;
}

static int par_prune(ApeList * master_cont, ApeListIterator * master_itor, ApeList * par_cont, ApeListIterator * par_itor,
  char master_more_recent) {
  int status = eOK;
  if (*master_itor != ape_list_end(master_cont) && (*par_itor == ape_list_end(par_cont))) {
    ApePar * par = (ApePar *) ape_list_get(*master_itor);

    /* The master container has a parameter the other container does not have. Delete the one in the master container. */
    *master_itor = ape_list_remove_entry(master_cont, *master_itor);
    if (*master_itor != ape_list_end(master_cont)) *master_itor = ape_list_next(*master_itor);
    ape_par_destroy(par);
  }
  return status;
}

static int par_synch(ApeList * master_cont, ApeListIterator * master_itor, ApeList * par_cont, ApeListIterator * par_itor,
  char master_more_recent) {
  int status = eOK;
  if (*master_itor != ape_list_end(master_cont) && (*par_itor != ape_list_end(par_cont))) {
    ApePar * master = (ApePar *) ape_list_get(*master_itor);
    ApePar * other = (ApePar *) ape_list_get(*par_itor);
    char master_type[APE_PAR_TYPE_CODE_LEN] = "";
    char other_type[APE_PAR_TYPE_CODE_LEN] = "";
    char master_mode[APE_PAR_MODE_CODE_LEN] = "";
    char other_mode[APE_PAR_MODE_CODE_LEN] = "";
    char ** master_enum = 0;
    char ** other_enum = 0;
    char * master_max = 0;
    char * other_max = 0;

    /* Check the overall state of the local parameter. If it's corrupt somehow, just use the master parameter for everything. */
    status = ape_par_check(other, 0);
    if (eOK != status) {
      ApePar * new_other = 0;
      status = ape_par_clone(master, &new_other); 
      if (eOK == status) {
        ape_list_set(*par_itor, new_other);
        ape_par_destroy(other);
      }
      return status;
    }

    /* Compare types of the two parameters. */
    if (eOK == status) status = ape_par_get_type(master, master_type);
    if (eOK == status) {
      status = ape_par_get_type(other, other_type);
    }
    if (eOK == status && 0 != ape_util_strcmp(master_type, other_type, 1)) {
      char * master_type_string = 0;

      /* If any difference, the master parameter is considered correct, so set other's type to match. */
      status = ape_par_get_field(master, eType, &master_type_string);
      if (eOK == status) {
        status = ape_par_set_field(other, eType, master_type_string);
      }
      free(master_type_string); master_type_string = 0;
    }

    /* Compare modes of the two parameters. */
    status = ape_par_get_mode(master, master_mode);
    if (eOK == status) {
      status = ape_par_get_mode(other, other_mode);
    }
    if (eOK == status && 0 != ape_util_strcmp(master_mode, other_mode, 1)) {
      char * master_mode_string = 0;

      /* If any difference, the master parameter is considered correct, so set other's mode to match. */
      status = ape_par_get_field(master, eMode, &master_mode_string);
      if (eOK == status) {
        status = ape_par_set_field(other, eMode, master_mode_string);
      }
      free(master_mode_string); master_mode_string = 0;
    }

    /* Compare min of the two parameters; use the enum version which is more general. */
    status = ape_par_get_enum_string(master, &master_enum);
    if (eOK == status || eRangeNoEnum == status) {
      int other_status = ape_par_get_enum_string(other, &other_enum);
      if (status != other_status || 0 != ape_util_cmp_string_array((const char **) master_enum, (const char **) other_enum, 1)) {
        /* Although it was more convenient and more forgiving to compare using the array, it's easier
           to get/set a single string. */
        char * master_enum_string = 0;

        /* If any difference, the master parameter is considered correct, so set other's min to match. */
        status = ape_par_get_field(master, eMin, &master_enum_string);
        if (eOK == status) {
          status = ape_par_set_field(other, eMin, master_enum_string);
        }
        free(master_enum_string); master_enum_string = 0;
      }
    }

    /* Clean up. */
    ape_util_free_string_array(other_enum); other_enum = 0;
    ape_util_free_string_array(master_enum); master_enum = 0;

    /* Compare maxes of the two parameters. */
    status = ape_par_get_field(master, eMax, &master_max);
    if (eOK == status) {
      status = ape_par_get_field(other, eMax, &other_max);
    }
    if (eOK == status && 0 != ape_util_strcmp(master_max, other_max, 1)) {
      char * master_max_string = 0;

      /* If any difference, the master parameter is considered correct, so set other's max to match. */
      status = ape_par_get_field(master, eMax, &master_max_string);
      if (eOK == status) {
        status = ape_par_set_field(other, eMax, master_max_string);
      }
      free(master_max_string); master_max_string = 0;
    }

    /* Clean up. */
    free(other_max); other_max = 0;
    free(master_max); master_max = 0;

    /* Check how to merge the values. The logic for this is different from the other fields. If
       the master file has a more recent time stamp than the other file, hidden parameter values only should
       be copied from the master to the other file. Otherwise, parameter values should be copied
       from the other file to the master file. */
    if (eOK == status) {
      char * value = 0;
      if (0 != master_more_recent && 'h' == *master_mode) {
        /* Copy value field from master parameter to other parameter. */
        status = ape_par_get_field(master, eValue, &value);
        if (eOK == status) {
          status = ape_par_set_field(other, eValue, value);
        }
      } else {
        /* Copy value field from other parameter to master parameter. */
        status = ape_par_get_field(other, eValue, &value);
        if (eOK == status) {
          status = ape_par_set_field(master, eValue, value);
        }
      }
      free(value); value = 0;
    }

  } else if (*master_itor != ape_list_end(master_cont)) {
    /* The master container has a parameter the other container does not have. Flag this so that the master
       container will be used (with modification) instead of the other container. */
    status = eNewParameter;
  } else if (*par_itor != ape_list_end(par_cont)) {
    /* The other container has a parameter the master container does not have. Delete the one in the other container. */
    ApePar * par = (ApePar *) ape_list_get(*par_itor);
    ape_list_remove_entry(par_cont, *par_itor);
    ape_par_destroy(par);
  }
  return status;
}

/* Compare parameter groups; check whether the parameters come in the same order and have the
   same names, types, modes, min and max. Also check the time stamps. Skip comments and blank lines. */
/* TODO Clean this up. Add some helpers to ape_par interface to make it easier to compare parameters
   and avoid some of the repetitious and spaghetti-like logic below. */
static int rough_compare_files(ApeParFile * par_file1, ApeParFile * par_file2) {
  int status = eOK;
  ApeList * par_cont1 = 0;
  ApeList * par_cont2 = 0;
  ApeListIterator itor1 = 0;
  ApeListIterator itor2 = 0;
  ApeListIterator end1 = 0;
  ApeListIterator end2 = 0;

  status = ape_io_get_par_cont(par_file1, &par_cont1);
  if (eOK == status) {
    status = ape_io_get_par_cont(par_file2, &par_cont2);
  }

  /* The master file is par_file1. Check whether it is more recent than the other file and return non-zero
     status if so. */
  if (eOK == status && par_file1->mod_time > par_file2->mod_time) {
    status = -1;
  }

  if (eOK != status) {
    status = -1;
    return status;
  }

  itor1 = ape_list_begin(par_cont1);
  itor2 = ape_list_begin(par_cont2);

  end1 = ape_list_end(par_cont1);
  end2 = ape_list_end(par_cont2);

  /* Go through containers looking for parameters which are different in some significant way. */
  while (eOK == status && (itor1 != end1 || itor2 != end2)) {
    static ParFieldId s_field[] = { eType, eMode, eMin, eMax };
    static size_t s_num_field = sizeof(s_field) / sizeof(s_field[0]);
    int status1 = eOK;
    int status2 = eOK;
    ApePar * par1 = itor1 != end1 ? (ApePar *) ape_list_get(itor1) : 0;
    ApePar * par2 = itor2 != end2 ? (ApePar *) ape_list_get(itor2) : 0;
    size_t ii = 0;
    char * field1 = 0;
    char * field2 = 0;

    /* Skip parameter if its iterator already points to the end of the container. */
    if (eOK == status && itor1 != end1) {
      /* Check parameter to see if it has a name. */
      status1 = ape_par_get_name(par1, &field1);
      if (eUnnamedPar == status1) {
        /* No name so do not compare further; just go on. */
        free(field1); field1 = 0;
        itor1 = ape_list_next(itor1);
        continue;
      }
    }
    if (eOK == status && itor2 != end2) {
      /* Check parameter to see if it has a name. */
      status2 = ape_par_get_name(par2, &field2);
      if (eUnnamedPar == status2) {
        /* No name so do not compare further; just go on. */
        free(field2); field2 = 0;
        free(field1); field1 = 0;
        itor2 = ape_list_next(itor2);
        continue;
      }
    }
    /* Make sure that if there was a problem both pars had it. */
    if (eOK == status && status1 != status2) status = -1;

    /* Parameters are considered different if one of them is at the end of its container
       and the other is not, or if they point to parameters with different names. */
    if (eOK == status && (itor1 == end1 || itor2 == end2 || 0 != ape_util_strcmp(field1, field2, 1))) status = -1;

    /* Clean up name field. */
    free(field2); field2 = 0;
    free(field1); field1 = 0;

    for (ii = 0; eOK == status && eOK == status1 && eOK == status2 && ii != s_num_field; ++ii) {
      /* At this point, both parameters have names that agree, so insist that the remaining
         critical fields also match. */
      status1 = ape_par_get_field(par1, s_field[ii], &field1);
      status2 = ape_par_get_field(par2, s_field[ii], &field2);
      /* Only compare fields if both were obtained cleanly. */
      if (eOK == status1 && eOK == status2) {
        if (0 != ape_util_strcmp(field1, field2, 1)) status = -1;
      } else if (status1 != status2) {
        /* If the same error was not encountered for both fields, there is disagreement. */
        status = -1;
      }
      free(field2); field2 = 0;
      free(field1); field1 = 0;
    }

    /* For better or worse, done comparing these two parameters so go on. */
    if (itor1 != end1) itor1 = ape_list_next(itor1);
    if (itor2 != end2) itor2 = ape_list_next(itor2);
  }
  return status;
}

int ape_io_merge_par_files(ApeParFile * sys_par_file, ApeParFile * loc_par_file, ApeParFile ** merged_par_file) {
  int status = eOK;
  ApeParFile *loc_ptr = 0;
  ApeParFile *sys_ptr = 0;

  /* Check output variable. */
  if (0 == merged_par_file) {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Initialize output. */
    *merged_par_file = 0;

    /* Check input. */
    if (0 == sys_par_file && 0 == loc_par_file) {
      status = eNullPointer;
    }
  }

  if (eOK == status && 0 != loc_par_file) {
    /* Make copy of input local par file. */
    status = ape_io_clone_file(loc_par_file, &loc_ptr);
  }

  if (eOK == status) {
    status = rough_compare_files(sys_par_file, loc_par_file);

    if (eOK == status) {
      *merged_par_file = loc_ptr;
      return status;
    } else if (-1 == status) {
      status = eOK;
    }
  }

  if (eOK == status && 0 != sys_par_file) {
    /* Make copy of input system par file. */
    status = ape_io_clone_file(sys_par_file, &sys_ptr);
  }

  if (eOK == status) {
    if (0 != loc_ptr && 0 != sys_ptr) {
      /* Synchronize the two cloned parameter files. */
      status = apply_operation(sys_ptr, loc_ptr, par_synch);
      if (eNewParameter == status) {
        /* Swap file names, because merged file should have the same name as the local file. */
        char * sys_name = sys_ptr->file_name;
        sys_ptr->file_name = loc_ptr->file_name;
        loc_ptr->file_name = sys_name;

        /* Assign to merged file the system file (with local name), and destroy the local file. */
        *merged_par_file = sys_ptr; sys_ptr = 0;
        ape_io_destroy_file(loc_ptr); loc_ptr = 0;

        /* Reset status, since caller doesn't care which merged file we're passing back. */
        status = eOK;
      } else if (eOK == status) {
        /* Remove from the local clone any/all parameters which are not present in the system file. */
        status = apply_operation(loc_ptr, sys_ptr, par_prune);
        if (eOK == status) {
          *merged_par_file = loc_ptr; loc_ptr = 0;
          ape_io_destroy_file(sys_ptr); sys_ptr = 0;
        }
      }
    } else if (0 != loc_ptr) {
      *merged_par_file = loc_ptr; loc_ptr = 0;
    } else if (0 != sys_ptr) {
      *merged_par_file = sys_ptr; sys_ptr = 0;
      (*merged_par_file)->read_only = 1;
    }
  }

  if (eOK != status) {
    *merged_par_file = 0;
    ape_io_destroy_file(sys_ptr); sys_ptr = 0;
    ape_io_destroy_file(loc_ptr); loc_ptr = 0;
  }
  return status;
}

/* TODO Find a way to achieve this without creating another static variable. */
static char s_default_mode[APE_PAR_MODE_CODE_LEN] = "";

static int par_revert(ApeList * master_cont, ApeListIterator * master_itor, ApeList * par_cont, ApeListIterator * par_itor,
  char master_more_recent) {
  int status = eOK;
  if (*master_itor != ape_list_end(master_cont)) {
    ApePar * master = (ApePar *) ape_list_get(*master_itor);
    ApePar * par = 0;
    char mode[APE_PAR_MODE_CODE_LEN] = "";
    /* Determine the effective mode of the parameter. */
    status = ape_par_get_eff_mode(master, s_default_mode, mode);
    if (eOK == status && 'l' != mode[1]) {
      /* Parameter is not learned, so restore previous value. */
      if (*par_itor != ape_list_end(par_cont)) {
        char * value = 0;
        par = (ApePar *) ape_list_get(*par_itor);
        status = ape_par_get_field(par, eValue, &value);
        if (eOK == status) {
          status = ape_par_set_field(master, eValue, value);
        }
        free(value); value = 0;
      } else {
        status = eParNotFound;
      }
    }
  }
  return status;
  
}

int ape_io_revert_unlearned(ApeParFile * current, ApeParFile * previous) {
  int status = eOK;

  /* Check arguments. */
  if (0 == current || 0 == previous) status = eNullPointer;

  if (eOK == status) {
    /* Get default mode from current parameter file, and assume learned behavior if default mode is not defined. */
    status = ape_io_get_default_mode(current, s_default_mode);
    if (eParNotFound == status) {
      status = eOK;
      strcpy(s_default_mode, "ql");
    }
  }

  if (eOK == status) {
    status = apply_operation(current, previous, par_revert);
  }

  return status;
}

int ape_io_get_pfiles(char ** pfiles) {
  int status = eOK;

  if (0 != pfiles) {
    *pfiles = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == s_par_path.init) {
    status = create_paths(0);
  }

  if (eOK == status) {
    status = ape_util_copy_string(s_par_path.pfiles, pfiles);
  }

  return status;
}

int ape_io_set_pfiles(const char * pfiles) {
  int status = eOK;
  /* Clear up previous value for pfiles, and for associated path structures. */
  destroy_paths();

  /* Create new pfiles. */
  create_paths(pfiles);

  return status;
}

int ape_io_get_loc_path(ApeList ** loc_path) {
  int status = eOK;

  if (0 != loc_path) {
    *loc_path = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == s_par_path.init) {
    status = create_paths(0);
  }

  if (eOK == status) {
    *loc_path = s_par_path.loc_path;
  }

  return status;
}

int ape_io_get_sys_path(ApeList ** sys_path) {
  int status = eOK;

  if (0 != sys_path) {
    *sys_path = s_par_path.sys_path;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == s_par_path.init) {
    status = create_paths(0);
  }

  if (eOK == status) {
    *sys_path = s_par_path.sys_path;
  }

  return status;
}

static time_t get_mod_time(const char * file_name) {
  time_t mod_time = 0;
  int status = 0;
#ifdef WIN32
  struct _stat file_stat;
  status = _stat(file_name, &file_stat);
#else
  struct stat file_stat;
  status = stat(file_name, &file_stat);
#endif
  if (0 != status) {
    /* If stat failed, ignore modification time (set to 0), but warn about it. */
    ape_msg_debug("get_mod_time: unable to get modification time of file \"%s\", errno was %d.\n", file_name, errno);
  } else {
    mod_time = file_stat.st_mtime;
  }
  return mod_time;
}

#define APE_BUF_SIZE (8192u)

int ape_io_read(const char * file_name, ApeParFile ** par_file) {
  ApeList ** cont = 0;
  int status = eOK;
  int parse_status = eOK;
  FILE * fp = 0;

  /* Make sure output parameter file pointer is valid. */
  if (0 != par_file) {
    /* Initialize output parameter file object pointer. */
    *par_file = 0;
  } else {
    ape_msg_debug("ape_io_read: output parameter file pointer is null.\n");
    status = eInvalidArgument;
  }

  /* Make sure file name is valid. */
  if (eOK == status && 0 != file_name) {
    /* Skip leading white space. */
    while (0 != isspace(*file_name)) ++file_name;
  } else {
    ape_msg_debug("ape_io_read: parameter file name pointer is null.\n");
    status = eNullPointer;
  }

  /* Make sure file name is not empty. */
  if (eOK == status && '\0' == *file_name) {
    ape_msg_debug("ape_io_read: parameter file name is blank.\n");
    status = eInvalidArgument;
  }

  if (eOK == status) {
    /* Open the file. */
    fp = fopen(file_name, "r");
    /* Make sure file actually opened. */
    if (0 == fp) {
      ape_msg_debug("ape_io_read: unable to open file \"%s\" for reading.\n", file_name);
      status = eFileReadError;
    }
  }

  if (eOK == status) {
    /* Create output file object. */
    *par_file = (ApeParFile *) calloc(1, sizeof(ApeParFile));
    if (0 == *par_file) {
      ape_msg_debug("ape_io_read: unable to allocate parameter file object.\n");
      status = eDynAllocFailed;
    }
    cont = &(*par_file)->par_cont;
  }

  if (eOK == status) {
    (*par_file)->mod_time = get_mod_time(file_name);
  }

  if (eOK == status) {
    /* Save file name. */
    status = ape_util_copy_string(file_name, &(*par_file)->file_name);
  }

  if (eOK == status) {
    /* Create output container. */
    *cont = ape_list_create();
    /* Make sure container was created. */
    if (0 == *cont) {
      ape_msg_debug("ape_io_read: unable to allocate container of parameters.\n");
      status = eDynAllocFailed;
    }
  }

  /* Iterate over file. */
  while (eOK == status && 0 == feof(fp) && 0 == ferror(fp)) {
    char line[APE_BUF_SIZE] = "";
    static const size_t ape_buf_size = APE_BUF_SIZE;

    /* Read a line. If this returns 0, the read was not useful (error or nothing left in the file). */
    if (0 != fgets(line, (int) ape_buf_size, fp)) {
      static const char newline = '\n';
      static const char linefeed = 10;
      static const char carriage_return = 13;
      ApePar * par = 0;
      char * cp = 0;
      char newline_found = 0;
      size_t line_len = strlen(line);

      /* Look for and wipe out newlines (and other line terminators, in case the parameter file was
         written on a different OS. */
      for (cp = line + line_len; cp > line && (*(cp-1) == newline || *(cp-1) == linefeed || *(cp-1) == carriage_return); --cp) {
        *(cp-1) = '\0';
        newline_found = 1;
      }
      /* If no newline was found, and the buffer is full, the line is too long. This is a problem but just
         flag it and keep going. */
      if ((line_len == ape_buf_size - 1) && 0 == newline_found) {
        char next_char = '\0';
        char start_line[64] = "";
        strncpy(start_line, line, sizeof(start_line)/sizeof(start_line[0]) - 1);
        ape_msg_debug("ape_io_read: line starting with %s is longer than maximum length (%u), truncating.\n",
          start_line, ape_buf_size);
        /* Skip over everything on the line until the next line. */
        while ((0 == feof(fp)) && (0 == ferror(fp)) &&
          (newline != next_char) && (linefeed != next_char) && (carriage_return != next_char)) {
          next_char = fgetc(fp);
        }
        parse_status = eLineTooLong;
      }

      /* Extract the parameter from the container. */
      status = ape_par_create(&par, line);
      /* Make sure parameter was created. */
      if (eOK == status) {
        /* Add parameter to container. */
        ApeListIterator itor = ape_list_append(*cont, par);
        if (ape_list_end(*cont) == itor) {
          ape_msg_debug("ape_io_read: unable to add to container parameter \"%s\" .\n", line);
          status = eListError;
        }
      } else {
        ape_msg_debug("ape_io_read: unable to parse line \"%s\" into a parameter.\n", line);
        /* Status was already set to non-0 value by ape_par_create, so don't change it. */
      }
    }
  }
  if (eOK == status && 0 != ferror(fp)) {
    ape_msg_error("error reading from file \"%s\".\n", file_name);
    status = eFileReadError;
  }
  /* Clean up. */
  if (eOK != status) {
    ape_io_destroy_file(*par_file); *par_file = 0;
  }
  if(0 != fp) fclose(fp);

  /* Report parse errors but don't destroy the output. */
  if (eOK == status) status = parse_status;

  return status;
}

static char * append_path(const char * dir_name, const char * file_name) {
  char * dir_plus_sep = 0;
  char * full_file_name = 0;
  int status = eOK;

  if (0 != dir_name && 0 != file_name) {
    /* TODO: Prevent adding extra separators, and doing 2 concatenations? */
    status = ape_util_cat_string(dir_name, QD, &dir_plus_sep);
    if (eOK == status) {
      status = ape_util_cat_string(dir_plus_sep, file_name, &full_file_name);
    }
    free(dir_plus_sep);
  }
  return full_file_name;
}

static char is_blank(const char * s) {
  if (0 == s) return 1;
  for (; '\0' != *s; ++s) {
    if (0 == isspace(*s)) return 0;
  }
  return 1;
}

static int apply_operation(ApeParFile * master_file, ApeParFile * par_file, par_operation_type operation) {
  /* Note: no checks of inputs are performed; this is purely internal for now. */
  int status = eOK;
  int operation_status = eOK;
  ApeList * master_cont = master_file->par_cont;
  ApeList * par_cont = par_file->par_cont;
  char master_more_recent = 1;
  ApeListIterator master_itor = 0;

  /* If master file and other file have valid modification times, use them to determine whether
     master file is more recent. */
  if ((0 != master_file->mod_time) && (0 != par_file->mod_time)) {
    if (master_file->mod_time < par_file->mod_time) master_more_recent = 0;
  }

  /* Update the virtual file modification times so that future operations will work correctly. */
  if (master_more_recent) {
    par_file->mod_time = master_file->mod_time;
  } else {
    master_file->mod_time = par_file->mod_time;
  }

  /* Iterate over the master parameter container. */
  master_itor = ape_list_end(master_cont);
  while (master_itor != ape_list_begin(master_cont)) {
    char * name = 0;
    ApePar * par = 0;

    /* Get the current parameter. */
    master_itor = ape_list_prev(master_itor);
    par = (ApePar *) ape_list_get(master_itor);

    /* Get the name of the parameter. */
    status = ape_par_get_field(par, eName, &name);

    /* If the name is not blank (i.e. if this is an actual parameter, not a comment or blank line),
       look it up in the other parameter file. */
    if (eOK == status && 0 == is_blank(name)) {
      ApeListIterator par_itor = 0;
      /* Look for this parameter in the other list. Ignore the status. */
      ape_io_find_par(name, par_file, &par_itor);

      /* Call the operation for this set of iterators, even if the matching parameter was not found in the other file. */
      status = operation(master_cont, &master_itor, par_cont, &par_itor, master_more_recent);

      /* Preserve the first non-OK status in operation_status. */
      operation_status = (eOK == operation_status) ? status : operation_status;
    }

    /* Clean up. */
    free(name); name = 0;
  }

  return operation_status;
}

int ape_io_read_file_path(const char * file_name, ApeList * path, ApeParFile ** par_file) {
  int status = eFileNotFound;
  if (0 != par_file) {
    /* Start by initializing output to 0, in case we can't find the file. */
    *par_file = 0;

    if (0 != file_name) {
      ApeListIterator itor = 0;
/* TODO: Why does this change the unit test output? */
/*      ApeListIterator end = ape_list_end(path); */

      /* Iterate over the paths, try to open the given file in each directory until there is success. */
      for (itor = ape_list_begin(path); eOK != status && itor != ape_list_end(path); itor = ape_list_next(itor)) {
        /* Get name of next directory from the path. */
        const char * dir = (const char *) ape_list_get(itor);
        char * full_file_name;

        /* It is an error if at any point a component of a path is null. */
        if (0 == dir) {
          status = eNullPointer;
          break;
        } else if (0 != is_blank(dir)) {
          /* Blank path components are OK -- just skip them. */
          status = eFileNotFound;
          continue;
        }

        /* Compose the full file name from the directory name and the given file_name. */
        full_file_name = append_path(dir, file_name);

        /* Attempt to read the list of parameters from this file. */
        status = ape_io_read(full_file_name, par_file);

        /* Clean up. */
        free(full_file_name); full_file_name = 0;
      }
    } else {
      status = eNullPointer;
    }
  } else {
    status = eNullPointer;
  }
  return status;
}

static char * create_tmp_file_name(const char * file_name) {
  char * tmp_file_name = 0;
  char host_name[FILENAME_MAX] = "";
  char pid[64] = "";
  if (0 != file_name) {
#ifdef WIN32
    sprintf(pid, "%lu", GetCurrentProcessId());
#else
    sprintf(pid, "%lu", (unsigned long) getpid());

/* On Cygwin, gethostname tries to go out to the network and hangs if firewall installed. */
#ifndef __CYGWIN__
    gethostname(host_name, FILENAME_MAX);
#endif

#endif 
    /* Allocate space for base name + host name + process id + 2 dashes + terminating 0. */
    tmp_file_name = (char *) calloc(strlen(file_name) + strlen(host_name) + strlen(pid) + 3, sizeof(char));
    if (0 != tmp_file_name)
      sprintf(tmp_file_name, "%s-%s-%s", file_name, host_name, pid);
  }
  return tmp_file_name;
}

int ape_io_write(ApeParFile * par_file, char force_write) {
  int status = eOK;
  ApeList * cont = 0;
  const char * file_name = 0;

  /* Check argument. */
  if (0 == par_file || 0 == par_file->file_name) {
    status = eNullPointer;
  } else if (0 != is_blank(par_file->file_name)) {
    status = eInvalidArgument;
  } else {
    file_name = par_file->file_name;
    cont = par_file->par_cont;
  }

  if (eOK == status && (0 != force_write || 0 == par_file->read_only)) {
    /* Skip leading white space. */
    while (0 != isspace(*file_name)) ++file_name;

    /* Check whether file name is empty. */
    if ('\0' != *file_name) {
      FILE * fp = 0;

      /* Devise name of temporary file from base file name, hostname and pid. */
      char * tmp_file_name = create_tmp_file_name(file_name);
      if (0 == tmp_file_name) status = eDynAllocFailed;

      if (eOK == status) {
        /* Open the file. */
        fp = fopen(tmp_file_name, "w");
        if (0 == fp) {
          ape_msg_error("unable to open file \"%s\" for writing.\n", tmp_file_name);
          status = eFileWriteError;
        }
      }

      if (eOK == status) {
        if (0 != cont) {
          /* Container contains some data, so write it to the file. */
          ApeListIterator itor = 0;
          ApeListIterator end = ape_list_end(cont);

          /* Iterate over container. */
          for (itor = ape_list_begin(cont); itor != end; itor = ape_list_next(itor)) {
            /* Extract the parameter from the container. */
            ApePar * par = (ApePar *) ape_list_get(itor);

            if (0 != par) {
              /* Put the parameter into string format. */
              char * line = ape_par_get_line(par);

              /* Write the parameter into the file. */
              status = fprintf(fp, "%s\n", line);

              /* Free the parameter representation. */
              free(line);

              /* Check for error while writing to file. */
              if (EOF != status) {
                status = eOK;
              } else {
                ape_msg_error("error writing to file \"%s\".\n", tmp_file_name);
                status = eFileWriteError;
                break;
              }
            } else {
              ape_msg_debug("ape_io_write: a list element is a null pointer.\n");
              status = eNullPointer;
            }
          }
        } else {
          ape_msg_debug("ape_io_write: container of parameters is empty.\n");
          status = eInvalidArgument;
        }
      }
      if (0 != fp) fclose(fp); fp = 0;

      if (0 != tmp_file_name) {
        /* If all went well thus far, clobber the output file. */
        if (eOK == status) remove(file_name);

        /* Atomically move the temporary file over top of the real desired file name. Report errors. */
        if (eOK == status && 0 != rename(tmp_file_name, file_name)) status = eFileRenameError;

        /* In case rename failed in such a way that tmp_file_name still exists, remove tmp_file_name. */
        remove(tmp_file_name);
      }

      /* Clean up memory. */
      free(tmp_file_name); tmp_file_name = 0;

    } else {
      ape_msg_debug("ape_io_write: parameter file name is blank.\n");
      status = eInvalidArgument;
    }
  }

  return status;
}

/* Codes to classify arguments. */
enum { UNKNOWN = 0, PAR = 1, EQUALS = 2, VALUE = 4, PLUS = 8, MINUS = 16, DONE = 32 };

static int classify_arg(char ** par_name, int end_par, int arg_index, char * arg_p, unsigned char * type, int * par_index,
  char ** value) {
  int index;
  int status = eOK;
  int match = 0;
  char * arg_save = 0;
  char * arg_orig = 0;

  /* Skip leading white space in argument. */
  while (0 != isspace(*arg_p)) ++arg_p;

  /* Save original argument for use (if needed) in error message later. */
  arg_orig = arg_p;

  /* Initially classify as unknown, with parameter index one-past-last par, and value is the first significant character. */
  *type = UNKNOWN;
  *par_index = end_par;
  *value = arg_p;

  /* Check argument against array of parameters to see if argument is of one of the forms PAR, PAR =, or PAR = VALUE. */
  for (index = 0; index < end_par; ++index) {
    char * par_p = 0;

    /* Get pointer to the name of this parameter. */
    par_p = par_name[index];

    /* Skip empty parameters (with no name). */
    if (0 == par_p || '\0' == *par_p) continue;

    /* Always start looking for parameter at the first significant character. */
    arg_p = *value;

    /* Compare parameter to argument. Break when a difference is detected or the end of either string. */
    for (; '\0' != *arg_p && '\0' != *par_p && toupper(*arg_p) == toupper(*par_p); ++arg_p, ++par_p) {}

    /* If (all or some of) the parameter name was consumed, this argument *may* match this parameter. */
    if ('\0' == *par_p || par_p != par_name[index]) {
      /* To match for sure the next character in the argument must be either white space, equals sign, plus or minus. */
      /* TODO check at this point for dereferencing parameter fields, e.g. p_name, p_type etc. */
      while (0 != isspace(*arg_p)) ++arg_p;
      if ('\0' == *arg_p || '=' == *arg_p || '+' == *arg_p || '-' == *arg_p) {
        /* This argument holds the name of a parameter. */
        *type |= PAR;
        *par_index = index;
        /* The value should be reset to point to the first significant character after the parameter name. */
        arg_save = arg_p;
        /* Break if an exact match has been found, otherwise add to the number of potential matches. */
        if ('\0' == *par_p) { match = 1; break; } else ++match;
      }
    }
  }

  /* If at least one match was found, restore value & arg_p and proceed.
     If more than one possible match was found, flag as ambiguous. */
  if (0 != match) {
    *value = arg_save; arg_p = arg_save;
    if (1 < match) status = eAmbiguousParName;
  }

  /* Look for equals, plus or minus sign. */
  if ('=' == *arg_p) {
    *type |= EQUALS;

    /* Look past equals sign as well as following white space. */
    ++arg_p;
    while (0 != isspace(*arg_p)) ++arg_p;
  } else if ('+' == *arg_p) {
    *type |= PLUS;

    /* Look past plus sign as well as following white space. */
    ++arg_p;
    while (0 != isspace(*arg_p)) ++arg_p;
  } else if ('-' == *arg_p) {
    *type |= MINUS;

    /* Look past minus sign as well as following white space. */
    ++arg_p;
    while (0 != isspace(*arg_p)) ++arg_p;
  }

  /* If there are any significant characters in the argument, the argument contains a value. */
  if ('\0' != *arg_p) *type |= VALUE;
  return status;
}

/* Return flag indicating whether this character type contains a non-trivial string
   which is *not* a parameter name. */
static char is_significant_value(unsigned char type) {
  return EQUALS == type || (EQUALS | VALUE) == type || VALUE == type;
}

/* \brief Check for invalid parameter assignments, that is constructions of the form par = value where
   par is not a valid parameter name. This assumes any matches to this pattern are invalid, i.e. that valid
   matches for the given command line arguments have been flagged (using the "type" parameter).
   \param argc The number of command line arguments.
   \param argv The command line arguments.
   \param type Previous classifications of command line arguments including screening for *valid* parameter assignments.
*/
static int check_par_assignment(int argc, char ** argv, unsigned char * type) {
  int status = eOK;
  int end_arg = argc;
  int idx = 0;
  char * begin_name = 0;
  char * end_name = 0;

  /* Check arguments. */
  if (0 > argc) {
    status = eInvalidArgument;
  } else if (0 < argc && (0 == argv || 0 == type)) {
    status = eNullPointer;
  } else if (0 == argc) {
    return status;
  }

  /* Iterate over arguments. */
  for (idx = 0; eOK == status && idx != end_arg; ++idx) {
    char * cp = 0;

    /* If this argument was already handled, its "type" is DONE. In that case, the par = value pattern is broken, so
       go on to next argument and start over looking for the pattern. */
    if (DONE == type[idx]) {
      begin_name = 0;
      end_name = 0;
      continue;
    } else if ('=' != *argv[idx]) {
      /* First character is not equals, so start over looking for the pattern. */
      begin_name = 0;
      end_name = 0;
    }

    /* If par = value pattern not yet detected, start looking for a match to that pattern in this argument. */
    if (0 == begin_name) {
      begin_name = argv[idx];
      end_name = begin_name;
    }

    /* Go through this argument, looking for an unescaped equals, signifying the start to par = value pattern. */
    for (cp = argv[idx]; eOK == status && (0 != isalnum(*cp) || '_' == *cp || '-' == *cp || '.' == *cp); ++cp) {
      /* Build up the parameter name for error reporting purposes, excluding trailing white space. */
      if (begin_name == argv[idx] && 0 == isspace(*cp)) end_name = cp;
    }

    /* If the first character encountered that was not part of a valid parameter name was an equals which is
       not followed by a second equals, the par = value pattern was matched, and therefore we have an invalid
       parameter name. */
    if ('=' == *cp && '=' != *(cp + 1)) status = eInvalidParName;

    /* Final update to end_name to ensure the whole parameter name is reported. This only matters if the par = value
       pattern extends to more than one argument. */
    if (begin_name == argv[idx] && 0 == isspace(*end_name) && '=' != *end_name) ++end_name;

    if (eInvalidParName == status) {
      if (begin_name == end_name) {
        ape_msg_error("Empty parameter name in argument: \"%s\"\n", argv[idx]);
      } else {
        char * name = 0;
        /* Copy the parameter name for reporting purposes. */
        int local_status = ape_util_copy_range(begin_name, end_name, &name);
        if (eOK == local_status) {
          ape_msg_error("Invalid parameter name in command line assignment: \"%s\"\n", name);
        } else {
          ape_msg_error("Invalid parameter name in command line argument: \"%s\"\n", begin_name);
        }

        /* Clean up. */
        free(name); name = 0;
      }
    }
  }

  return status;
}

int ape_io_apply_command_line(ApeParFile * par_file, int argc, char ** argv) {
  int status = eOK;
  int end_par = 0;
  ApePar ** par = 0;
  char ** name = 0;
  char ** mode = 0;
  int idx = 0;
  
  /* Check arguments. */
  if (0 == par_file || 0 == par_file->par_cont) {
    status = eNullPointer;
  } else if (0 > argc) {
    status = eInvalidArgument;
  } else if (0 < argc && 0 == argv) {
    status = eNullPointer;
  } else if (0 == argc) {
    return status;
  }

  if (eOK == status) {
    /* Note this breaks if there are more than 2 billion parameters. */
    end_par = (int) ape_list_get_size(par_file->par_cont);
  }

  if (eOK == status) {
    /* Make an array holding each par from the list. Throw in one extra to make sure the array is not size 0. */
    par = (ApePar **) calloc(end_par + 1, sizeof(ApePar *));
    if (0 != par) {
      ApeListIterator itor = ape_list_begin(par_file->par_cont);
      /* Populate the array with list elements. */
      for (idx = 0; idx != end_par; ++idx) {
        par[idx] = (ApePar *) ape_list_get(itor);
        itor = ape_list_next(itor);
      }
    } else {
      status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* Make an array holding name of each par from the list. Throw in one extra to make sure the array is not size 0. */
    name = (char **) calloc(end_par + 1, sizeof(char *));
    if (0 != name) {
      ApeListIterator itor = ape_list_begin(par_file->par_cont);
      /* Populate the array with names of list elements. */
      for (idx = 0; eOK == status && idx != end_par; ++idx) {
        ApePar * par = (ApePar *) ape_list_get(itor);
        status = ape_par_get_field(par, eName, name + idx);
        itor = ape_list_next(itor);
      }
    } else {
      status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    /* Make an array holding mode of each par from the list. Throw in one extra to make sure the array is not size 0. */
    mode = (char **) calloc(end_par + 1, sizeof(char *));
    if (0 != mode) {
      ApeListIterator itor = ape_list_begin(par_file->par_cont);
      /* Populate the array with mode of list elements. */
      for (idx = 0; idx != end_par; ++idx) {
        mode[idx] = calloc(APE_PAR_MODE_CODE_LEN, sizeof(char));
        if (0 != mode[idx]) {
          ApePar * par = (ApePar *) ape_list_get(itor);
          ape_par_get_mode(par, mode[idx]);
        } else {
          status = eDynAllocFailed;
        }
        itor = ape_list_next(itor);
      }
    } else {
      status = eDynAllocFailed;
    }
  }

  if (eOK == status) {
    int end_arg = argc;

    /* Classifications of each argument. Include room for two extra arguments, because explicit parameters
       can involve an argument plus up to two arguments after it. This is initialized to 0 == UNKNOWN type. */
    unsigned char * arg_type = (unsigned char *) calloc(end_arg + 2, sizeof(unsigned char));
    /* Argument number indexed by parameter number. These must be initialized! */
    int * arg_par = (int *) calloc(end_par, sizeof(int));
    /* Parameter number indexed by argument number. These must be initialized! */
    int * par_arg = (int *) calloc(end_arg, sizeof(int));
    /* Pointers to the values. These are initialized to 0. */
    char ** value = (char **) calloc(end_arg, sizeof(char *));

    if (0 == arg_type || 0 == arg_par || 0 == par_arg || 0 == value) {
      free(value);
      free(par_arg);
      free(arg_par);
      free(arg_type);
      status = eDynAllocFailed;
    } else {
      int arg_index = 0;
      int par_index = 0;

      /* Initialize arguments indexed by parameter number and parameters indexed by argument number to be
         one-past-last of each list. (Default is no matching has occurred. */
      for (arg_index = 0; arg_index < end_arg; ++arg_index) par_arg[arg_index] = end_par;
      for (par_index = 0; par_index < end_par; ++par_index) arg_par[par_index] = end_arg;

      /* Where possible, match up sequences of arguments to form expressions of the form
         PAR EQUALS VALUE. Any argument which fails to form part of such an expression is re-classified
         as a pure value. Go backwards, because that way a re-classified value may still be assigned to
         a valid expression from earlier in the argument list. */
      for (arg_index = end_arg - 1; arg_index >= 0 && eOK == status; --arg_index) {
        /* Classify argument into a form like PAR=VALUE, PAR=, =VALUE, etc. based only on the contents
           of the argument and the set of parameters in this file. */
        status = classify_arg(name, end_par, arg_index, argv[arg_index], arg_type + arg_index, par_arg + arg_index, value + arg_index);

        /* If an argument appears to involve a parameter name, further checking is needed to make sure
           that makes sense in context. */
        if (0 != (PAR & arg_type[arg_index])) {
          if ((PAR | PLUS) == arg_type[arg_index]) {
            value[arg_index] = "yes";
          } else if ((PAR | MINUS) == arg_type[arg_index]) {
            value[arg_index] = "no";
          } else if ((PAR | EQUALS | VALUE) == arg_type[arg_index] || ((PAR | EQUALS) == arg_type[arg_index] &&
            0 == is_significant_value(arg_type[arg_index + 1]))) {
            /* PAR=VALUE. Adjust value to point to the first significant character after EQUALS. */
            ++value[arg_index];
            while (0 != isspace(*value[arg_index])) ++value[arg_index];
          } else if ((PAR | EQUALS) == arg_type[arg_index] && 0 != is_significant_value(arg_type[arg_index + 1])) {
            /* PAR= VALUE or PAR= =VALUE or PAR= =. */
            arg_type[arg_index + 1] = DONE;
            value[arg_index] = value[arg_index + 1];
          } else if (PAR == arg_type[arg_index] && (EQUALS | VALUE) == arg_type[arg_index + 1]) {
            /* PAR =VALUE. */
            ++value[arg_index + 1];
            while (0 != isspace(*value[arg_index + 1])) ++value[arg_index + 1];
            arg_type[arg_index + 1] = DONE;
            value[arg_index] = value[arg_index + 1];
          } else if (PAR == arg_type[arg_index] && EQUALS == arg_type[arg_index + 1] &&
            0 != is_significant_value(arg_type[arg_index + 2])) {
            /* PAR = VALUE or PAR = =VALUE or PAR = =. */
            arg_type[arg_index + 1] = DONE;
            arg_type[arg_index + 2] = DONE;
            value[arg_index] = value[arg_index + 2];
          } else if ((PAR | EQUALS) == arg_type[arg_index]) {
            /* PAR= (blank). */
            ++value[arg_index];
            while (0 != isspace(*value[arg_index])) ++value[arg_index];
          } else if (PAR == arg_type[arg_index] && EQUALS == arg_type[arg_index + 1]) {
            /* PAR = (blank). */
            ++value[arg_index + 1];
            while (0 != isspace(*value[arg_index + 1])) ++value[arg_index + 1];
            arg_type[arg_index + 1] = DONE;
            value[arg_index] = value[arg_index + 1];
          } else if (0 != (PAR & arg_type[arg_index])) {
            /* None of the possible explicit cases were matched, so re-classify this parameter as a value. */
            arg_type[arg_index] = VALUE;
            arg_par[par_arg[arg_index]] = end_arg;
            par_arg[arg_index] = end_par;
            value[arg_index] = argv[arg_index];
            /* Ambiguity isn't an issue for a value. */
            if (eAmbiguousParName == status) status = eOK;
            continue;
          } else {
            continue;
          }
          /* Ambiguity *is* an issue for parameters. */
          if (0 != (PAR & arg_type[arg_index]) && eAmbiguousParName == status) {
	    ape_msg_error("Argument \"%s\" contains an ambiguous parameter name.\n",argv[arg_index]);
            continue;
	  }
          arg_type[arg_index] = DONE;
          arg_par[par_arg[arg_index]] = arg_index;
        }
      }

      /* Check for invalid explicit parameter assignments. */
      if (eOK == status) status = check_par_assignment(argc, argv, arg_type);

      for (arg_index = 0; eOK == status && arg_index < end_arg; ++arg_index) {
        /* Check if parameter was already assigned explicitly; if so, one of the other par_arg's which came before
           will equal the current par_arg. */
        for (idx = 0; idx < arg_index; ++idx) {
          if (par_arg[arg_index] == par_arg[idx] && end_par != par_arg[arg_index]) {
            ape_msg_error("Parameter \"%s\" duplicated on command line.\n", name[par_arg[arg_index]]);
            status = eParameterDuplicated;
            break;
          }
        }
      }

      /* Interpret any remaining arguments as positional parameters. Iterate over all arguments and parameters. */
      for (arg_index = 0, par_index = 0; eOK == status && arg_index != end_arg && par_index != end_par; ) {
        if (DONE == arg_type[arg_index]) {
          /* This argument was already consumed above as part of an explicit assignment. */
          ++arg_index;
        } else if (end_arg != arg_par[par_index]) {
          /* This parameter was already assigned above. */
          if ('a' == *mode[par_index] || 'q' == *mode[par_index]) {
            /* Non-hidden parameter assigned positionally as well. */
            status = eParameterDuplicated;
            ape_msg_error("Parameter \"%s\" duplicated on command line.\n", name[par_index]);
          } else {
            /* Skip this parameter. */
            ++par_index;
          }
        } else {
          if ('a' == *mode[par_index] || 'q' == *mode[par_index]) {
            /* Non-hidden parameter: assign positionally. */
            arg_type[arg_index] = DONE;
            arg_par[par_index] = arg_index;
            par_arg[arg_index] = par_index;
            value[arg_index] = argv[arg_index];
            ++arg_index;
          }
          ++par_index;
        }
      }

      if (eOK == status && arg_index != end_arg) {
        status = eTooManyArguments;
      }

      if (eOK == status) {
        /* Match up positional arguments with left-over (unassigned) parameters. */
        /* Find first eligible parameter by iterating over all parameters and skipping those which were not
           already assigned, or which are hidden. */
        for (par_index = 0; par_index != end_par &&
          (end_arg != arg_par[par_index] || ('a' != *mode[par_index] && 'q' != *mode[par_index])); ++par_index) {}

        /* Iterate over arguments, and assign each one which was not already assigned to a positional parameter. */
        for (arg_index = 0; arg_index < end_arg && par_index < end_par; ++arg_index) {
          if (DONE != arg_type[arg_index]) {
            arg_type[arg_index] = DONE;
            arg_par[par_index] = arg_index;
            par_arg[arg_index] = par_index;
            value[arg_index] = argv[arg_index];
            /* Find next unassigned parameter. */
            for (; par_index != end_par &&
              (end_arg != arg_par[par_index] || ('a' != *mode[par_index] && 'q' != *mode[par_index])); ++par_index) {}
          }
        }
      }

      /* Interpretation is done. Perform assignments of arguments to parameters. */
      for (arg_index = 0; eOK == status && arg_index < end_arg; ++arg_index) {
        /* If this argument should be assigned to a parameter, but was not already assigned, assign it. */
        if (end_par != par_arg[arg_index]) {
          /* int local_status = PIL_OK;
          PIL_free(fp->pp[par_arg[arg_index]].strvalue);
          local_status = PIL_dup(&fp->pp[par_arg[arg_index]].strvalue, value[arg_index]);
          if (PIL_OK != local_status) {
            status = PIL_OK == status ? local_status : status;
          } */
          int local_status = ape_par_set_value_string(par[par_arg[arg_index]], value[arg_index]);
          status = eOK == status ? local_status : status;

          /* Flag this value as being set on the command line. */
          local_status = ape_par_flag_cmd_line(par[par_arg[arg_index]]);
          status = eOK == status ? local_status : status;
        }
      }

      if (eTooManyArguments == status)
        ape_msg_error("Invalid command line: too many arguments.\n");

      /* Clean up. */
      free(value);
      free(par_arg);
      free(arg_par);
      free(arg_type);
    }
  }
  ape_util_free_string_array(mode); mode = 0;
  ape_util_free_string_array(name); name = 0;
  free(par); par = 0;
  return status;
}

int ape_io_check_file_format(ApeParFile * par_file, char check_value) {
  int status = eOK;

  /* Check argument. */
  if (0 == par_file || 0 == par_file->par_cont) {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Extract the parameter from the file. */
    ApeList * par_cont = par_file->par_cont;
    char auto_mode = 0; /* Flag indicating one or more parameters used "auto" mode. */

    ApeListIterator itor = ape_list_begin(par_cont);
    ApeListIterator end = ape_list_end(par_cont);
    for (; itor != end; itor = ape_list_next(itor)) {
      ApePar * par = (ApePar *) ape_list_get(itor);

      /* Check each parameter individually. */
      int local_status = ape_par_check(par, check_value);

      /* Ignore complaints that a parameter is or isn't enumerated. */
      if (eRangeEnum == local_status || eRangeNoEnum == local_status) local_status = eOK;

      /* Status of the file as a whole is the status of the first parameter error. */
      status = eOK == status ? local_status : status;

      /* Check whether this parameter used "auto" mode. */
      if (0 == auto_mode) {
        char mode[APE_PAR_MODE_CODE_LEN] = "";
        local_status = ape_par_get_mode(par, mode);
        if (eOK == local_status && 'a' == *mode) auto_mode = 1;
      }
    }

    for (itor = ape_list_begin(par_cont); itor != end; itor = ape_list_next(itor)) {
      ApeListIterator other_itor = ape_list_next(itor);
      char * name = 0;
      /* Get name of the current parameter. */
      int local_status = ape_par_get_field((ApePar *) ape_list_get(itor), eName, &name);

      /* Proceed if there was no error and the parameter name is not blank. */
      if (eOK == local_status && 0 == is_blank(name)) {
        for (; other_itor != end; other_itor = ape_list_next(other_itor)) {
          char * other_name = 0;

          /* Get name of the other parameter. */
          local_status = ape_par_get_field((ApePar *) ape_list_get(other_itor), eName, &other_name);

          /* If there was no error and the parameter names match, case insensitively, flag error. */
          if (eOK == local_status && 0 == ape_util_strcmp(name, other_name, 1)) {
            ape_msg_debug("Parameter \"%s\" is repeated in parameter file", name);
            if (0 == is_blank(par_file->file_name)) ape_msg_debug(" \"%s\".\n", par_file->file_name);
            else ape_msg_debug(".\n");
            local_status = eParameterDuplicated;
          }
          status = eOK == status ? local_status : status;
          free(other_name); other_name = 0;
        }
      }
      status = eOK == status ? local_status : status;
      free(name); name = 0;
    }

    /* If the "auto" mode was used, check for the default mode. */
    if (0 != auto_mode) {
      char mode[APE_PAR_MODE_CODE_LEN] = "";
      int local_status = ape_io_get_default_mode(par_file, mode);
      if (eParNotFound == local_status) local_status = eNoModePar;
      if (eOK != local_status)
        ape_msg_error("One or more parameters has mode \"a\", but cannot get \"mode\" parameter.\n");
      status = eOK == status ? local_status : status;
    }
  }

  return status;
}

static void test_pfiles(const char * pfiles, const char ** loc_path, size_t num_loc, const char ** sys_path, size_t num_sys) {
  int status = eOK;
  ApeList * list = 0;
  ApeListIterator itor = 0;
  size_t idx = 0;

  /* Present a banner for the test to give a consistent context. */
  if (0 != pfiles)
    ape_msg_info(0, "Testing ape_io pfile/path functions for pfiles == \"%s\"\n", pfiles);
  else
    ape_msg_info(0, "Testing ape_io pfile/path functions for pfiles == 0.\n");

  /* Test ape_io_set_pfiles. */
  status = ape_io_set_pfiles(pfiles);
  ape_test_cmp_string("ape_io_set_pfiles", 0, 0, status, eOK);

  /* Test ape_io_get_loc_path. */
  status = ape_io_get_loc_path(&list);
  ape_test_cmp_string("ape_io_get_loc_path", 0, 0, status, eOK);

  /* Compare results of local path. */
  for (itor = ape_list_begin(list), idx = 0; itor != ape_list_end(list); itor = ape_list_next(itor), ++idx) {
    /* Make sure we don't dereference beyond the end of loc_path. */
    const char * path = idx < num_loc ? loc_path[idx] : 0;
    /* Already reported status above, so short-circuit it here so we don't have to see the message repeated. */
    ape_test_cmp_string("ape_io_get_loc_path", (const char *) ape_list_get(itor), path, eOK, eOK);
  }

  /* Confirm that all directories in path were checked. */
  ape_test_cmp_ulong("number of paths in ape_io_get_loc_path", (unsigned long) idx, (unsigned long) num_loc, eOK, eOK);

  /* Test ape_io_get_sys_path. */
  status = ape_io_get_sys_path(&list);
  ape_test_cmp_string("ape_io_get_sys_path", 0, 0, status, eOK);

  /* Compare results of system path. */
  for (itor = ape_list_begin(list), idx = 0; itor != ape_list_end(list); itor = ape_list_next(itor), ++idx) {
    /* Make sure we don't dereference beyond the end of sys_path. */
    const char * path = idx < num_sys ? sys_path[idx] : 0;
    /* Already reported status above, so short-circuit it here so we don't have to see the message repeated. */
    ape_test_cmp_string("ape_io_get_sys_path", (const char *) ape_list_get(itor), path, eOK, eOK);
  }

  /* Confirm that all directories in path were checked. */
  ape_test_cmp_ulong("number of paths in ape_io_get_sys_path", (unsigned long) idx, (unsigned long) num_sys, eOK, eOK);
}

/* Test ape_io_read_file_path:
   hint Contextual hint about this particular test.
   file_name is the base name of the parameter file to seek.
   expected_found is the full path to the found file.
   path is the path on which to search.
   expected_status is the expected status that the function is expected to return.
*/
static void test_read_path(const char * hint, const char * file_name, const char * expected_found, const char ** path,
  int expected_status) {
  int status = eOK;
  ApeParFile * par_file = 0;
  const char * found_file = 0;
  ApeList * path_list = 0;

  /* Form input path list from given path strings. */
  if (0 != path) {
    const char ** dir = 0;
    path_list = ape_list_create();
    for (dir = path; 0 != *dir; ++dir) {
      ape_list_append(path_list, (void *) *dir);
    }
  }

  /* Read file in from the given path. */
  status = ape_io_read_file_path(file_name, path_list, &par_file);

  /* Clean up path list. */
  ape_list_destroy(path_list); path_list = 0;

  /* See if the name of the file and status returned by read_file_path is as expected. */
  ape_io_get_file_name(par_file, &found_file);
  ape_test_cmp_string(hint, found_file, expected_found, status, expected_status);

  /* Clean up par file structure. */
  if (0 != par_file) {
    destroy_par_list(par_file->par_cont); par_file->par_cont = 0;

    /* Clean up file name. */
    free(par_file->file_name); par_file->file_name = 0;
  }
  free(par_file);
}

/* Sanity check for testing apply_operation, by requiring identical agreement from the two parameter containers. */
static int test_apply_op(ApeList * master_cont, ApeListIterator * master_itor, ApeList * par_cont, ApeListIterator * par_itor,
  char master_more_recent) {
  int status = eOK;

  /* Extract parameters from the master container and the other container. */
  ApePar * master_par = (ApePar *) ape_list_get(*master_itor);
  ApePar * other_par = (ApePar *) ape_list_get(*par_itor);
  char * master_line = ape_par_get_line(master_par);
  char * other_line = ape_par_get_line(other_par);

  ape_msg_debug("test_apply_op found line \"%s\".\n", master_line);
  ape_test_cmp_ptr("test_apply_op comparing pars", master_par, other_par, 0, 0);

  /* Compare the two lines; they should be identical. */
  status = ape_util_strcmp(other_line, master_line, 0);
  ape_test_cmp_string("test_apply_op", other_line, master_line, status, 0);

  /* Clean up. */
  free(other_line); other_line = 0;
  free(master_line); master_line = 0;

  return status;
}

static int create_test_file(const char * file_name, const char ** par_string, ApeParFile ** out_par_file) {
  int status = eOK;
  ApeParFile * par_file = 0;

  /* Check arguments. */
  if (0 != out_par_file) {
    *out_par_file = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && (0 == file_name || 0 == par_string)) {
    status = eNullPointer;
  }

  /* Allocate par file object. */
  if (eOK == status) {
    par_file = (ApeParFile *) calloc(1, sizeof(ApeParFile));
    if (0 == par_file) {
      status = eDynAllocFailed;
    }
  }

  /* Use now for time stamp. */
  par_file->mod_time = time(0);

  /* Copy the parameter name. */
  if (eOK == status) {
    status = ape_util_copy_string(file_name, &par_file->file_name);
  }

  if (eOK == status) {
    /* Create the container of parameters. */
    par_file->par_cont = ape_list_create();
    if (0 == par_file->par_cont) {
      status = eDynAllocFailed;
    }
  }
  if (eOK == status) {
    ApeList * cont = par_file->par_cont;

    /* Convert the given strings into parameters. */
    for (; eOK == status && 0 != *par_string; ++par_string) {
      ApePar * par = 0;
      status = ape_par_create(&par, *par_string);
      if (eOK == status) {
        if (ape_list_end(cont) == ape_list_append(cont, par)) {
          status = eListError;
        }
      }
    }
  }

  if (eOK == status) {
    /* Return the successfully created file. */
    *out_par_file = par_file;
  } else {
    /* Clean up. */
    ape_io_destroy_file(par_file);
  }

  par_file = 0;

  return status;
}

static int get_par_strings(ApeParFile * par_file, char *** out_par_text) {
  int status = eOK;
  char ** par_text = 0;

  if (0 != out_par_text) {
    *out_par_text = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == par_file) {
    status = eNullPointer;
  }

  if (eOK == status) {
    ApeList * par_cont = par_file->par_cont;
    size_t num_par = ape_list_get_size(par_cont) + 1; /* 1 for terminating null. */
    par_text = (char **) calloc(num_par, sizeof(char *));
    if (0 == par_text) {
      status = eDynAllocFailed;
    }
    if (eOK == status) {
      ApeListIterator itor = ape_list_begin(par_cont);
      ApeListIterator end = ape_list_end(par_cont);
      size_t idx = 0;
      for (; eOK == status && itor != end; itor = ape_list_next(itor), ++idx) {
        ApePar * par = (ApePar *) ape_list_get(itor);
        par_text[idx] = ape_par_get_line(par);
        if (0 == par_text[idx]) status = eDynAllocFailed;
      }
    }
  }
  if (eOK == status) {
    *out_par_text = par_text;
  } else {
    ape_util_free_string_array(par_text);
  }

  return status;
}

static void test_apply_command_line(void) {
  const char * par_name[] = { "par0", "par1", "par2", "par3", "par4", "parf", "parfive" };
  const char * default_par_value[] = { "default-val0", "-1", "default-val2", "default-val3", "no", "no", "no" };
  const char * par_value[sizeof(default_par_value)/sizeof(const char *)] = { 0, 0, 0, 0, 0, 0, 0 };

  /* Test nothing on command line. */
  { char * argv[] = { "ape_test" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(0, argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test single positional parameter. */
  { char * argv[] = { "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test two positional parameters. */
  { char * argv[] = { "val0", "1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter. */
  { char * argv[] = { "par0=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter which is not the first parameter in the file. */
  { char * argv[] = { "par1=1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test two explicit parameters. */
  { char * argv[] = { "par0=val0", "par1=1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test two explicit parameters in reverse order. */
  { char * argv[] = { "par1=1", "par0=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one positional, one explicit parameter. */
  { char * argv[] = { "val0", "par1=1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one positional, one explicit parameter in reverse order. */
  { char * argv[] = { "par1=1", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter with simulated space after the =. */
  { char * argv[] = { "par0=", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter with simulated space before the =. */
  { char * argv[] = { "par0", "=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter with simulated space before and after the =. */
  { char * argv[] = { "par0", "=", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter with simulated space before and after the =, followed by positional parameter. */
  { char * argv[] = { "par1", "=", "1", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one explicit parameter with space before and after the = (double quotes around whole argument). */
  { char * argv[] = { "par0  = \tval0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one positional, one explicit parameter with simulated space after the =. */
  { char * argv[] = { "val0", "par1=", "1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one positional, one explicit parameter with simulated space after the =, positional parameter is not the first
     one in the file. */
  { char * argv[] = { "1", "par0=", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eParameterDuplicated);
  }
  /* Test one positional, one explicit parameter with simulated space before the =, positional parameter is not the first
     one in the file. */
  { char * argv[] = { "par1", "=1", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test one positional, one explicit parameter with simulated space before the =, positional parameter is not the first
     one in the file and reverse order of arguments. */
  { char * argv[] = { "par0", "=val0", "1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eParameterDuplicated);
  }
  /* Test explicit parameters with all blank values. */
  { char * argv[] = { "par0=", "par2=", "par3=" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "";
    par_value[2] = "";
    par_value[3] = "";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test explicit parameters with all values assigned, in variety of flavors, with extra spaces around =. */
  { char * argv[] = { "par1 = ", "1", "par3", " = ", "val3", "par0", " = val0", "par2 = val2" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    par_value[2] = "val2";
    par_value[3] = "val3";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test positional parameters whose values happen to be the same as parameter names. */
  { char * argv[] = { "par0", "par2", "=", "par2", "1", "par3" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par0";
    par_value[1] = "1";
    par_value[2] = "par2";
    par_value[3] = "par3";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test explicit assignment from a value which happens to be the same as other parameter name. */
  { char * argv[] = { "par0", "=", "par1" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par1";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test explicit assignments from a value which happens to be an equals sign. */
  { char * argv[] = { "par0", "=", "=" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "=";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  { char * argv[] = { "par0=", "=" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "=";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  { char * argv[] = { "par0", "==" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "=";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  { char * argv[] = { "par0==" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "=";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test parameters which happen to have equals signs in them. */
  { char * argv[] = { "par2 =", "= val2", "[=]", "1", "par3==val3" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "[=]";
    par_value[1] = "1";
    par_value[2] = "= val2";
    par_value[3] = "=val3";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test positional parameters which have equals signs in them. */
  { char * argv[] = { "[=]" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "[=]";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test explicit assignment with equals signs. */
  { char * argv[] = { "par0", "=", "=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "=val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test mixed command line containing unmatched bracket error with positional parameters which have equals signs in them. */
  { char * argv[] = { "par0", "par2=val2", "=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eInvalidParName);
  }
  /* Test unmatched bracket error with positional parameters which have equals signs in them. */
  { char * argv[] = { "][=" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "][=";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test matched bracket with positional parameters which have equals signs in them. */
  { char * argv[] = { "[[=]][=]" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "[[=]][=]";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test matched bracket with non-trivial positional parameters which have equals signs in them. */
  { char * argv[] = { "somefile.fits[#row=8]" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "somefile.fits[#row=8]";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test parameter duplicated positionally. */
  { char * argv[] = { "0", "par0=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eParameterDuplicated);
  }
  /* Test bare = sign. */
  { char * argv[] = { "=" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eInvalidParName);
  }
  /* Test trailing = sign 4. */
  { char * argv[] = { "==" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "==";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test = sign in value. */
  { char * argv[] = { "=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eInvalidParName);
  }
  /* Test correct handling of == . */
  { char * argv[] = { "par0==val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "=val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test correct handling of == . */
  { char * argv[] = { "pr0==val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "pr0==val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test correct handling of == . */
  { char * argv[] = { "pr0=val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eInvalidParName);
  }
  /* Test a common typo. */
  { char * argv[] = { "val0", "1", "pr2=val2" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eInvalidParName);
  }
  /* Test another common typo. */
  { char * argv[] = { "val0", "1", "val3", "pr2=val2" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eInvalidParName);
  }
  /* Test a typo which can't be detected by ape. */
  { char * argv[] = { "val0", "1", "pr2<=val2" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    par_value[3] = "pr2<=val2";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eOK);
  }
  /* Test another common mistake: parameter duplicated. */
  { char * argv[] = { "par0=val0A", "par1=1", "par0=", "val0B" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eParameterDuplicated);
  }
  /* Test two explicit parameters, simulating different amounts of whitespace. */
  { char * argv[] = { "val0", "par1", "=", "1", "par2=", "val2" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[1] = "1";
    par_value[2] = "val2";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test too many arguments on command line. */
  { char * argv[] = { "val0", "val1", "val3", "no", "val5" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eTooManyArguments);
  }
  /* Test boolean parameter with + by itself. */
  { char * argv[] = { "par4+" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[4] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test boolean parameter with + and a positional parameter. */
  { char * argv[] = { "par4+", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "val0";
    par_value[4] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test boolean parameter with + in the wrong place and a positional parameter. */
  { char * argv[] = { "par4", "+", "val0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par4";
    par_value[1] = "+";
    par_value[3] = "val0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test boolean parameter with +, as well as a near-miss with embedded + in value. */
  { char * argv[] = { "par4+", "par3+val3" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par3+val3";
    par_value[4] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test two parameters with +/-, and some other variations. */
  { char * argv[] = { "par4+", "par3-", "par2++", "+" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par2++";
    par_value[1] = "+";
    par_value[3] = "no";
    par_value[4] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test boolean parameter with + along with many others. */
  { char * argv[] = { "par4+", "par3", "+2", "par2++" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par3";
    par_value[1] = "+2";
    par_value[3] = "par2++";
    par_value[4] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test partial matching. */
  { char * argv[] = { "parfi=yes" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[6] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test potential - but not actual - ambiguity. */
  { char * argv[] = { "parf=yes" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[5] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  /* Test ambiguity. */
  { char * argv[] = { "par=yes" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eAmbiguousParName);
  }
  { char * argv[] = { "par = yes" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eAmbiguousParName);
  }
  { char * argv[] = { "par", "=", "yes" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value,
      eAmbiguousParName);
  }
  { char * argv[] = { "par" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  { char * argv[] = { "par", "2", "par", "yes" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par";
    par_value[1] = "2";
    par_value[3] = "par";
    par_value[4] = "yes";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
  { char * argv[] = { "par0" };
    memcpy((void *) par_value, default_par_value, sizeof(default_par_value));
    par_value[0] = "par0";
    test_one_command_line(sizeof(argv) / sizeof(char *), argv, sizeof(par_name) / sizeof(char *), par_name, par_value, eOK);
  }
}

static void test_one_command_line(int argc, char ** argv, int num_par, const char ** par_name, const char ** par_value,
  int expected_status) {
  int status = eOK;
  char line[APE_BUF_SIZE] = "";
  const char * par_file_text[] = {
    "par0    , s, a, \"default-val0\", , , \"Enter par0\"",
    "par1    , r, a, \"-1\", , , \"Enter par1\"",
    "par2    , s, h, \"default-val2\", , , \"Enter par2\"",
    "par3    , s, a, \"default-val3\", , , \"Enter par3\"",
    "par4    , b, a, no, , , \"Enter par4\"",
    "parf    , b, h, no, , , \"Enter parf\"",
    "parfive , b, h, no, , , \"Enter parfive\"",
    "mode    , s, h, ql, , , \"Mode\"",
    0
  };
  ApeParFile * par_file = 0;
  status = create_test_file("", par_file_text, &par_file);
  if (eOK == status) {
    int ii = 0;

    cat_cmd_line(argc, argv, line);
    ape_msg_debug("test_one_command_line: testing line |%s|\n", line);

    status = ape_io_apply_command_line(par_file, argc, argv);
    if (expected_status != status) {
      ape_test_failed("test_one_command_line: ape_io_apply_command_line returned status %d, not %d as expected.\n", status,
        expected_status);
    }
    for (ii = 0; ii < num_par; ++ii) {
      ApePar * par = 0;
      ApeListIterator itor = 0;
      char * value = 0;
      status = ape_io_find_par(par_name[ii], par_file, &itor);
      if (eOK == status) {
        par = (ApePar *) ape_list_get(itor);
        status = ape_par_get_field(par, eValue, &value);
        if (eOK != status) {
          ape_test_failed("test_one_command_line: ape_par_get_field failed for parameter \"%s\".\n", par_name[ii]);
          free(value); value = 0;
          break;
        }
      } else {
        ape_test_failed("test_one_command_line: ape_io_find_par unexpectedly failed to find parameter \"%s\".\n", par_name[ii]);
      }
      if (eOK == status) {
        if (0 != strcmp(value, par_value[ii])) {
          ape_test_failed("test_one_command_line: after applying command line arguments, parameter \"%s\" is \"%s\" not \"%s\" "
            "as expected.\n", par_name[ii], value, par_value[ii]);
        }
      }
      free(value); value = 0;
    }
  } else {
    ape_test_failed("ape_io_test could not set up to test ape_io_apply_command_line (status is %d).\n", status);
  }
  ape_io_destroy_file(par_file); par_file = 0;

}

static void cat_cmd_line(int argc, char ** argv, char * line) {
  if (0 >= argc) {
    *line = '\0';
  } else {
    int ii;

    strcpy(line, "\"");
    strcat(line, argv[0]);
    strcat(line, "\"");

    for (ii = 1; ii < argc; ++ii) {
      strcat(line, " \"");
      strcat(line, argv[ii]);
      strcat(line, "\"");
    }
  }
}

static void test_revert_unlearned(const char * msg, const char ** current_string, const char ** previous_string,
  const char ** expected) {
  int status = eOK;
  ApeParFile * current = 0;
  ApeParFile * previous = 0;
  
  status = create_test_file("", current_string, &current);
  if (eOK == status) {
    status = create_test_file("", previous_string, &previous);
  }
  if (eOK == status) {
    status = ape_io_revert_unlearned(current, previous);
    if (eOK == status) {
      ApeList * par_cont = current->par_cont;
      ApeListIterator itor = 0;
      for (itor = ape_list_begin(par_cont); itor != ape_list_end(par_cont) && 0 != *expected; itor = ape_list_next(itor)) {
        ApePar * par = (ApePar *) ape_list_get(itor);
        char * par_line = ape_par_get_line(par);
        ape_test_cmp_string(msg, par_line, *expected, eOK, eOK);
        free(par_line); par_line = 0;
        ++expected;
      }
      if (ape_list_end(par_cont) != itor) ape_test_failed("%s current has more parameters than expected.\n", msg);
      if (0 != *expected) ape_test_failed("%s current has fewer parameters than expected.\n", msg);
    } else {
      ape_test_failed("%s returned status %d, not %d as expected.\n", msg, status, eOK);
    }
  } else {
    ape_test_failed("ape_io_test could not set up to test ape_io_revert_unlearned (status is %d).\n", status);
  }
  ape_io_destroy_file(previous); previous = 0;
  ape_io_destroy_file(current); current = 0;
}

void ape_io_test(void) {
  ApeParFile * par_file = 0;
  const char * file_name = 0;
  int status = eOK;

  ape_msg_debug("ape_io_test: beginning tests.\n");

  /* Attempt to read a non-existent parameter file. */
  ape_msg_debug("ape_io_test: verify attempting to read a non-existent parameter file fails as expected.\n");
  status = ape_io_read("non-existent.par", &par_file);
  if (eFileReadError != status ) {
    ape_test_failed("ape_io_read(\"non-existent.par\") returned status %d, not %d as expected.\n", status,
      eFileReadError);
  } else if (0 != par_file) {
    ape_test_failed("ape_io_read(\"non-existent.par\") returned non-0 par file pointer.\n");
  }
  ape_msg_debug("ape_io_test:\n");

  /* Fabricate a test file containing as many legal variations as practical. */
  { const char * par_string[] = {
      "# Include all distinct modes.",
      "sa, s, a, , , , string automatic",
      "sal, s, al, , , , string automatic learned",
      "sh, s, h, , , , string hidden",
      "shl, s, hl, , , , string hidden learned",
      "sq, s, q, , , , string queried",
      "sql, s, ql, , , , string queried learned",
      "",
      "# White space variations:",
      "# Empty, blank",
      "",
      " ",
      "# Comment variations:",
      "# Empty comment, blank comment.",
      "#",
      "# ",
      "# Empty comment, blank comment with leading white space.",
      "	#",
      "	#	",
      "	# Comment with leading white space.",
      "# Comment at end of a line.",
      "sh_comment, s, h, , , , string hidden with comment # this is not part of the prompt",
      "#",
      "# All distinct parameter types without value, min or max.",
      "bh, b, h, , , , boolean hidden",
      "dh, d, h, , , , double hidden",
      "fh, f, h, , , , file name hidden",
      "frh, fr, h, , , , file name readable hidden",
      "fwh, fw, h, , , , file name writable hidden",
      "ih, i, h, , , , int hidden",
      "rh, r, h, , , , real hidden",
      "#sh, s, h, , , , string hidden # already appeared above.",
      "#",
      "# All distinct parameter types without value or max.",
      "bhmin, b, h, , no, , boolean hidden",
      "dhmin, d, h, , 0., , double hidden",
      "fhmin, f, h, , M, , file name hidden",
      "frhmin, fr, h, , M, , file name readable hidden",
      "fwhmin, fw, h, , M, , file name writable hidden",
      "ihmin, i, h, , 0, , int hidden",
      "rhmin, r, h, , 0., , real hidden",
      "shmin, s, h, , M, , string hidden",
      "#",
      "# All distinct parameter types with enumerated range.",
      "bhenum, b, h, , no|yes,, boolean hidden",
      "dhenum, d, h, , 0.|1.|2., , double hidden",
      "fhenum, f, h, , ape_test.par|., , file name hidden",
      "frhenum, fr, h, , ape_test.par|., , file name readable hidden",
      "fwhenum, fw, h, , ape_test.par|., , file name writable hidden",
      "ihenum, i, h, , 0|1|2|3, , int hidden",
      "rhenum, r, h, , 0.|1.|2.|3., , real hidden",
      "shenum, s, h, , S0|S1|S2, , string hidden",
      "#",
      "# All distinct parameter types without value or min.",
      "bhmax, b, h, , , yes, boolean hidden",
      "dhmax, d, h, , , 10., double hidden",
      "fhmax, f, h, , , m, file name hidden",
      "frhmax, fr, h, , , m, file name readable hidden",
      "fwhmax, fw, h, , , m, file name writable hidden",
      "ihmax, i, h, , , 10, int hidden",
      "rhmax, r, h, , , 10., real hidden",
      "shmax, s, h, , , m, string hidden",
      "#",
      "# All distinct parameter types without value.",
      "bhrange, b, h, , no, yes, boolean hidden",
      "dhrange, d, h, , 0., 10., double hidden",
      "fhrange, f, h, , M, m, file name hidden",
      "frhrange, fr, h, , M, m, file name readable hidden",
      "fwhrange, fw, h, , M, m, file name writable hidden",
      "ihrange, i, h, , 0, 10, int hidden",
      "rhrange, r, h, , 0., 10., real hidden",
      "shrange, s, h, , M, m, string hidden",
      "#",
      "# All distinct parameter types with OK values.",
      "bhvalid, b, h, no, no, yes, boolean hidden",
      "dhvalid, d, h, 5., 0., 10., double hidden",
      "fhvalid, f, h, Z, M, m, file name hidden",
      "frhvalid, fr, h, Z, M, m, file name readable hidden",
      "fwhvalid, fw, h, Z, M, m, file name writable hidden",
      "ihvalid, i, h, 5, 0, 10, int hidden",
      "rhvalid, r, h, 5., 0., 10., real hidden",
      "shvalid, s, h, Z, M, m, string hidden",
      "#",
      "# All distinct parameter types with values out of range low.",
      "bhlow, b, h, maybe, no, yes, boolean hidden",
      "dhlow, d, h, -1., 0., 10., double hidden",
      "fhlow, f, h, A, M, m, file name hidden",
      "frhlow, fr, h, A, M, m, file name readable hidden",
      "fwhlow, fw, h, A, M, m, file name writable hidden",
      "ihlow, i, h, -1, 0, 10, int hidden",
      "rhlow, r, h, -1., 0., 10., real hidden",
      "shlow, s, h, A, M, m, string hidden",
      "#",
      "# All distinct parameter types with values out of range high.",
      "bhhigh, b, h, zaybe, no, yes, boolean hidden",
      "dhhigh, d, h, 11., 0., 10., double hidden",
      "fhhigh, f, h, z, M, m, file name hidden",
      "frhhigh, fr, h, z, M, m, file name readable hidden",
      "fwhhigh, fw, h, z, M, m, file name writable hidden",
      "ihhigh, i, h, 11, 0, 10, int hidden",
      "rhhigh, r, h, 11., 0., 10., real hidden",
      "shhigh, s, h, z, M, m, string hidden",
      "#",
      "# Numerical parameter types with invalid values.",
      "dhinvalid, d, h, five, 0., 10., double hidden",
      "ihinvalid, i, h, 5., 0, 10, int hidden",
      "rhinvalid, r, h, 5. s, 0., 10., real hidden",
      "#",
      "# Numerical parameter types with special values.",
      "ihindef, i, h, inDef, , , int hidden indef",
      "ihinf, i, h, Inf, , , int hidden inf",
      "ihinfinity, i, h, infiNity, , , int hidden infinity",
      "ihnan, i, h, NaN, , , int hidden nan",
      "ihnone, i, h, none, , , int hidden none",
      "ihundef, i, h, unDef, , , int hidden undef",
      "ihundefined, i, h, unDefined, , , int hidden undefined",
      "rhindef, r, h, inDef, , , real hidden indef",
      "rhinf, r, h, Inf, , , real hidden inf",
      "rhinfinity, r, h, infiNity, , , real hidden infinity",
      "rhnan, r, h, NaN, , , real hidden nan",
      "rhnone, r, h, none, , , real hidden none",
      "rhundef, r, h, unDef, , , real hidden undef",
      "rhundefined, r, h, unDefined, , , real hidden undefined",
      "#",
      "mode, s, h, ql, , , \"Mode for automatic parameters\"",
      0
    };
    ApeParFile * read_file = 0;

    file_name = "ape_test.par";
    ape_msg_debug("ape_io_test: in preparation for tests, creating test parameter file object.\n");

    status = create_test_file(file_name, par_string, &par_file);
    if (eOK != status) {
      ape_test_failed("ape_io_test: failed to set up ApeParFile test object.\n");
      return;
    }

    ape_msg_debug("ape_io_test:\n");

    /* Test writing these parameters to output file. */
    file_name = par_file->file_name;
    ape_msg_debug("ape_io_test: testing ape_io_write by writing file \"%s\".\n", file_name);
    status = ape_io_write(par_file, 0);
    if (eOK != status) {
      ape_test_failed("ape_io_test: ape_io_write for file \"%s\" returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test:\n");

    /* Test reading these parameters from test unix input file. */
    file_name = "ape_test_unix.par";
    ape_msg_debug("ape_io_test: testing ape_io_read by reading file \"%s\".\n", file_name);
    status = ape_io_read(file_name, &read_file);
    if (eOK != status) {
      ape_test_failed("ape_io_test: ape_io_read(\"%s\", &read_file) returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test:\n");

    /* Verify the read parameters by writing them back out. */
    file_name = "ape_test_unix_copy.par";
    free(read_file->file_name); read_file->file_name = 0;
    status = ape_util_copy_string(file_name, &read_file->file_name);
    if (eOK != status) {
      ape_test_failed("ape_io_test: copy_string for file \"%s\" returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test: testing ape_io_write by writing file \"%s\".\n", file_name);
    status = ape_io_write(read_file, 0);
    if (eOK != status) {
      ape_test_failed("ape_io_test: ape_io_write for file \"%s\" returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test:\n");
    ape_io_destroy_file(read_file); read_file = 0;

    /* Test reading these parameters from test windows input file. */
    file_name = "ape_test_win32.par";
    ape_msg_debug("ape_io_test: testing ape_io_read by reading file \"%s\".\n", file_name);
    status = ape_io_read(file_name, &read_file);
    if (eOK != status) {
      ape_test_failed("ape_io_test: ape_io_read(\"%s\", &read_file) returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test:\n");

    /* Verify the read parameters by writing them back out. */
    file_name = "ape_test_win32_copy.par";
    free(read_file->file_name); read_file->file_name = 0;
    status = ape_util_copy_string(file_name, &read_file->file_name);
    if (eOK != status) {
      ape_test_failed("ape_io_test: copy_string for file \"%s\" returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test: testing ape_io_read by writing file \"%s\".\n", file_name);
    status = ape_io_write(read_file, 0);
    if (eOK != status) {
      ape_test_failed("ape_io_test: ape_io_write(read_file) returned status = %d, not %d as expected.\n",
        file_name, status, eOK);
      return;
    }
    ape_msg_debug("ape_io_test:\n");
    ape_io_destroy_file(read_file); read_file = 0;
  }

  /* Clean up parameter file. */
  ape_io_destroy_file(par_file); par_file = 0;

  /* Test ape_io_set_pfiles, ape_io_get_loc_path and ape_io_get_sys_path. */
  { const char * pfiles = 0; /* Assume for this test that PFILES is not set. */
    const char ** loc = 0;
    const char ** sys = 0;
    test_pfiles(pfiles, loc, 0, sys, 0);
  }
  { const char * pfiles = "";
    const char ** loc = 0;
    const char ** sys = 0;
    test_pfiles(pfiles, loc, 0, sys, 0);
  }
  { const char * pfiles = ".";
    const char * loc[] = { "." };
    const char ** sys = 0;
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, 0);
  }
  { const char * pfiles = "."QS;
    const char * loc[] = { "." };
    const char ** sys = 0;
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, 0);
  }
  { const char * pfiles = QS".";
    const char ** loc = 0;
    const char * sys[] = { "." };
    test_pfiles(pfiles, loc, 0, sys, sizeof(sys) / sizeof(const char *));
  }
  { const char * pfiles = "."QS".";
    const char * loc[] = { "." };
    const char * sys[] = { "." };
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, sizeof(sys) / sizeof(const char *));
  }
  { const char * pfiles = "/local/pfiles1"QC"/local/pfiles2"QS"/system/pfiles1"QC"/system/pfiles2";
    const char * loc[] = { "/local/pfiles1", "/local/pfiles2" };
    const char * sys[] = { "/system/pfiles1", "/system/pfiles2" };
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, sizeof(sys) / sizeof(const char *));
  }
  { const char * pfiles = "."QS QC;
    const char * loc[] = { "." };
    const char * sys[] = { "", "" };
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, sizeof(sys) / sizeof(const char *));
  }
  { const char * pfiles = "."QS QC"/system/pfiles1";
    const char * loc[] = { "." };
    const char * sys[] = { "", "/system/pfiles1" };
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, sizeof(sys) / sizeof(const char *));
  }
  { const char * pfiles = "."QS QC"/system/pfiles1"QC QC"/system/pfiles2";
    const char * loc[] = { "." };
    const char * sys[] = { "", "/system/pfiles1", "", "/system/pfiles2" };
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, sizeof(sys) / sizeof(const char *));
  }
  { const char * pfiles = QC"."QC QC"/local/pfiles1"QC " "QC QS QC"/system/pfiles1"QC QC"/system/pfiles2"QC;
    const char * loc[] = { "", ".", "", "/local/pfiles1", " ", "" };
    const char * sys[] = { "", "/system/pfiles1", "", "/system/pfiles2", "" };
    test_pfiles(pfiles, loc, sizeof(loc) / sizeof(const char *), sys, sizeof(sys) / sizeof(const char *));
  }

  /* Test ape_io_file_get_name. */
  { ApeParFile * file = 0;
    const char * expected = 0;
    const char * name = 0;
    status = ape_io_get_file_name(file, &name);
    ape_test_cmp_string("ape_io_get_file_name called for a null ApeParFile pointer", name, expected, status, eOK);
  }
  { ApeParFile file = { "spud" };
    const char * expected = 0;
    const char ** name = 0;
    status = ape_io_get_file_name(&file, name);
    ape_test_cmp_string("ape_io_get_file_name called with a null name pointer", expected, expected,
      status, eNullPointer);
  }
  { ApeParFile file = { "spud" };
    const char * expected = "spud";
    const char * name = 0;
    status = ape_io_get_file_name(&file, &name);
    ape_test_cmp_string("ape_io_get_file_name called for a legit ApeParFile pointer", name, expected, status, eOK);
  }
  { ApeParFile file = { 0 };
    const char * expected = 0;
    const char * name = 0;
    status = ape_io_get_file_name(&file, &name);
    ape_test_cmp_string("ape_io_get_file_name called for an ApeParFile pointer with null name", name, expected,
      status, eOK);
  }

  /* Test ape_io_read_file_path. */
  /* In preparation for this, reset pfiles to its default value. */
  status = ape_io_set_pfiles(0); if (eOK == status) {
    { const char * file_name = "ape_test.par";
      const char * expected_found = 0;
      const char ** path = 0;
      test_read_path("ape_io_read_file_path with null path", file_name, expected_found, path, eFileNotFound);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = 0;
      const char * path[] = { 0 };
      test_read_path("ape_io_read_file_path with empty path", file_name, expected_found, path, eFileNotFound);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = 0;
      const char * path[] = { "", 0 };
      test_read_path("ape_io_read_file_path with empty string in path", file_name, expected_found, path, eFileNotFound);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = "."QD"ape_test.par";
      const char * path[] = { ".", 0 };
      test_read_path("ape_io_read_file_path with path = { . }", file_name, expected_found, path, eOK);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = "."QD"ape_test.par";
      const char * path[] = { "non-existent-dir", ".", 0 };
      test_read_path("ape_io_read_file_path with path = { non-existent-dir, . }", file_name, expected_found,
        path, eOK);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = "."QD"ape_test.par";
      const char * path[] = { ".", "non-existent-dir", 0 };
      test_read_path("ape_io_read_file_path with path = { ., non-existent-dir }", file_name, expected_found,
        path, eOK);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = "."QD"."QD"ape_test.par";
      const char * path[] = { "."QD".", ".", 0 };
      test_read_path("ape_io_read_file_path with path = { ./., . }", file_name, expected_found, path, eOK);
    }
    { const char * file_name = 0;
      const char * expected_found = 0;
      const char * path[] = { ".", 0 };
      test_read_path("ape_io_read_file_path with valid path but null par file name", file_name, expected_found,
        path, eNullPointer);
    }
    { const char * file_name = "";
      const char * expected_found = 0;
      const char * path[] = { ".", 0 };
      test_read_path("ape_io_read_file_path with valid path but empty par file name", file_name, expected_found,
        path, eFileReadError);
    }
    { const char * file_name = "ape_test.par";
      const char * expected_found = "."QD"ape_test.par";
      const char * path[] = { "", " ", ".", 0 };
      test_read_path("ape_io_read_file_path with path containing blank directory names", file_name, expected_found,
        path, eOK);
    }
  }

  /* Do some tests which require a parameter file. */
  { ApeParFile * par_file = 0;
    ApeListIterator par_itor = 0;
    status = ape_io_read("ape_test_unix.par", &par_file);
    if (eOK == status) {
      /* Test ape_io_find_par. */

      /* Test for a parameter which is present, with matching case. */
      const char * expected_par_name = "sql";
      status = ape_io_find_par(expected_par_name, par_file, &par_itor);
      if (eOK == status) {
        char * par_name = 0;
        ApePar * par = (ApePar *) ape_list_get(par_itor);
        status = ape_par_get_field(par, eName, &par_name);
        ape_test_cmp_string("ape_io_find_par(\"sql\", ...)", par_name, expected_par_name, status, eOK);
        free(par_name); par_name = 0;
      } else {
        ape_test_cmp_string("ape_io_find_par(\"sql\", ...) did not find par", 0, 0, status, eOK);
      }

      /* Test for a parameter which is present, with mixed case. */
      status = ape_io_find_par("sqL", par_file, &par_itor);
      if (eOK == status) {
        char * par_name = 0;
        ApePar * par = (ApePar *) ape_list_get(par_itor);
        status = ape_par_get_field(par, eName, &par_name);
        ape_test_cmp_string("ape_io_find_par(\"sqL\", ...)", par_name, expected_par_name, status, eOK);
        free(par_name); par_name = 0;
      } else {
        ape_test_cmp_string("ape_io_find_par(\"sqL\", ...) did not find par", 0, 0, status, eOK);
      }

      /* Test for a parameter which is *not* present. */
      expected_par_name = "non-exist";
      status = ape_io_find_par(expected_par_name, par_file, &par_itor);
      ape_test_cmp_ptr("ape_io_find_par(\"non-exist\", ...)", par_itor, ape_list_end(par_file->par_cont), status, eParNotFound);

      /* Test apply_operation. */
      status = apply_operation(par_file, par_file, &test_apply_op);

    } else {
      ape_test_failed("ape_io_test could not set up to run tests using ape_test_unix.par.\n");
    }
    ape_io_destroy_file(par_file);
  }

  /* Test ape_io_merge_par_files, first for case where the local file values should trump system, except system file is
     more recent, causing all hidden parameters to keep their system values. */
  { const char * local_text[] = {
      "# Common file layout.",
      "hidden, i, h, 100, 0|10|100, , \"prompt\"",
      "automatic, i, a, 100, 0|10|100, , \"prompt\"",
      "mode, s, h, \"ql\", , , \"prompt\"",
      0
    };
    char * system_text[] = {
      "# Common file layout.",
      "hidden, i, h, 10, 0|10|100, , \"prompt\"",
      "automatic, i, a, 10, 0|10|100, , \"prompt\"",
      "mode, s, h, \"ql\", , , \"prompt\"",
      0
    };
    const char ** system_text_cp = (const char **) system_text;

    ApeParFile * local = 0;
    ApeParFile * system = 0;
    file_name = "";
    /* Creating the files in this order means the system file will have a later modification time. */
    status = create_test_file("local", local_text, &local);
    ape_util_sleep(1);
    if (eOK == status) {
      status = create_test_file("system", system_text_cp, &system);
    }
    if (eOK == status) {
      ApeParFile * merged = 0;
      status = ape_io_merge_par_files(system, local, &merged);
      if (eOK == status) {
        const char * expected_text[] = {
          "# Common file layout.",
          "hidden, i, h, 10, 0|10|100, , \"prompt\"",
          "automatic, i, a, 100, 0|10|100, , \"prompt\"",
          "mode, s, h, \"ql\", , , \"prompt\"",
          0
        };
        char ** merged_text = 0;
        status = get_par_strings(merged, &merged_text);
        ape_test_cmp_string_array("ape_io_merge_par_files test 0A", merged_text, expected_text, status, eOK);
        ape_util_free_string_array(merged_text);
      } else {
        ape_test_failed("ape_io_test: ape_io_merge_par_files test 0A returned status %d, not %d as expected.\n", status,
          eOK);
      }
      ape_io_destroy_file(merged); merged = 0;
    } else {
      ape_test_failed("ape_io_test could not set up to test ape_io_merge_par_files (status is %d).\n", status);
    }
    ape_io_destroy_file(system); system = 0;
    ape_io_destroy_file(local); local = 0;
  }

  /* Test ape_io_merge_par_files, first for case where the local file mostly survives. */
  { const char * local_text[] = {
      "Changeenum, i, h, 10, 0|10|100, , \"local\"",
      "Changemax, i, \"h\", 5, 0, 100, \"local\"",
      "# Local file layout.",
      "Changemin, i, h, 5, -10, 10, \"local\"",
      "Changemode, i, a, 5, 0, 10, \"local\"",
      "Changetype, r, h, 5, 0, 10, \"local\"",
      "Changevalue, r, h, 5, 0, 10, \"local\"",
      "Changenothing, r, h, 5, 0, 10, \"local\"",
      "Changepar, s, h, \"unterminated quote, min, max, \"local\"",
      0
    };
    char * system_text[] = {
      "changepar,s,h,,min,max,\"system\"",
      "changenothing,r,h,5,0,10,\"system\"",
      "# System file layout.",
      "changevalue,r,h,0,0,10,\"system\"",
      "changetype,i,h,0,0,10,\"system\"",
      "changemode,i,h,0,0,10,\"system\"",
      "changemin,i,h,0,0,10,\"system\"",
      "changemax,i,h,0,0,10,\"system\"",
      "changeenum,i,h,0,0|10,,\"system\"",
      0, /* Leave room for an additional parameter for a future test. */
      0
    };
    const char ** system_text_cp = (const char **) system_text;

    ApeParFile * local = 0;
    ApeParFile * system = 0;
    file_name = "";
    /* Creating the files in this order means the system file will have a later modification time. */
    status = create_test_file("local", local_text, &local);
    ape_util_sleep(1);
    if (eOK == status) {
      status = create_test_file("system", system_text_cp, &system);
    }
    if (eOK == status) {
      ApeParFile * merged = 0;
      status = ape_io_merge_par_files(system, local, &merged);
      if (eOK == status) {
        const char * expected_text[] = {
          "Changeenum, i, h, 0, 0|10, , \"local\"",
          "Changemax, i, \"h\", 0, 0, 10, \"local\"",
          "# Local file layout.",
          "Changemin, i, h, 0, 0, 10, \"local\"",
          "Changemode, i, h, 0, 0, 10, \"local\"",
          "Changetype, i, h, 0, 0, 10, \"local\"",
          "Changevalue, r, h, 0, 0, 10, \"local\"",
          "Changenothing, r, h, 5, 0, 10, \"local\"",
          "changepar,s,h,,min,max,\"system\"",
          0
        };
        char ** merged_text = 0;
        status = get_par_strings(merged, &merged_text);
        ape_test_cmp_string_array("ape_io_merge_par_files test 0", merged_text, expected_text, status, eOK);
        ape_util_free_string_array(merged_text);
      } else {
        ape_test_failed("ape_io_test: ape_io_merge_par_files test 0 returned status %d, not %d as expected.\n", status,
          eOK);
      }
      ape_io_destroy_file(merged); merged = 0;
    } else {
      ape_test_failed("ape_io_test could not set up to test ape_io_merge_par_files (status is %d).\n", status);
    }
    ape_io_destroy_file(system); system = 0;
    ape_io_destroy_file(local); local = 0;

    /* Creating the files in this order means the local file will have a later modification time. */
    status = create_test_file("system", system_text_cp, &system);
    ape_util_sleep(1);
    if (eOK == status) {
      status = create_test_file("local", local_text, &local);
    }
    if (eOK == status) {
      ApeParFile * merged = 0;
      status = ape_io_merge_par_files(system, local, &merged);
      if (eOK == status) {
        const char * expected_text[] = {
          "Changeenum, i, h, 10, 0|10, , \"local\"",
          "Changemax, i, \"h\", 5, 0, 10, \"local\"",
          "# Local file layout.",
          "Changemin, i, h, 5, 0, 10, \"local\"",
          "Changemode, i, h, 5, 0, 10, \"local\"",
          "Changetype, i, h, 5, 0, 10, \"local\"",
          "Changevalue, r, h, 5, 0, 10, \"local\"",
          "Changenothing, r, h, 5, 0, 10, \"local\"",
          "changepar,s,h,,min,max,\"system\"",
          0
        };
        char ** merged_text = 0;
        status = get_par_strings(merged, &merged_text);
        ape_test_cmp_string_array("ape_io_merge_par_files test 1", merged_text, expected_text, status, eOK);
        ape_util_free_string_array(merged_text);
      } else {
        ape_test_failed("ape_io_test: ape_io_merge_par_files test 1 returned status %d, not %d as expected.\n", status,
          eOK);
      }
      ape_io_destroy_file(merged); merged = 0;
    } else {
      ape_test_failed("ape_io_test could not set up to test ape_io_merge_par_files (status is %d).\n", status);
    }
    ape_io_destroy_file(local); local = 0;
    ape_io_destroy_file(system); system = 0;

    /* Prepare for the next test by adding a new parameter to the system par file. This should result in
       the merged file looking like system par file in formatting, but with local par file values, except
       for the new parameter. Add the parameter in the middle, since this seems more likely to fail somehow. */
    system_text[9] = system_text[8];
    system_text[8] = system_text[7];
    system_text[7] = "newpar,i,h,5,0,10,\"system\"";
    status = create_test_file(file_name, system_text_cp, &system);
    ape_util_sleep(1);
    if (eOK == status) {
      status = create_test_file(file_name, local_text, &local);
    }
    if (eOK == status) {
      ApeParFile * merged = 0;
      status = ape_io_merge_par_files(system, local, &merged);
      if (eOK == status) {
        const char * expected_text[] = {
          "changepar,s,h,,min,max,\"system\"",
          "changenothing,r,h,5,0,10,\"system\"",
          "# System file layout.",
          "changevalue,r,h,5,0,10,\"system\"",
          "changetype,i,h,5,0,10,\"system\"",
          "changemode,i,h,5,0,10,\"system\"",
          "changemin,i,h,5,0,10,\"system\"",
          "newpar,i,h,5,0,10,\"system\"",
          "changemax,i,h,5,0,10,\"system\"",
          "changeenum,i,h,10,0|10,,\"system\"",
          0
        };
        char ** merged_text = 0;
        status = get_par_strings(merged, &merged_text);
        ape_test_cmp_string_array("ape_io_merge_par_files test 2", merged_text, expected_text, status, eOK);
        ape_util_free_string_array(merged_text);
      } else {
        ape_test_failed("ape_io_test: ape_io_merge_par_files test 2 returned status %d, not %d as expected.\n", status,
          eOK);
      }
      ape_io_destroy_file(merged); merged = 0;
    } else {
      ape_test_failed("ape_io_test could not set up to test ape_io_merge_par_files (status is %d).\n", status);
    }
    ape_io_destroy_file(local); local = 0;
    ape_io_destroy_file(system); system = 0;

    /* Prepare for the next test by removing some parameters from the system par file. This should result in
       the merged file looking like local par file in formatting, but with system par file type, mode etc.
       The old parameters should not be present in the merged file. */
    system_text[0] = system_text[2];
    system_text[1] = system_text[3];
    system_text[2] = system_text[6];
    system_text[3] = 0;
    status = create_test_file(file_name, system_text_cp, &system);
    ape_util_sleep(1);
    if (eOK == status) {
      status = create_test_file(file_name, local_text, &local);
    }
    if (eOK == status) {
      ApeParFile * merged = 0;
      status = ape_io_merge_par_files(system, local, &merged);
      if (eOK == status) {
        const char * expected_text[] = {
          "# Local file layout.",
          "Changemin, i, h, 5, 0, 10, \"local\"",
          "Changevalue, r, h, 5, 0, 10, \"local\"",
          0
        };
        char ** merged_text = 0;
        status = get_par_strings(merged, &merged_text);
        ape_test_cmp_string_array("ape_io_merge_par_files test 3", merged_text, expected_text, status, eOK);
        ape_util_free_string_array(merged_text);
      } else {
        ape_test_failed("ape_io_test: ape_io_merge_par_files test 3 returned status %d, not %d as expected.\n", status, eOK);
      }
      ape_io_destroy_file(merged); merged = 0;
    } else {
      ape_test_failed("ape_io_test could not set up to test ape_io_merge_par_files (status is %d).\n", status);
    }
    ape_io_destroy_file(local); local = 0;
    ape_io_destroy_file(system); system = 0;

    /* Prepare for the next test by appending one parameter to the system par file. The merged file should include the
       new parameter. Start out with an exact copy of the local parameter file. */
    memcpy(system_text, local_text, sizeof(local_text));

    /* Append one parameter to the system version. */
    system_text[sizeof(local_text)/sizeof(local_text[0]) - 1] = "newstring,s,a,,,,\"new string\"";
    system_text[sizeof(local_text)/sizeof(local_text[0])] = 0;

    status = create_test_file(file_name, system_text_cp, &system);
    if (eOK == status) {
      status = create_test_file(file_name, local_text, &local);
    }
    if (eOK == status) {
      ApeParFile * merged = 0;
      status = ape_io_merge_par_files(system, local, &merged);
      if (eOK == status) {
        const char * expected_text[] = {
          "Changeenum, i, h, 10, 0|10|100, , \"local\"",
          "Changemax, i, \"h\", 5, 0, 100, \"local\"",
          "# Local file layout.",
          "Changemin, i, h, 5, -10, 10, \"local\"",
          "Changemode, i, a, 5, 0, 10, \"local\"",
          "Changetype, r, h, 5, 0, 10, \"local\"",
          "Changevalue, r, h, 5, 0, 10, \"local\"",
          "Changenothing, r, h, 5, 0, 10, \"local\"",
          "Changepar, s, h, \"unterminated quote, min, max, \"local\"",
          "newstring,s,a,,,,\"new string\"",
          0
        };
        char ** merged_text = 0;
        status = get_par_strings(merged, &merged_text);
        ape_test_cmp_string_array("ape_io_merge_par_files test 4", merged_text, expected_text, status, eOK);
        ape_util_free_string_array(merged_text);
      } else {
        ape_test_failed("ape_io_test: ape_io_merge_par_files test 4 returned status %d, not %d as expected.\n", status, eOK);
      }
      ape_io_destroy_file(merged); merged = 0;
    } else {
      ape_test_failed("ape_io_test could not set up to test ape_io_merge_par_files (status is %d).\n", status);
    }
    ape_io_destroy_file(local); local = 0;
    ape_io_destroy_file(system); system = 0;
  }

  /* Test ape_io_revert_unlearned. */
  { const char * current_string[] = {
      " a, s,  a, current, , , prompt",
      "al, s, al, current, , , prompt",
      " h, s,  h, current, , , prompt",
      "hl, s, hl, current, , , prompt",
      " q, s,  q, current, , , prompt",
      "ql, s, ql, current, , , prompt",
      "mode, s, h, hl, , , \"Mode for automatic parameters\"",
      0
    };
    const char * previous_string[] = {
      " a, s,  a, previous, , , prompt",
      "al, s, al, previous, , , prompt",
      " h, s,  h, previous, , , prompt",
      "hl, s, hl, previous, , , prompt",
      " q, s,  q, previous, , , prompt",
      "ql, s, ql, previous, , , prompt",
      "mode, s, h, h, , , \"Mode for automatic parameters\"",
      0
    };
    const char * expected[] = {
      " a, s,  a, current, , , prompt",
      "al, s, al, current, , , prompt",
      " h, s,  h, previous, , , prompt",
      "hl, s, hl, current, , , prompt",
      " q, s,  q, previous, , , prompt",
      "ql, s, ql, current, , , prompt",
      "mode, s, h, h, , , \"Mode for automatic parameters\"",
      0
    };
    test_revert_unlearned("ape_io_revert_unlearned: previous mode is h, current mode is hl", current_string, previous_string,
      expected);
  }
  { const char * current_string[] = {
      " a, s,  a, current, , , prompt",
      "al, s, al, current, , , prompt",
      " h, s,  h, current, , , prompt",
      "hl, s, hl, current, , , prompt",
      " q, s,  q, current, , , prompt",
      "ql, s, ql, current, , , prompt",
      "mode, s, h, h, , , \"Mode for automatic parameters\"",
      0
    };
    const char * previous_string[] = {
      " a, s,  a, previous, , , prompt",
      "al, s, al, previous, , , prompt",
      " h, s,  h, previous, , , prompt",
      "hl, s, hl, previous, , , prompt",
      " q, s,  q, previous, , , prompt",
      "ql, s, ql, previous, , , prompt",
      "mode, s, h, hl, , , \"Mode for automatic parameters\"",
      0
    };
    const char * expected[] = {
      " a, s,  a, previous, , , prompt",
      "al, s, al, current, , , prompt",
      " h, s,  h, previous, , , prompt",
      "hl, s, hl, current, , , prompt",
      " q, s,  q, previous, , , prompt",
      "ql, s, ql, current, , , prompt",
      "mode, s, h, hl, , , \"Mode for automatic parameters\"",
      0
    };
    test_revert_unlearned("ape_io_revert_unlearned: previous mode is hl, current mode is h", current_string, previous_string,
      expected);
  }
  { const char * current_string[] = {
      " a, s,  a, current, , , prompt",
      "al, s, al, current, , , prompt",
      " h, s,  h, current, , , prompt",
      "hl, s, hl, current, , , prompt",
      " q, s,  q, current, , , prompt",
      "ql, s, ql, current, , , prompt",
      0
    };
    const char * previous_string[] = {
      " a, s,  a, previous, , , prompt",
      "al, s, al, previous, , , prompt",
      " h, s,  h, previous, , , prompt",
      "hl, s, hl, previous, , , prompt",
      " q, s,  q, previous, , , prompt",
      "ql, s, ql, previous, , , prompt",
      0
    };
    const char * expected[] = {
      " a, s,  a, current, , , prompt",
      "al, s, al, current, , , prompt",
      " h, s,  h, previous, , , prompt",
      "hl, s, hl, current, , , prompt",
      " q, s,  q, previous, , , prompt",
      "ql, s, ql, current, , , prompt",
      0
    };
    test_revert_unlearned("ape_io_revert_unlearned: no previous mode, no current mode", current_string, previous_string,
      expected);
  }
  test_apply_command_line();

  /* Test ape_io_check_file_format. */
  { const char * par_string[] = {
      "p0, s, h, , , , parameter 0",
      "p0, s, h, , , , parameter 0",
      "p1, s, h, , , , parameter 1",
      "p2, s, h, , , , parameter 2",
      "p1, s, h, , , , parameter 1 duplicated",
      "p3, s, h, , , , parameter 3",
      "p2, s, h, , , , parameter 2",
      0
    };
    ApeParFile * file = 0;
    status = create_test_file("dups.par", par_string, &file);
    if (eOK == status) {
      status = ape_io_check_file_format(file, 0);
      if (eParameterDuplicated != status) {
        ape_test_failed("ape_io_check_file_format(file containing duplicates) returned status %d, not %d as expected.\n",
          status, eParameterDuplicated);
      }
    } else {
      ape_test_failed("Unable to set up to test ape_io_check_file_format (status was %d.)\n", status);
    }
    ape_io_destroy_file(file); file = 0;
  }

  /* Test ability to handle a file larger than the buffer size. */
  { const char * par_string[] = {
      "p0, s, h, \"A value for p0\", , , parameter 0",
      "p1, s, h, \"A value for p1\", , , parameter 1",
      "p2, s, h, \"A value for p2\", , , parameter 2",
      "p3, s, h, \"A value for p3\", , , parameter 3",
      0
    };
    char large_value[APE_BUF_SIZE] = "";
    size_t ii = 0;
    ApeParFile * file = 0;

    /* Create a large value string. */
    for (ii = 0; ii < APE_BUF_SIZE - 1; ++ii) {
      large_value[ii] = 'x';
    }
    large_value[APE_BUF_SIZE - 1] = '\0';

    status = create_test_file("ape_large.par", par_string, &file);
    if (eOK != status) {
      ape_test_failed("Unable to set up to test files larger than buffer size (status was %d.)\n", status);
    }
    if (eOK == status) {
      ApeListIterator par_itor = 0;
      ApePar * par = 0;
      status = ape_io_find_par("p2", file, &par_itor);
      if (eOK == status) {
        par = (ApePar *) ape_list_get(par_itor);
        if (0 == par) {
          status = eNullPointer;
          ape_test_failed("Unable to get ApePar* parameter p2 when testing files larger than buffer size (status was %d.)\n",
            status);
        }
      } else {
        ape_test_failed("Unable to find parameter p2 when testing files larger than buffer size (status was %d.)\n", status);
      }

      if (eOK == status) {
        /* Assign a large value to parameter p2. */
        status = ape_par_set_value_string(par, large_value);
        if (eOK != status) {
          ape_test_failed("Unable to assign string larger than buffer size to parameter (status was %d.)\n", status);
        }
      }

      /* Verify that the string value is correct. */
      if (eOK == status) {
        char * set_value = 0;
        status = ape_par_get_string_case(par, &set_value, eDefaultCase);
        if (eOK != status) {
          ape_test_failed("Unable to read string parameter larger than buffer size (status was %d.)\n", status);
        }
        if (eOK == status && 0 != strncmp(large_value, set_value, APE_BUF_SIZE)) {
          ape_test_failed("Value read disagreed with value written for parameter larger than buffer size (status was %d.)\n",
            status);
        }
        free(set_value); set_value = 0;
      }
    }

    /* Write the file with the large line. */
    status = ape_io_write(file, 1);
    if (eOK != status) {
      ape_test_failed("Unable to write file with parameter value larger than buffer size (status was %d.)\n", status);
    }
    ape_io_destroy_file(file); file = 0;

    /* Now read the file with the large line. */
    status = ape_io_read("ape_large.par", &file);
    if (eLineTooLong != status) {
      ape_test_failed("ape_io_read for a file with parameter value larger than buffer size returned status "
        "%d, not %d as expected.\n", status, eLineTooLong);
    } else {
      /* Expected problem, so ignore it and continue. The other lines should be correct. */
      status = eOK;
    }

    /* Compare strings in written file to the expected values. */
    if (eOK == status) {
      const char * expected_text[] = {
        "p0, s, h, \"A value for p0\", , , parameter 0",
        "p1, s, h, \"A value for p1\", , , parameter 1",
        0,
        "p3, s, h, \"A value for p3\", , , parameter 3",
        0
      };
      char ** actual_text = 0;

      /* Construct expected value for the long parameter (p2). */
      char line[APE_BUF_SIZE] = "";
      strcpy(line, "p2, s, h, \"");
      strncat(line, large_value, APE_BUF_SIZE - strlen(line) - 1);

      /* Assign expected line for p2. */
      expected_text[2] = line;

      status = get_par_strings(file, &actual_text);
      ape_test_cmp_string_array("Comparing read to written parameter value larger than buffer size", actual_text, expected_text,
        status, eOK);
      ape_util_free_string_array(actual_text);
    }
    ape_io_destroy_file(file); file = 0;
  }
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_io.c,v $
 * Revision 1.81  2011/03/02 22:18:58  irby
 * Fixed previous revision, which did not allow for the scenario in which a
 * command line *value* ambiguously matches a parameter name.  Added several
 * more tests checking ambiguity.
 *
 * Revision 1.80  2011/02/18 19:56:29  irby
 * In classify_arg(), check for ambiguity in parameters entered on the
 * command line: allow partial matching when there is no ambiguity, or
 * return a non-zero status (eAmbiguousParName) when there is.  Add tests
 * for partial matching and ambiguity.
 *
 * Revision 1.79  2010/11/22 16:24:34  jpeachey
 * Add parsing support for + and - at end of parameter names. These are
 * interpreted as yes and no, i.e. parameter+ means the same as parameter=yes.
 *
 * Revision 1.78  2009/07/21 20:59:31  peachey
 * Improve on encapsulation of stat: add function get_mod_timeto get the modification time in OS independent way.
 *
 * Revision 1.77  2009/07/21 19:08:00  peachey
 * Use _stat instead of stat on Windows.
 *
 * Revision 1.76  2009/07/09 21:02:01  peachey
 * Add cast to prevent warnings converting size_t to int.
 *
 * Revision 1.75  2009/07/08 19:37:38  peachey
 * Add tests for handling lines longer than the maximum length
 * (8192) and correct code to detect and deal with such lines properly.
 *
 * Revision 1.74  2009/06/23 17:53:37  peachey
 * Add a test for ape_io_merge_par_files in which the local and system
 * parameters agree in type, min, max, and position, but system parameter
 * file is more recent, causing local hidden parameters to be trumped by
 * the system file. This is to catch the corner case in which
 * rough_compare_files would otherwise conclude the files did not need the
 * full merge behavior.
 *
 * Revision 1.73  2009/06/11 19:13:13  peachey
 * Check modification times of files and overwrite local values
 * with system values if the system parameter file is more recent.
 *
 * Revision 1.72  2008/07/17 14:55:07  peachey
 * In ape_io_read_file_path, relax restriction that directories in the
 * path be non-blank. Skip over blank directories instead.
 *
 * Revision 1.71  2007/11/16 18:18:23  peachey
 * Add check for automatic parameters without a mode parameter to the
 * function ape_io_check_file_format. However, do not automatically enable
 * debugging messages in ape_io_check_file_format.
 *
 * Revision 1.70  2007/11/12 16:54:52  peachey
 * Check status more specifically when calling ape_io_get_default_mode.
 * Ignore error only if the mode parameter was not found.
 *
 * Revision 1.69  2007/09/18 16:40:56  peachey
 * Add more casts to satisfy gcc.
 *
 * Revision 1.68  2007/09/14 18:35:26  peachey
 * Change pointer type to non-const so that it can be used as first argumentto memset without a warning from Visual Studio.
 *
 * Revision 1.67  2007/08/27 13:53:58  peachey
 * Fix bug in the merging optimization code that prevented it
 * from detecting situations where the system parameter file had extra
 * parameters appended only.
 *
 * Revision 1.66  2007/05/15 20:57:58  peachey
 * Plug a memory leak in a newer piece of optimized code.
 *
 * Revision 1.65  2007/05/15 20:52:48  peachey
 * Additional speed optimizations:
 * o Reduce number of calls to ape_list_end by calling it once and storing
 *   the returned iterator for use in loop logic when it is safe to do so.
 * o Make ape_list_destroy more efficient by iterating through the list
 *   back to front once and freeing the elements instead of repeatedly
 *   calling ape_list_remove_entry.
 *
 * Revision 1.64  2007/05/03 19:28:24  peachey
 * Slight clean-up/readability improvement.
 *
 * Revision 1.63  2007/05/03 19:24:03  peachey
 * Ignore differences in comments and blank lines when checking for consistency.
 *
 * Revision 1.62  2007/05/03 17:20:39  peachey
 * Improve efficiency of code that merges local and system parameter files.
 * Check for exact consistency between the two, and if consistent, skip the
 * more complicated and time-consuming merging process.
 *
 * Revision 1.61  2007/02/21 18:52:48  peachey
 * Tweak cast to silence warnings on Unix, without creating new warnings
 * on Windows.
 *
 * Revision 1.60  2007/02/21 18:42:44  peachey
 * Add a cast and change a declaration to silence warnings from Visual Studio.
 *
 * Revision 1.59  2007/02/01 18:55:41  peachey
 * In ape_io_check_file_format, do not let fact that a parameter
 * is or is not enumerated constitute an error.
 *
 * Revision 1.58  2006/11/24 19:50:51  peachey
 * On windows, rename fails if target file is present, so remove it first.
 *
 * Revision 1.57  2006/11/24 15:57:28  peachey
 * Correct a comment.
 *
 * Revision 1.56  2006/11/24 15:54:40  peachey
 * Change command line parsing: instead of a specific exclusion for the case
 * of = signs inside square brackets, simply restrict error checking for invalid
 * parameter = value patterns to legal parameter names (alphanumeric, _, - and .).
 *
 * Revision 1.55  2006/11/08 22:19:20  peachey
 * Revamp command line handling to detect invalid parameter names, and
 * to detect parameter duplications when positional and explicit parameters
 * specify the same parameter twice or more.
 *
 * Revision 1.54  2006/11/03 14:17:08  peachey
 * Update comments.
 *
 * Revision 1.53  2006/08/25 19:43:44  peachey
 * Add tests for undefined and infinite numerical parameters.
 * Be more careful about returning if test set-up steps fail.
 *
 * Revision 1.52  2006/07/05 20:05:38  peachey
 * Add more specific error messages for command line parsing problems.
 *
 * Revision 1.51  2006/06/23 01:57:51  peachey
 * Use temporary file when writing output, then move file over top of real par file.
 *
 * Revision 1.50  2006/06/23 01:13:01  peachey
 * Add auto-correct feature to ape_io_merge_par_files; if a local parameter is
 * corrupt somehow, destroy it and substitute the system parameter of the same name.
 *
 * Revision 1.49  2006/06/20 02:54:00  peachey
 * Change test to reflect fact that white space is now preserved by
 * assignments.
 *
 * Revision 1.48  2006/06/16 01:19:11  peachey
 * Enable debug mode when checking parameter file.
 *
 * Revision 1.47  2006/06/08 02:25:01  peachey
 * Change annoying error message when any file can't be opened to a
 * debug statement.
 *
 * Revision 1.46  2006/06/06 13:29:24  peachey
 * Explicitly flag parameters whose value was set on command line.
 *
 * Revision 1.45  2006/05/31 01:41:56  peachey
 * Rename ape_par_set_string to ape_par_set_field.
 *
 * Revision 1.44  2006/05/31 01:36:51  peachey
 * Rename ape_par_get_string to ape_par_get_field and ape_par_get_string_array to
 * ape_par_get_field_array.
 *
 * Revision 1.43  2006/05/23 16:26:57  peachey
 * Add read-only flag to par file structure. Add force_write arguement to
 * ape_io_write for overriding read-only status.
 *
 * Revision 1.42  2006/05/19 17:37:19  peachey
 * Update TODO which was handled.
 *
 * Revision 1.41  2006/05/18 03:15:44  peachey
 * Add ape_io_set_file_name to allow the file name used by a
 * file object to be modified.
 *
 * Revision 1.40  2006/05/17 02:19:49  peachey
 * Add ape_io_get_par_cont, for getting the container of parameters from a file.
 *
 * Revision 1.39  2006/05/16 19:46:17  peachey
 * Add ape_io_get_pfiles, and automatically perform default
 * initialization if ape_io_get_loc_path or ape_io_get_sys_path are called
 * and PFILES is not defined.
 *
 * Revision 1.38  2006/05/16 17:22:55  peachey
 * Check for duplicated parameters in ape_io_check_file_format.
 *
 * Revision 1.37  2006/05/16 15:17:35  peachey
 * Be more careful about file names when merging parameter files, and when writing
 * output parameter files.
 *
 * Revision 1.36  2006/05/16 14:55:20  peachey
 * Write a more detailed test parameter file.
 *
 * Revision 1.35  2006/05/12 17:22:55  peachey
 * Add ape_io_revert_unlearned, to forget unlearned parameters when requested.
 *
 * Revision 1.34  2006/05/12 03:30:18  peachey
 * Add function ape_io_get_default_mode, for getting mode parameter value from file.
 *
 * Revision 1.33  2006/05/12 00:24:37  peachey
 * When applying command line argument, use ape_par_set_value_string to apply the value,
 * so that it will be flagged as modified.
 *
 * Revision 1.32  2006/05/10 16:48:50  peachey
 * Generalize format tests to allow them to be run with or without
 * value checking, to let them be used both in contexts where the whole parameter
 * needs to be checked (checking user input) and contexts in which a bad value
 * is not a show-stopper, because the user or client may yet correct the problem.
 *
 * Revision 1.31  2006/05/09 19:31:41  peachey
 * Add ape_io_check_file_format, for checking the whole parameter file.
 *
 * Revision 1.30  2006/05/06 00:40:30  peachey
 * Correct bug in which blank lines and comment-only were not
 * properly skipped when interpreting command line arguments.
 *
 * Revision 1.29  2006/05/02 23:48:29  peachey
 * Add ape_io_apply_command_line for handling command line parameters,
 * directly imported from pil's PIL_allow_spaces_cmndarg2value.
 *
 * Revision 1.28  2006/05/01 19:28:25  peachey
 * Fix accidental mixed declarations and code.
 *
 * Revision 1.27  2006/05/01 19:13:58  peachey
 * Make tests for ape_io_merge_files more demanding.
 *
 * Revision 1.26  2006/05/01 15:13:45  peachey
 * Fix some memory leaks.
 *
 * Revision 1.25  2006/05/01 14:51:57  peachey
 * Change signature of functions applied by apply_operation so that they
 * can modify the iterators which are supplied to the function. Complete
 * implementation and tests for par_synch and par_prune, and use them to
 * finish implementation of ape_io_merge_files.
 *
 * Revision 1.24  2006/04/28 14:40:41  peachey
 * Add ape_io_merge_par_files and beginnings of a test for it.
 * Add helper function create_test_file for creating files from arrays of strings.
 *
 * Revision 1.23  2006/04/28 01:24:04  peachey
 * Add a couple casts to get the constness right when calling free.
 * Enable al mode, and remove eFileOpenError, using eFileReadError in lieu thereof.
 * The latter change was to make error messages consistent between Windows and
 * Unix, because although the two both fail to open a directory as if it were
 * a file, they fail at different points in the process of opening the file.
 *
 * Revision 1.22  2006/04/24 13:28:44  peachey
 * Add ape_io_clone_file and the beginnings of ape_io_merge_par_files.
 *
 * Revision 1.21  2006/04/22 01:38:26  peachey
 * Add internal function apply_operation, a utility to help with merging par files.
 *
 * Revision 1.20  2006/04/21 14:52:49  peachey
 * Change signature of ape_io_find_par so that it returns an iterator
 * instead of the parameter itself.
 *
 * Revision 1.19  2006/04/21 14:26:50  peachey
 * Do case-insensitive matches in ape_io_find_par.
 *
 * Revision 1.18  2006/04/21 01:34:49  peachey
 * Patch memory leaks in unit test.
 *
 * Revision 1.17  2006/04/21 01:30:59  peachey
 * Add and test ape_io_find_par, to locate a parameter by name.
 *
 * Revision 1.16  2006/04/19 15:41:32  peachey
 * Use ape_util_atexit instead of directly calling atexit.
 *
 * Revision 1.15  2006/04/19 03:06:16  peachey
 * Add name of function to debug message.
 *
 * Revision 1.14  2006/04/19 02:54:57  peachey
 * Fix memory leak.
 *
 * Revision 1.13  2006/04/19 02:43:28  peachey
 * Removed array of parameters which is no longer needed after refactoring
 * ape_io_write.
 *
 * Revision 1.12  2006/04/19 02:30:10  peachey
 * Refactor ape_io_write to use ApeParFile instead of ApeList.
 *
 * Revision 1.11  2006/04/19 01:33:03  peachey
 * Change signature of ape_io_read so that it returns an ApeParFile rather than
 * an ApeList of parameters.
 *
 * Revision 1.10  2006/04/19 00:46:15  peachey
 * Removed unused variable.
 *
 * Revision 1.9  2006/04/19 00:38:02  peachey
 * Make ape_io_read able to read parameter files written on either
 * Unix or Windows. Improve reading code and unit test.
 *
 * Revision 1.8  2006/04/15 03:06:26  peachey
 * Add ape_io_destroy_file for cleaning up ApeParFiles.
 *
 * Revision 1.7  2006/04/15 02:47:39  peachey
 * Revamp ape_io_read to reduce nesting levels.
 *
 * Revision 1.6  2006/04/15 00:47:07  peachey
 * Add more support for dealing with parameter files. Specifically, a
 * struct ApeParFile, which contains a file name and a container of
 * parameters, and two new functions:
 * o ape_io_get_file_name: retrieves the file name from an ApeParFile.
 * o ape_io_read_file_path: searches a path for a file. If it is found,
 *   open it and load its parameters.
 *
 * Revision 1.5  2006/04/14 14:24:36  peachey
 * Add cast to clarify to Windows compiler that in this case it's OK to convert
 * size_t to unsigned long.
 *
 * Revision 1.4  2006/04/14 14:16:05  peachey
 * Add void to empty paren in destroy_paths, just to be extra clear.
 *
 * Revision 1.3  2006/04/14 03:53:48  peachey
 * Arrange for auto-cleanup of pfile/path formulation.
 *
 * Revision 1.2  2006/04/14 03:43:20  peachey
 * Add facilities for handling PFILES: ape_io_set_pfiles, ape_io_get_loc_path,
 * ape_io_get_sys_path.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
