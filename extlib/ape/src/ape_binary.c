/** \file ape_binary.c
    \brief Implementation of facilities to support binaries using ape, including especially pget, pset etc..
    \author James Peachey, HEASARC/EUD/GSFC.
*/

#include "ape/ape_binary.h"
#include "ape/ape_error.h"
#include "ape/ape_io.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_test.h"
#include "ape/ape_trad.h"
#include "ape/ape_util.h"

#include <limits.h>
#include <stddef.h>
#include <stdio.h>
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

static int check_args(int argc, char ** argv, int min_argc, int max_argc) {
  int status = eOK;

  if (0 == argv) status = eNullPointer;

  if (eOK == status) {
    if (0 > min_argc) min_argc = INT_MIN;
    if (0 > max_argc) max_argc = INT_MAX;

    if (min_argc > argc || max_argc < argc) status = eInvalidArgument;
  }

  return status;
}

static void usage(const char * usage) {
  if (0 == usage) usage = "";
  ape_msg_error("Usage: %s\n", usage);
}

static int adjust_cmd(int * argc, char *** argv, char * is_file) {
  int status = eOK;

  if (0 != is_file) *is_file = 0;
  else status = eNullPointer;

  if (eOK == status && (0 == argc || 0 == argv)) status = eNullPointer;
  if (eOK == status && 1 > *argc) status = eInvalidArgument;

  if (eOK == status) {
    /* Check whether there are at least two arguments, the first of which is -f and second of which is
       presumably a parameter file name. */
    if (0 == strcmp((*argv)[0], "-f")) {
      if (1 < *argc) {
        --*argc; ++*argv; *is_file = 1;
      } else {
        status = eInvalidArgument;
      }
    } else {
      const size_t dot_par_len = strlen(".par");
      const size_t arg_len = strlen((*argv)[0]);
      /* Check if the first argument ends in ".par", indicating it may be a parameter file name. */
      if (arg_len >= dot_par_len && 0 == strcmp((*argv)[0] + arg_len - dot_par_len, ".par")) {
        /* Require also that the parameter file name includes a slash or backslash, i.e. a absolute or relative path. */
        if (0 != strstr((*argv)[0], QD) || 0 != strstr((*argv)[0], QD2))
          *is_file = 1;
      }
    }
  }

  return status;
}

static int read_par_file(int argc, char ** argv, char is_file, ApeParFile ** par_file) {
  int status = eOK;

  if (0 != par_file) *par_file = 0;
  else status = eNullPointer;

  if (eOK == status && 0 == argv) status = eNullPointer;
  if (eOK == status && 1 > argc) status = eInvalidArgument;

  if (eOK == status) {
    if (0 == is_file) {
      /* argv[0] is a binary name, not a parameter file name, so perform standard initialization. */
      status = ape_trad_init(argc, argv);

      /* Get the current parameter file, even if there was an error. */
      ape_trad_get_current(par_file);
    } else {
      /* argv[0] is a parameter file name, not a binary name, so use it to read parameters, and apply command line explicitly. */
      status = ape_io_read(argv[0], par_file);
      if (eOK == status) {
        /* Modify the parameter file using the command line arguments. */
        status = ape_io_apply_command_line(*par_file, argc - 1, argv + 1);
      }
    }
  }
  return status;
}

static void redirect_prompts(char open) {
  static char s_redirected = 0;
  static FILE * s_prompt_stream = 0;
  if (0 != open) {
    if (0 == s_redirected) {
#ifdef WIN32
      s_prompt_stream = stdout;
#else
      /* All prompts should go to the terminal since so that utilities
         will work via backticks from Perl scripts. */
      s_prompt_stream = fopen("/dev/tty","w");
#endif
      if (0 != s_prompt_stream) {
        ape_par_redirect_prompt_stream(s_prompt_stream);
        s_redirected = 1;
      }
    }
  } else {
    if (0 != s_redirected) {
      ape_par_redirect_prompt_stream(stdout);
      if (0 != s_prompt_stream) fclose(s_prompt_stream);
      s_redirected = 0;
    }
  }
}

#ifdef __cplusplus
extern "C" {
#endif

int ape_binary_run(int argc, char ** argv, ApeBinary * binary) {
  int status = eOK;

  /* Check arguments. */
  if (0 == argv || 0 == binary || 0 == binary->name || 0 == binary->usage) {
    status = eNullPointer;
  }

  if (eOK == status && 1 > argc) {
    status = eInvalidArgument;
  }

  if (eOK == status) {
    /* Interpret environment. */
    status = ape_util_interpret_env();
  }

  if (eOK == status) {
    /* Skip the name of the executable. */
    --argc; ++argv;
  }

  /* Check command line arguments. */
  if (eOK == status && (binary->min_argc > argc || (0 < binary->max_argc && binary->max_argc < argc))) {
    ape_msg_error("Usage: %s %s\n", binary->name, binary->usage);
    status = eInvalidArgument;
  } else {
    /* Make sure prompts go to the right place. */
    redirect_prompts(1);

    status = binary->func(argc, argv);

    /* Close down prompt stream. */
    redirect_prompts(0);
  }

  return eOK != status ? 1 : eOK;
}

int ape_binary_pget(int argc, char ** argv) {
  int status = eOK;
  char is_file = 0;
  ApeParFile * par_file = 0;

  /* Interpret environment. */
  status = ape_util_interpret_env();

  /* Interpret (and adjust if necessary) the command line to see if a parameter file name was given. */
  if (eOK == status) status = adjust_cmd(&argc, &argv, &is_file);

  /* Check arguments for basic validity. */
  if (eOK == status) status = check_args(argc, argv, 1, -1);

  /* Issue usage message and return if basic conditions were not met. */
  if (eOK != status) {
    usage("pget [ -f ] tool-or-par-file-name [ parameter-name-list ]");
    return status;
  }

  /* Load parameters from file(s). Do not apply other command line arguments, which name parameters to get. */
  if (eOK == status) status = read_par_file(1, argv, is_file, &par_file);
  if (eOK != status) {
    ape_msg_error("pget: unable to initialize tool/parameter file \"%s\".\n", argv[0]);
  }

  if (eOK == status) {
    int idx = 0;
    /* Iterate over command line arguments. */
    for (idx = 1; idx != argc; ++idx) {
      ApeListIterator itor = 0;
      ApePar * par = 0;
      char * value = 0;
      /* Get parameter named by this argument. */
      int local_status = ape_io_find_par(argv[idx], par_file, &itor);
      if (eOK == local_status) {
        par = (ApePar *) ape_list_get(itor);
      }
      if (eOK == local_status) {
        local_status = ape_par_get_string(par, &value);
      }
      if (eOK == local_status) {
        /* Display parameter value. */
        ape_msg_out("%s\n", value);
      } else {
        /* Display error message for parameters which could not be obtained. */
        ape_msg_error("pget: could not get parameter \"%s\"\n", argv[idx]);
      }

      /* If there was a problem getting any parameters, change the overall status. */
      status = eOK == status ? local_status : status;

      /* Clean up. */
      free(value); value = 0;
    }
  }

  /* Clean up without saving anything. */
  if (0 == is_file) ape_trad_close(0);
  else ape_io_destroy_file(par_file);

  return status;
}

int ape_binary_plist(int argc, char ** argv) {
  int status = eOK;
  char is_file = 0;
  ApeParFile * par_file = 0;
  const char * file_name = 0;
  ApeList * par_cont = 0;

  /* Interpret environment. */
  status = ape_util_interpret_env();

  /* Interpret (and adjust if necessary) the command line to see if a parameter file name was given. */
  if (eOK == status) status = adjust_cmd(&argc, &argv, &is_file);

  /* Check arguments for basic validity. */
  if (eOK == status) status = check_args(argc, argv, 1, 1);

  /* Issue usage message and return if basic conditions were not met. */
  if (eOK != status) {
    usage("plist [ -f ] tool-or-par-file-name");
    return status;
  }

  /* Load parameters from file(s). Do not apply other command line arguments, if any. */
  if (eOK == status) status = read_par_file(1, argv, is_file, &par_file);
  if (eOK != status) {
    ape_msg_error("plist: unable to initialize tool/parameter file \"%s\".\n", argv[0]);
  }

  if (eOK == status) {
    status = ape_io_get_file_name(par_file, &file_name);
  }

  if (eOK == status) {
    ape_msg_out("Parameters for %s\n", file_name);
  }

  if (eOK == status) {
    status = ape_io_get_par_cont(par_file, &par_cont);
  }

  if (eOK == status) {
    ApeListIterator itor = 0;
    for (itor = ape_list_begin(par_cont); itor != ape_list_end(par_cont); itor = ape_list_next(itor)) {
      char * name = 0;
      char * value = 0;
      char * prompt = 0;
      char * new_name = 0;
      char * new_value = 0;
      char auto_mode[APE_PAR_MODE_CODE_LEN] = "";
      char mode[APE_PAR_MODE_CODE_LEN] = "";
      ApePar * par = (ApePar *) ape_list_get(itor);

      /* Get name field. */
      status = ape_par_get_field(par, eName, &name);
      if (eOK == status) {
        /* Get value field. */
        status = ape_par_get_field(par, eValue, &value);
      }
      if (eOK == status) {
        /* Get prompt field. */
        status = ape_par_get_field(par, ePrompt, &prompt);
      }
      if (eOK == status) {
        /* Get mode parameter. */
        status = ape_io_get_default_mode(par_file, auto_mode);
        if (eParNotFound == status) {
          status = eOK;
          strcpy(auto_mode, "ql");
        }
      }
      if (eOK == status) {
        /* Get mode. */
        status = ape_par_get_eff_mode(par, auto_mode, mode);
      }
      if (eOK == status && 'h' == *mode) {
        /* For hidden mode use parentheses. */
        status = ape_util_cat_string("(", name, &new_name);
        if (eOK == status) {
          status = ape_util_cat_string(value, ")", &new_value); 
        }
      }
      if (eOK == status) {
        if (0 != new_value) {
          ape_msg_out("%13s = %-16s %s\n", new_name, new_value, prompt);
        } else {
          ape_msg_out("%13s = %-16s %s\n", name, value, prompt);
        }
      }

      /* Clean up. */
      free(new_value); new_value = 0;
      free(new_name); new_name = 0;
      free(prompt); prompt = 0;
      free(value); value = 0;
      free(name); name = 0;
    }
  }

  /* Clean up without saving anything. */
  if (0 == is_file) ape_trad_close(0);
  else ape_io_destroy_file(par_file);

  return status;
}

int ape_binary_pquery(int argc, char ** argv) {
  int status = eOK;
  int new_argc = 0;
  char ** new_argv = 0;
  char * par_name = 0;
  char is_file = 0;
  ApeParFile * par_file = 0;
  ApeParFile * orig_par_file = 0;
  ApeList * par_cont = 0;
  ApePar * par = 0;
  char * value = 0;
  char default_mode[APE_PAR_MODE_CODE_LEN] = "";

  /* Interpret environment. */
  status = ape_util_interpret_env();

  /* Interpret (and adjust if necessary) the command line to see if a parameter file name was given. */
  if (eOK == status) status = adjust_cmd(&argc, &argv, &is_file);

  /* Check arguments for basic validity. */
  if (eOK == status) status = check_args(argc, argv, 2, -1);

  /* Issue usage message and return if basic conditions were not met. */
  if (eOK != status) {
    usage("pquery2 [ -f ] tool-or-par-file-name parameter-name [ command line arguments ]");
    return status;
  }

  /* Make sure prompts go to the right place. */
  redirect_prompts(1);

  /* Make new copy of command line arguments. */
  if (eOK == status) {
    new_argv = (char **) calloc(argc - 1, sizeof(char *));
    if (0 == new_argv) status = eDynAllocFailed;
  }

  /* First, load parameters from file(s) without applying command line arguments. */
  if (eOK == status) {
    new_argv[0] = argv[0];
    new_argc = 1;
    status = read_par_file(new_argc, new_argv, is_file, &par_file);
    if (eOK != status)
      ape_msg_error("pquery2: unable to initialize tool/parameter file \"%s\".\n", new_argv[0]);
  }

  /* Clone the original parameter file. */
  if (eOK == status) {
    status = ape_io_clone_file(par_file, &orig_par_file);
  }

  /* Clean up without saving anything. */
  if (0 == is_file) ape_trad_close(0);
  else ape_io_destroy_file(par_file);

  if (eOK == status) {
    int idx = 0;

    new_argc = argc - 1;

    /* Rearrange arguments: argv[0] == name of executable/pfile. argv[1] == name of parameter for which to pquery. */
    new_argv[0] = argv[0];
    par_name = argv[1];

    /* Remaining arguments are passed to the command line handler for the parameter file. */
    for (idx = 2; idx < argc; ++idx) {
      new_argv[idx - 1] = argv[idx];
    }
  }

  /* Load parameters from file(s), applying command line arguments as needed. */
  if (eOK == status) {
    status = read_par_file(new_argc, new_argv, is_file, &par_file);
    if (eOK != status)
      ape_msg_error("pquery2: unable to initialize tool/parameter file \"%s\".\n", new_argv[0]);
  }

  if (eOK == status) {
    status = ape_io_get_par_cont(par_file, &par_cont);
  }

  if (eOK == status) {
    ApeListIterator itor = 0;
    status = ape_io_find_par(par_name, par_file, &itor);
    if (eOK == status) {
      par = (ApePar *) ape_list_get(itor);
    } else {
      ape_msg_error("pquery2: unable to prompt for parameter \"%s\".\n", par_name);
    }
  }

  if (eOK == status) {
    status = ape_io_get_default_mode(par_file, default_mode);
    if (eParNotFound == status) {
      status = eOK;
      strcpy(default_mode, "ql");
    }
  }

  if (eOK == status) {
    status = ape_par_query(par, default_mode);
  }

  if (eOK == status) {
    status = ape_par_get_string(par, &value);
  }

  if (eOK == status) {
    ape_msg_out("%s\n", value);
  }

  if (eOK == status) {
    /* Revert unlearned parameters. This gives the "pquery2" behavior, whereas pquery does write
       even unlearned parameters. */
    status = ape_io_revert_unlearned(par_file, orig_par_file);
  }

  /* Write parameter file with final values. */
  if (eOK == status) status = ape_io_write(par_file, 1);

  /* Clean up. */
  free(value); value = 0;
  free(new_argv); new_argv = 0;

  /* Clean up without saving anything. */
  if (0 == is_file) ape_trad_close(0);
  else ape_io_destroy_file(par_file);
  par_file = 0;

  /* Clean up copy of parameter file used to revert hidden parameters to previous values. */
  ape_io_destroy_file(orig_par_file); orig_par_file = 0;

  /* Close down prompt stream. */
  redirect_prompts(0);

  return status;
}

int ape_binary_pset(int argc, char ** argv) {
  int status = eOK;
  char is_file = 0;
  ApeParFile * par_file = 0;

  /* Interpret environment. */
  status = ape_util_interpret_env();

  /* Interpret (and adjust if necessary) the command line to see if a parameter file name was given. */
  if (eOK == status) status = adjust_cmd(&argc, &argv, &is_file);

  /* Check arguments for basic validity. */
  if (eOK == status) status = check_args(argc, argv, 1, -1);

  /* Issue usage message and return if basic conditions were not met. */
  if (eOK != status) {
    usage("pset [ -f ] tool-or-par-file-name [ parameter-args ]");
    return status;
  }

  /* Make sure prompts go to the right place. */
  redirect_prompts(1);

  /* Load parameters from file(s), applying command line arguments as needed. */
  if (eOK == status) status = read_par_file(argc, argv, is_file, &par_file);
  if (eOK != status) {
    ape_msg_error("pset: unable to initialize tool/parameter file \"%s\".\n", argv[0]);
  }

  /* Check command line arguments. */
  /* Check whether only 1 argument was supplied, in which case it is necessary to prompt for each parameter in turn. */
  if (eOK == status && 1 == argc) {
    ApeList * par_cont = 0;
    ApeListIterator itor = 0;

    /* Get the container of parameters from this file. */
    status = ape_io_get_par_cont(par_file, &par_cont);
    if (eOK == status) {
      /* Iterate over all paramaters. */
      for (itor = ape_list_begin(par_cont); eOK == status && itor != ape_list_end(par_cont); itor = ape_list_next(itor)) {
        ApePar * par = (ApePar *) ape_list_get(itor);
        /* Prompt for each parameter, regardless of mode, and rewrite parameter file each time if OK. */
        status = ape_par_prompt(par);
        if (eOK == status) {
          status = ape_io_write(par_file, 1);
        } else {
          char * par_name = 0;
          if (eOK == ape_par_get_field(par, eName, &par_name))
            ape_msg_error("pset: problem prompting for parameter \"%s\".\n", par_name);
          free(par_name); par_name = 0;
        }
      }
    }
    if (eOK != status) ape_msg_error("pset: problem(s) with one or more parameter prompts.\n");
  }

  /* Write parameter file with final values if no error. */
  if (eOK == status) status = ape_io_write(par_file, 1);

  if (0 == is_file) ape_trad_close(0); /* Don't write while closing. */
  else ape_io_destroy_file(par_file);

  /* Close down prompt stream. */
  redirect_prompts(0);

  return status;
}

int ape_binary_punlearn(int argc, char ** argv) {
  int status = eOK;
  const char * usage_msg = "punlearn [ -f ] tool-or-par-file-name [ [ -f ] tool-or-par-file-name ... ] ";

  /* Interpret environment. */
  status = ape_util_interpret_env();

  if (0 >= argc) usage(usage_msg);

  while (eOK == status && 0 < argc) {
    char is_file = 0;
    ApeParFile * sys_pfile = 0;
    const char * sys_name = 0;
    const char * unlearn_name = 0;

    /* Interpret (and adjust if necessary) the command line to see if the first argument is the name of a
       parameter file. */
    if (eOK == status) status = adjust_cmd(&argc, &argv, &is_file);

    /* Check arguments for basic validity, and display usage if conditions not met. */
    if (eOK == status) status = check_args(argc, argv, 1, -1);
    if (eOK != status) usage(usage_msg);

    /* Perform standard initialization, without arguments since they are all names of par files/tools.
       Ignore status in case there is an error opening the local parameter file. */
    if (eOK == status) ape_trad_init(1, argv);

    /* Get the system parameter file. */
    if (eOK == status) {
      status = ape_trad_get_sys_pfile(&sys_pfile);
      if (eOK != status) ape_msg_error("punlearn: cannot find a system parameter file for \"%s\".\n", argv[0]);
    }

    /* Get the name of the system parameter file. */
    if (eOK == status) status = ape_io_get_file_name(sys_pfile, &sys_name);

    /* Get the name of the parameter file to "unlearn". */
    if (eOK == status) {
      /* If the current first argument is the name of a parameter file, 0th argument is the file name. */
      if (0 != is_file) {
        unlearn_name = argv[0];
      } else {
        ApeParFile * current = 0;
        /* Saved file name is the name of the current parameter file, as opened by the standard initialization. */
        status = ape_trad_get_current(&current);
        if (eOK == status) {
          status = ape_io_get_file_name(current, &unlearn_name);
        }
      }
    }

    if (eOK == status) {
      /* Compare the names of the system and the "unlearn" file. If they are distinct, proceed with unlearn. */
      if (0 != strcmp(sys_name, unlearn_name)) {
        ApeParFile * unlearn_pfile = 0;
        /* Clone the system parameter file to create the file to be unlearned with values taken from the system file. */
        if (eOK == status) status = ape_io_clone_file(sys_pfile, &unlearn_pfile);

        /* Rename the cloned file to the name of the file to be unlearned. */
        if (eOK == status) status = ape_io_set_file_name(unlearn_pfile, unlearn_name);

        /* Save the cloned file no matter what. */
        if (eOK == status) status = ape_io_write(unlearn_pfile, 1);

        /* Clean up. */
        ape_io_destroy_file(unlearn_pfile);
      }
    }

    /* Clean up without saving current parameter file. */
    ape_trad_close(0);

    /* Go on to the next task in the list, if any. */
    --argc; ++argv;
  }

  return status;
}

int ape_binary_pcheck(int argc, char ** argv) {
  int status = eOK;
  int check_status = eOK;
  char is_file = 0;
  ApeParFile * par_file = 0;

  /* Interpret environment. */
  status = ape_util_interpret_env();

  /* Interpret (and adjust if necessary) the command line to see if a parameter file name was given. */
  if (eOK == status) status = adjust_cmd(&argc, &argv, &is_file);

  /* Check arguments for basic validity. */
  if (eOK == status) status = check_args(argc, argv, 1, -1);

  /* Issue usage message and return if basic conditions were not met. */
  if (eOK != status) {
    usage("pcheck tool-name-or-par-file-name [ arguments ].\n");
    return status;
  }

  if (eOK == status) {
    int idx = 0;
    ape_msg_info(2, "pcheck: checking validity of \"%s", argv[0]);
    for (idx = 1; idx < argc; ++idx) {
      ape_msg_info(2, " %s", argv[idx]);
    }
    ape_msg_info(2, "\"\n");
  }

  /* Show debug messages indicating any/all problems with parameter file. */
  ape_msg_debug_enable(1);

  /* Load parameters from file(s), applying command line arguments as needed. */
  if (eOK == status) status = read_par_file(argc, argv, is_file, &par_file);

  /* Check the content of the parameter file to ensure its format is valid,
     including all parameter values. */
  check_status = ape_io_check_file_format(par_file, 1);
  status = eOK == status ? check_status : status;

  if (0 == is_file) ape_trad_close(0); /* Don't write while closing. */
  else ape_io_destroy_file(par_file);

  return status;
}

int ape_binary_performance(int argc, char ** argv) {
  int status = eOK;
  size_t ii = 0;
  const size_t num_run = 500;
  char ** new_argv = (char **) calloc(++argc, sizeof(char *));
  char * ape_test = "ape_test";

  /* Set up arguments. */
  new_argv[0] = ape_test;
  for (ii = 1; ii != argc; ++ii) new_argv[ii] = argv[ii - 1];

  for (ii = 0; eOK == status && ii != num_run; ++ii) {
    /* Perform a bunch of ape_trad_inits to see how fast they run. */
    status = ape_trad_init(argc, new_argv);
    /* Close the parameter file each time, discarding changes. */
    ape_trad_close(0);
  }
  if (eOK != status) ape_msg_error("On trial %d, ape_trad_init returned status %d.\n", ii, status);

  free(new_argv); new_argv = 0;

  return status;
}

void ape_binary_test(void) {
  char * orig_pfiles = 0;
  int status = eOK;

  /* Interpret environment. */
  status = ape_util_interpret_env();

  if (eOK == status) status = ape_io_get_pfiles(&orig_pfiles);

  if (eOK == status) {
    /* Set pfiles in such a way that values will not be saved. */
    status = ape_io_set_pfiles(QS".");
  }

  if (eOK != status) {
    ape_test_failed("ape_binary_test was unable to set up to test ape binaries. (Status was %d.)\n", status);
  }

  /* Test pget. */
  if (eOK == status) {
    ApeBinary binary;
    binary.func = &ape_binary_pget;
    binary.name = "pget";
    binary.usage = "binary-or-par-file-name [ par1 par2 par3 ... ]";
    binary.min_argc = 1;
    binary.max_argc = -1; /* No maximum. */

    /* Test displaying usage. */
    { char * argv[] = { "pget", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pget\") returned %d, not %d as expected.\n", status, expected_status);
    }

    /* Test error for unknown parameter file. */
    { char * argv[] = { "pget", "fooble", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pget fooble\") returned %d, not %d as expected.\n", status, expected_status);
    }

    /* Test error for unknown parameter. */
    { char * argv[] = { "pget", "ape_test.par", "fooble", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pget ape_test.par fooble\") returned %d, not %d as expected.\n",
          status, expected_status);
    }

    /* Test success case. */
    { char * argv[] = { "pget", "ape_test.par", "dhvalid", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = eOK;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pget ape_test.par dhvalid\") returned %d, not %d as expected.\n",
          status, expected_status);
    }
  }

  /* Test plist. */
  if (eOK == status) {
    ApeBinary binary;
    binary.func = &ape_binary_plist;
    binary.name = "plist";
    binary.usage = "binary-or-par-file-name";
    binary.min_argc = 1;
    binary.max_argc = 1;

    /* Test displaying usage. */
    { char * argv[] = { "plist", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"plist\") returned %d, not %d as expected.\n", status, expected_status);
    }

    /* Test error for unknown parameter file. */
    { char * argv[] = { "plist", "fooble", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"plist fooble\") returned %d, not %d as expected.\n", status, expected_status);
    }

    /* Test error for unknown parameter. */
    { char * argv[] = { "plist", "ape_test.par", "fooble", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"plist ape_test.par fooble\") returned %d, not %d as expected.\n",
          status, expected_status);
    }

    /* Test success cases. */
    { char * argv[] = { "plist", "ape_test.par", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = eOK;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"plist ape_test.par\") returned %d, not %d as expected.\n",
          status, expected_status);
    }

    { char * argv[] = { "plist", "ape_test_no_mode.par", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int expected_status = eOK;
      int status = eOK;
      FILE * fp = fopen("ape_test_no_mode.par", "w");
      if (0 != fp) {
        fprintf(fp, "sql, s, ql, , , , \"sql parameter, no mode parameter should follow\"\n");
        fclose(fp); fp = 0;
      } else {
        ape_test_failed("failed to set up to run ape_binary_run(\"plist ape_test_no_mode.par\").\n",
          status, expected_status);
      }
      status = ape_binary_run(argc, argv, &binary);
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"plist ape_test_no_mode.par\") returned %d, not %d as expected.\n",
          status, expected_status);
      remove("ape_test_no_mode.par");
    }
  }

  /* Test pquery. */
  if (eOK == status) {
    ApeBinary binary;
    binary.func = &ape_binary_pquery;
    binary.name = "pquery";
    binary.usage = "binary-or-par-file-name parameter-name [ cmd-line-args ]";
    binary.min_argc = 2;
    binary.max_argc = -1;

    /* Test displaying usage. */
    { char * argv[] = { "pquery", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pquery\") returned %d, not %d as expected.\n", status, expected_status);
    }
    { char * argv[] = { "pquery", "ape_test", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pquery ape_test\") returned %d, not %d as expected.\n", status, expected_status);
    }

    /* Test error for unknown parameter file. */
    { char * argv[] = { "pquery", "fooble", "sh", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pquery fooble sh\") returned %d, not %d as expected.\n", status, expected_status);
    }

    /* Test error for unknown parameter. */
    { char * argv[] = { "pquery", "ape_test.par", "fooble", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = 1;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pquery ape_test.par fooble\") returned %d, not %d as expected.\n",
          status, expected_status);
    }

    /* Test success case. */
    { char * argv[] = { "pquery", "ape_test.par", "sh", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_run(argc, argv, &binary);
      int expected_status = eOK;
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pquery ape_test.par sh\") returned %d, not %d as expected.\n",
          status, expected_status);
    }

    /* Confirm that pquery rewrites parameter file with modified learned parameters, but unmodified hidden parameters. */
    { char * argv[] = { "pquery", "ape_test.par", "sql", "sql=pquery-learned", "sh=pquery-should-not-learn", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = eOK;
      int expected_status = eOK;
      char * sh_init_value = 0;
      ApeParFile * par_file = 0;
      char * orig_pfiles = 0;

      /* Test with PFILES set to just the current directory, to make sure this works even if only one par file is found. */
      ape_io_get_pfiles(&orig_pfiles);
      ape_io_set_pfiles(".");

      /* Read the parameter file. */
      status = read_par_file(1, argv + 1, 1, &par_file);
      if (eOK == status) {
        ApeListIterator itor = 0;
        /* Fetch initial values of sh parameter, that may be wrongly affected by this test. */
        status = ape_io_find_par("sh", par_file, &itor);
        if (eOK == status) {
          ApePar * par = (ApePar *) ape_list_get(itor);
          ape_par_get_string(par, &sh_init_value);
        }
        itor = 0;
      }

      /* Clean up. */
      ape_io_destroy_file(par_file); par_file = 0;

      status = ape_binary_run(argc, argv, &binary);
      if (expected_status != status)
        ape_test_failed("ape_binary_run(\"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\") "
          "returned %d, not %d as expected.\n", status, expected_status);
      status = read_par_file(1, argv + 1, 1, &par_file);
      if (expected_status != status) {
        ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
          "could not read parameter file.\n", status, expected_status);
      } else {
        ApeListIterator itor = 0;
        status = ape_io_find_par("sql", par_file, &itor);
        if (expected_status != status) {
          ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
            "could not find parameter \"sql\".\n", status, expected_status);
        } else {
          ApePar * par = (ApePar *) ape_list_get(itor);
          char * value = 0;
          status = ape_par_get_string(par, &value);
          if (expected_status != status) {
            ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
              "could not get value of parameter \"sql\".\n", status, expected_status);
          } else if (0 != strcmp(value, "pquery-learned")) {
            ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
              "value of \"sql\" was not \"pquery-learned\".\n", status, expected_status);
          }
          free(value); value = 0;
        }
        itor = 0;
        status = ape_io_find_par("sh", par_file, &itor);
        if (expected_status != status) {
          ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
            "could not find parameter \"sh\".\n", status, expected_status);
        } else {
          ApePar * par = (ApePar *) ape_list_get(itor);
          char * value = 0;
          status = ape_par_get_string(par, &value);
          if (expected_status != status) {
            ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
              "could not get value of parameter \"sh\".\n", status, expected_status);
          } else if (0 != strcmp(value, sh_init_value)) {
            ape_test_failed("after \"pquery ape_test.par sql sql=pquery-learned sh=pquery-should-not-learn\" "
              "value of \"sh\" was incorrect.\n", status, expected_status);
          }
          free(value); value = 0;
        }
      }
      ape_io_destroy_file(par_file); par_file = 0;
      free(sh_init_value); sh_init_value = 0;
      ape_io_set_pfiles(orig_pfiles);
      free(orig_pfiles); orig_pfiles = 0;
    }
  }
  /* Test pset. */
  if (eOK == status) {
    /* Test using name of binary only. */
    { char * argv[] = { "ape_test", "pset-string", "sh=pset-string", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_pset(argc, argv);
      int expected_status = eOK;
      if (expected_status != status)
        ape_test_failed("ape_binary_pset(\"ape_test.par pset-string sh=pset-string\") returned %d, not %d as expected.\n",
          status, expected_status);
    }
    /* Test using name of parameter file. */
    { char * argv[] = { "ape_test.par", "pset-string", "sh=pset-string", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_pset(argc, argv);
      int expected_status = eOK;
      if (expected_status != status)
        ape_test_failed("ape_binary_pset(\"ape_test.par pset-string sh=pset-string\") returned %d, not %d as expected.\n",
          status, expected_status);
    }
  }
  /* Test punlearn. */
  if (eOK == status) {
    const char new_pfiles[] = "pfiles"QS".";
    ape_io_set_pfiles(new_pfiles);
    /* Test unlearning ape_test.par. */
    { char * argv[] = { "ape_test", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_punlearn(argc, argv);
      int expected_status = eOK;
      if (expected_status != status)
        ape_test_failed("ape_binary_punlearn(\"ape_test\") with PFILES == \"%s\" returned %d, not %d as expected.\n",
          new_pfiles, status, expected_status);
    }
    /* Test unlearning invalid par file. */
    { char * argv[] = { "non-existent", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_punlearn(argc, argv);
      int expected_status = eUninitialized;
      if (expected_status != status)
        ape_test_failed("ape_binary_punlearn(\"non-existent\") with PFILES == \"%s\" returned %d, not %d as expected.\n",
          new_pfiles, status, expected_status);
    }
  }
  /* Test pcheck. */
  if (eOK == status) {
    const char new_pfiles[] = "pfiles"QS".";
    ape_io_set_pfiles(new_pfiles);
    /* Check ape_test.par. */
    { char * argv[] = { "ape_test", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_pcheck(argc, argv);
      int expected_status = eValueAboveMax;
      if (expected_status != status)
        ape_test_failed("ape_binary_pcheck(\"ape_test\") with PFILES == \"%s\" returned %d, not %d as expected.\n",
          new_pfiles, status, expected_status);
    }
  }

  /* Test pcheck for a file with auto parameters but no mode parameter. */
  if (eOK == status) {
    const char new_pfiles[] = "pfiles"QS".";
    FILE * fp = fopen("ape_test_no_mode.par", "w");
    int expected_status = eOK;
    if (0 != fp) {
      fprintf(fp, "sa, s, a, , , , \"sa parameter, no mode parameter, should be an error\"\n");
      fclose(fp); fp = 0;
    } else {
      ape_test_failed("failed to set up to run ape_binary_run(\"pcheck ape_test_no_mode.par\").\n",
        status, expected_status);
    }
    ape_io_set_pfiles(new_pfiles);
    /* Check ape_test.par. */
    { char * argv[] = { "ape_test_no_mode", 0 };
      int argc = sizeof(argv) / sizeof(argv[0]) - 1;
      int status = ape_binary_pcheck(argc, argv);
      int expected_status = eNoModePar;
      if (expected_status != status)
        ape_test_failed("ape_binary_pcheck(\"ape_test_no_mode\") with PFILES == \"%s\" returned %d, not %d as expected.\n",
          new_pfiles, status, expected_status);
    }
    remove("ape_test_no_mode.par");
  }

  /* Test adjust_cmd. */
  { int status = eOK;
    int argc = 1;
    char ** argv = (char **) calloc(2, sizeof(char *));
    char is_file = 1;

    if (0 == argv) {
      status = eDynAllocFailed;
      ape_test_failed("Unable to set up to test adjust_cmd in ape_io.\n");
    }
    
    if (eOK == status) {
      argv[0] = "filename.par";
      status = adjust_cmd(&argc, &argv, &is_file);
      ape_test_cmp_long("adjust_cmd(\"filename.par\")", is_file, 0, status, eOK);
    }
    free(argv); argv = 0;
  }
  { int status = eOK;
    int argc = 1;
    char ** argv = (char **) calloc(2, sizeof(char *));
    char is_file = 0;

    if (0 == argv) {
      status = eDynAllocFailed;
      ape_test_failed("Unable to set up to test adjust_cmd in ape_io.\n");
    }
    
    if (eOK == status) {
      argv[0] = "./filename.par";
      status = adjust_cmd(&argc, &argv, &is_file);
      ape_test_cmp_long("adjust_cmd(\"./filename.par\")", is_file, 1, status, eOK);
    }
    free(argv); argv = 0;
  }

  /* Clean up. */
  ape_io_set_pfiles(orig_pfiles);
  free(orig_pfiles); orig_pfiles = 0;
}


#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_binary.c,v $
 * Revision 1.30  2007/11/16 18:16:54  peachey
 * Fix a host of issues related to error checking and reporting. Only
 * pcheck now calls ape_io_check_file_format to check the parameter
 * file for validity and report errors. Removed unneccessary and
 * redundant checking from the pcheck code and streamlined pquery2
 * to avoid a duplicated error message.
 *
 * Revision 1.29  2007/11/12 19:46:43  peachey
 * In all the p* binaries, call ape_util_interpret_env before doing much else.
 *
 * Revision 1.28  2007/11/12 16:55:51  peachey
 * Use implicit default mode of ql when mode parameter is
 * missing in plist and pquery2. In pcheck, report this error. Test all the above.
 *
 * Revision 1.27  2007/10/09 16:45:48  peachey
 * Use ape_par_redirect_prompt_stream to handle prompt redirection. Clean up
 * logic in readline/non-readline/windows cases.
 *
 * Revision 1.26  2007/07/26 16:15:40  peachey
 * Add ape_binary_performance, for measuring Ape's speed.
 *
 * Revision 1.25  2007/07/26 16:11:01  peachey
 * Remove block that defines USE_READLINE, since that #define is now
 * provided by the Makefile.
 *
 * Revision 1.24  2006/11/30 16:40:45  peachey
 * Fix bug in test code: pfiles needs to be set correctly for Windows test.
 *
 * Revision 1.23  2006/11/24 22:25:28  peachey
 * Correct a typo in a unit test which prevented it from working as intended.
 *
 * Revision 1.22  2006/11/24 20:28:49  peachey
 * Modify pquery2 so that it reverts hidden parameters even in situations
 * in which there is no system parameter file.
 *
 * Revision 1.21  2006/11/03 21:32:17  peachey
 * If no system parameter file is available for unlearning, just skip that step.
 *
 * Revision 1.20  2006/10/30 17:56:54  peachey
 * Fix a bug: pquery2 was behaving like pquery, i.e. learning
 * all parameters, hidden or not. Changed to unlearn hidden parameters
 * before saving the file.
 *
 * Revision 1.19  2006/10/19 17:57:58  peachey
 * Save parameter file in pquery binary.
 *
 * Revision 1.18  2006/08/22 21:00:53  peachey
 * Rationalize pcheck behavior along the lines of the other utilities.
 *
 * Revision 1.17  2006/08/22 20:37:13  peachey
 * Move body of pcheck binary to a function ape_binary_pcheck. Close prompt
 * stream after redirecting it to prevent unneeded memory usage.
 *
 * Revision 1.16  2006/07/05 19:38:26  peachey
 * Redirect prompts to the terminal for use with backticks in perl scripts.
 *
 * Revision 1.15  2006/06/23 02:49:12  peachey
 * Require that an argument contain slash (and/or backslash on Windows)
 * before concluding it is a file (in addition to ending in .par.
 *
 * Revision 1.14  2006/06/20 19:02:52  peachey
 * In plist, include minimum of one space between parameter value and prompt.
 *
 * Revision 1.13  2006/06/20 17:22:03  peachey
 * Fix unit test to agree with change in error status caused by last change.
 *
 * Revision 1.12  2006/06/20 03:07:01  peachey
 * Be more tolerant of errors in the local parameter file. The point of
 * punlearn is after all to fix problems.
 *
 * Revision 1.11  2006/06/16 01:18:22  peachey
 * Check file format when argument is interpreted as a file name.
 *
 * Revision 1.10  2006/06/13 14:49:32  peachey
 * Tweak test of punlearn so that it works in HEAdas context. In HEAdas,
 * ape_test binary is right in the test directory, so previous version of the
 * test didn't work.
 *
 * Revision 1.9  2006/06/09 04:10:02  peachey
 * Rationalize pset and punlearn.
 *
 * Revision 1.8  2006/06/08 02:24:18  peachey
 * Rationalize pget, plist and pquery2 binaries to handle par files
 * the same way, or to use binary names and PFILES to find the parameter file(s).
 *
 * Revision 1.7  2006/06/07 15:49:43  peachey
 * Rework ape_binary_pget to handle correctly all possibilities for PFILES
 * and/or command lines.
 *
 * Revision 1.6  2006/05/31 01:36:51  peachey
 * Rename ape_par_get_string to ape_par_get_field and ape_par_get_string_array to
 * ape_par_get_field_array.
 *
 * Revision 1.5  2006/05/23 19:31:12  peachey
 * Add punlearn facilities through ape_binary_punlearn function.
 *
 * Revision 1.4  2006/05/23 16:23:36  peachey
 * Add pset facility through ape_binary_pset function.
 *
 * Revision 1.3  2006/05/22 17:36:17  peachey
 * Add pquery functionality.
 *
 * Revision 1.2  2006/05/22 01:16:58  peachey
 * Make plist output look right.
 *
 * Revision 1.1  2006/05/19 17:42:21  peachey
 * Add ape_binary module.
 *
*/
