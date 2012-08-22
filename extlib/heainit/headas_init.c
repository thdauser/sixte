/******************************************************************************
 *   File name: headas_init.c                                                 *
 *                                                                            *
 * Description: Universal initialization code for any HEADAS task. Normally   *
 *     this is called by headas_main.c, and need not be called explicitly.    *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: Mike Tripicco, for HEASARC/GSFC/NASA                          *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headas.h"
#include "headas_error.h"
#include "pil.h"
/******************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Static variable definitions.                                             *
   ****************************************************************************/
static char sPILInitSuccess = 0;
  /****************************************************************************/

  /****************************************************************************
   * Static function definitions.                                             *
   ****************************************************************************/
static int StandardPfile(int argc, char *argv[]);
static int hd_pil_err_logger(char *s);
static void ReportError(int taskStatus);
  /****************************************************************************/

  /****************************************************************************
   * Function definitions.                                                    *
   ****************************************************************************/
/* called by headas_main.c prior to tool subroutine */
int headas_init(int argc, char *argv[]){
    int status = 0;
    int tmpstatus = 0;

    headas_chatpar = -1;
    headas_clobpar = 0;

    /* Register defaults for toolname/version, used in error reporting.
       These functions cannot create an error condition. */
    set_toolname(hdbasename(argv[0]));
    set_toolversion("0.0");

    /* read env vars and setup output streams */
    tmpstatus = HDIO_init();
    if (tmpstatus) HDerror_throw(0, 0, 0, tmpstatus);

    status = status ? status : tmpstatus;

    /* Set up global error handler. This can/should be done even if an
       error occurred with the streams. */
    tmpstatus = HDerror_init(0);
    if (tmpstatus) HDerror_throw(0, 0, 0, tmpstatus);

    status = status ? status : tmpstatus;

    tmpstatus = StandardPfile(argc, argv);

    status = status ? status : tmpstatus;

    return status;
}

/* called by headas_main.c after tool subroutine */
int headas_close(int taskStatus) {
  int status = 0;

  /* Close PIL provided the taskStatus argument was not a PIL error. */
  if(sPILInitSuccess &&
      (PIL_ERR_BASE < taskStatus || PIL_ERR_MIN_IDX > taskStatus)) {
    if(PIL_OK != (status = PILClose(PIL_OK))) {
      HD_ERROR_THROW("PILClose failed", status);
    }
  }

  /* If the taskStatus argument is 0 (no error before this function
     was called) then set it to the local error status. */
  taskStatus = taskStatus ? taskStatus : status;

  /* Report a non-0 overall task error status. */
  if(taskStatus && !HDerror_dump_is_silent()) ReportError(taskStatus);

  return taskStatus;
}

int StandardPfile(int argc, char *argv[]) {
    int status;
    int headas_histpar = -1;

    if (PIL_OK != (status = PILInit(argc, argv))){
	HD_ERROR_THROW("PIL initialization failed", status);
	return status;
    }
    sPILInitSuccess = 1;

/* enable 'batch mode' via HEADASNOQUERY environment variable */
    if (HD_getenv("HEADASNOQUERY")) PILOverrideQueryMode(PIL_QUERY_OVERRIDE);

    status = PILSetLoggerFunction(hd_pil_err_logger);
    if(status != 0){
        HD_ERROR_THROW("PILSetLoggerFunction failed", status);
        return status;
    }

    /* Use heautils's function to check for file existence. */
    status = PILSetFileAccessFunction(&HDfile_check);
    if(status != 0){
        HD_ERROR_THROW("PILSetFileAccessFunction failed", status);
        return status;
    }

    /* unset temporarily to keep missing parameters     */
    /* (chatter/clobber/history) from throwing an error */
    if (PIL_OK != (status = PILSetLoggerFunction(NULL))) {
        HD_ERROR_THROW("PILSetLoggerFunction failed", status);
        return status;
    }

    if (PIL_OK != (status=PILGetInt("chatter", &headas_chatpar))) {
      /* OK not to have chatter parameter     */
      /* but we'll leave headas_chatpar at -1 */
	status = 0;
    }

    status=PILGetBool("clobber", &headas_clobpar);
    if(status != 0){ /* OK not to have clobber parameter */
	status = 0;
    }

    status=PILGetBool("history", &headas_histpar);
    if(status != 0){ /* OK not to have history parameter     */
	             /* but we'll leave headas_histpar at -1 */
	status = 0;
    }

    /* now reset to normal logging mode */
    status = PILSetLoggerFunction(hd_pil_err_logger);
    if(status != 0){
        HD_ERROR_THROW("PILSetLoggerFunction failed", status);
        return status;
    }

    /* Handle history parameter. */
    set_history(headas_histpar);

    return status;
}

static int hd_pil_err_logger(char *s)
{
   fprintf(heaerr, "%s\n", s);
   return(0);
}

static void ReportError(int taskStatus) {
  char task[73] = "";
  char vers[9] = "";
  const char *code;
  const char *pkgid;
  const char *text;
  char found;
  char doDump;
  char cfitsioMsg[FLEN_ERRMSG];

  /* First, dump the Cfitsio stack. */
  doDump = 1;
  while(HD_OK != fits_read_errmsg(cfitsioMsg)) {
    if(doDump) {
      fprintf(heaerr, "Dumping CFITSIO error stack:\n--------------------------------------------------\n");
      doDump = 0; /* Only write the header the first time. */
    }

    fprintf(heaerr, "%s\n", cfitsioMsg);
  }
  /* If doDump is not 1, the while loop above printed something, so
     print a footer. */
  if(0 == doDump) fprintf(heaerr, "--------------------------------------------------\nCFITSIO error stack dump complete.\n");

  /* If HEADAS error stack does not contain an error, throw one now, so
     that HDerror_dump can print the error message. */
  if(0 == HDerror_get()) HDerror_throw(NULL, NULL, 0, taskStatus);

  /* Dump the HEADAS message stack. */
  HDerror_dump(heaerr);

  /* Get error information. */
  HDerror_get_info(taskStatus, &found, &pkgid, &code, &text);

  /* If no info was available, get it from the system. */
  if(!found && 0 < taskStatus) {
    text = strerror(taskStatus);
    if(text) found = 1;
  }

  /* Print error information summary. */
  if(found && (pkgid || code || text)) {
     if(pkgid) fprintf(heaerr, "%s ", pkgid);
     fprintf(heaerr, "ERROR");
     if(code) fprintf(heaerr, " %s", code);
     if(text) fprintf(heaerr, ": %s", text);
     fprintf(heaerr, "\n");
  }

  /* Get information about the task. */
  get_toolname(task);
  get_toolversion(vers);

  /* Print information about the task. This is the final message
     written to heaerr. */
  fprintf(heaerr, "Task");
  if(task) fprintf(heaerr, " %s", task);
  if(vers) fprintf(heaerr, " %s", vers);
  fprintf(heaerr, " terminating with status %d\n", taskStatus);

  /* Make certain the returned status is in the range the shell can
     handle. */
  taskStatus = 201; /* Generic headas return value. */

  /* Give specific codes for known components. */
  if(pkgid) {
    if(strstr(pkgid, "CFITSIO")) taskStatus = 202;
    else if(strstr(pkgid, "PIL")) taskStatus = 203;;
  }
}
  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log: headas_init.c,v $
 * Revision 1.18  2006/06/16 00:21:01  peachey
 * Suppress attempting to print system error messages if error code is negative.
 *
 * Revision 1.17  2006/06/07 07:09:26  peachey
 * Do not call deprecated function PILVerifyCmdLine; PILInit does this.
 *
 * Revision 1.16  2004/08/17 19:30:57  peachey
 * Use new HDfile_check function to check for existence of files.
 *
 * Revision 1.15  2004/01/07 14:46:37  peachey
 * Check whether errors have been silenced before reporting them.
 *
 * Revision 1.14  2003/02/13 20:03:50  peachey
 * Refine error reporting.
 *
 * Revision 1.13  2003/02/13 18:45:55  peachey
 * Perform both error and stream initializations even if an error
 * occurs. This allows meaningful error reporting in such cases.
 *
 * Revision 1.12  2003/02/13 12:25:02  peachey
 * Do not close headas streams; they will be closed automatically on exit.
 * Otherwise, spurious errors can crop up from streams being closed more
 * than once.
 *
 * Revision 1.11  2003/02/12 23:27:45  peachey
 * hdIOInit code moved to headas_stdio.c, and renamed HDIO_init.
 *
 * Revision 1.10  2002/11/06 19:15:26  peachey
 * Initialize heain stream.
 *
 * Revision 1.9  2002/10/29 17:12:57  peachey
 * Do not put PIL error messages on the HEADAS stack; instead write them
 * directly to heaerr stream.
 *
 * Revision 1.8  2002/10/18 19:39:25  peachey
 * Reduce the number of dashes on the CFITSIO error stack banner.
 *
 * Revision 1.7  2002/10/09 16:26:39  peachey
 * Improve error handling in init/close. Use streams heaout, heaerr, heaprom
 * throughout. Move some code out of higher level functions into lower
 * level functions.
 *
 * Revision 1.6  2002/10/04 21:51:03  peachey
 * Changes to use the new error handling facility automatically.
 *
 ******************************************************************************/
