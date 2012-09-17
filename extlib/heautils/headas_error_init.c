/******************************************************************************
 *   File name:                                                               *
 *                                                                            *
 * Description:                                                               *
 *                                                                            *
 *      Author: James Peachey, LAC, for HEASARC/GSFC/NASA                     *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include <stdio.h>
#include <string.h>
#include "headas_error_internal.h"
#include "headas_error_cfitsio.h"
/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Constants.                                                                 *
 ******************************************************************************/
/******************************************************************************/

/******************************************************************************
 * Type definitions.                                                          *
 ******************************************************************************/
/******************************************************************************/

/******************************************************************************
 * Global variable definitions.                                               *
 ******************************************************************************/
/******************************************************************************/

/******************************************************************************
 * Static variable definitions.                                               *
 ******************************************************************************/
/******************************************************************************/

/******************************************************************************
 * Static function definitions.                                               *
 ******************************************************************************/
static const char* cfitsio_msg(hd_error_map* theMap, int errNum) {
  static char msg[FLEN_ERRMSG];

/* The following call gets cfitsio's standard message. */
  ffgerr(errNum, msg);
  return msg;
}

static const char* pil_msg(hd_error_map* theMap, int errNum) {
  static char msg[256];

/* The following call gets pil's standard message. */
  strcpy(msg, PIL_err_handler(errNum));
  return msg;
}
/******************************************************************************/

/******************************************************************************
 * Global function definitions.                                               *
 ******************************************************************************/
int HDerror_init(int status) {
  hd_error_manager* manager;
  hd_error_map* errMap;
  if(HD_OK != status) return status;

  do {
    status = HDerror_map_set_msg_func(&sCFITSIOerrorMap, &cfitsio_msg, status);

    status = HDerror_map_set_pkgid(&sCFITSIOerrorMap, "CFITSIO", status);

    status = HDerror_map_set_msg_func(&sPILerrorMap, &pil_msg, status);

    status = HDerror_map_set_pkgid(&sPILerrorMap, "PIL", status);

    manager = HDerror_manager_get_default();
    
    if(NULL == manager) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    status = HDerror_manager_get_map(manager, &errMap, status);

    status = HDerror_map_push(errMap, &sCFITSIOerrorMap, status);

    status = HDerror_map_push(errMap, &sPILerrorMap, status);

  } while(0);

  return status;
}
/******************************************************************************/

#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log: headas_error_init.c,v $
 * Revision 1.2  2002/10/04 21:41:46  peachey
 * Use new message function facility and package ID facility when setting
 * up to handle PIL and CFITSIO error messages.
 *
 * Revision 1.1  2002/09/16 13:45:41  peachey
 * Add headas error handling facilities.
 *
 ******************************************************************************/
