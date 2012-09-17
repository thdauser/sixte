/******************************************************************************
 *   File name: headas_error.c                                                *
 *                                                                            *
 * Description: Implementation of HEAdas basic error handling API functions.  *
 *                                                                            *
 *      Author: James Peachey, LAC, for HEASARC/GSFC/NASA                     *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include "headas_error_internal.h"
/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Constants.                                                                 *
 ******************************************************************************/
enum hd_error_manager_constants {
  HD_ERR_MNGR_ARRAY_SIZE = 32u, /* Default depth of stack. */
  HD_ERR_MNGR_BUF_SIZE = HD_ERR_MNGR_ARRAY_SIZE * 64u /* Arbitrary msg length */
};
/******************************************************************************/

/******************************************************************************
 * Static variable definitions.                                               *
 ******************************************************************************/
/* Declare default error map variables. */

/* Local short-cut for the longer macro used to define static entries. */
#define ERR HD_ERR_MAP_STATIC_ENTRY

/* Static map array of known HEADAS errors to standard messages. */
static hd_error_map_entry sHEADASerrors[] = {
  ERR(HD_ERR_DYN_ALLOC_FAIL, "Attempt to allocate dynamic memory failed."),
  ERR(HD_ERR_MNGR_ERROR, "HEADAS error manager improperly used."),
  ERR(HD_ERR_NULL_POINTER, "Null pointer passed where not expected."),
  ERR(HD_ERR_OUT_OF_BOUNDS, "Attempt to access array element past end of array.")
};

/* Pointer to the static error map. */
static hd_error_map sHEADASerrorMap = HD_ERR_MAP_INIT(sHEADASerrors);

/* Declare default error manager variables. */

/* Message array used by error manager to hold error messages and hints. */
static hd_error_msg sMsgArray[HD_ERR_MNGR_ARRAY_SIZE];

/* Character buffer used by error manager to hold all error text. */
static char sMsgBuf[HD_ERR_MNGR_BUF_SIZE];

/* The HEAdas default error manager. Normally, all error managers must
   be created using HDerror_manager_create. The default error manager is
   an exception, because it must be statically created to guarantee its
   availablity. */
static hd_error_manager sManager_s =
  { sMsgArray, sMsgArray, sMsgArray + HD_ERR_MNGR_ARRAY_SIZE,
    sMsgBuf, sMsgBuf, sMsgBuf + HD_ERR_MNGR_BUF_SIZE, &sHEADASerrorMap, HD_OK };

/* Pointer to the static error manager. */
static hd_error_manager* sManager = &sManager_s;
/******************************************************************************/

/******************************************************************************
 * Global function definitions.                                               *
 ******************************************************************************/
int HDerror_throw(const char* msg, const char* fileName, int line, int errNum) {
  return HDerror_manager_throw(sManager, msg, fileName, line, errNum);
}

int HDerror_hint(const char* msg, const char* fileName, int line, int errNum) {
  return HDerror_manager_hint(sManager, msg, fileName, line, errNum);
}

int HDerror_reset(void) {
  return HDerror_manager_reset(sManager);
}

int HDerror_get(void) {
/* For optimal speed, do not call HDerror_manager_get_err, but instead
   directly return the error number from the default error manager.
  return HDerror_manager_get_err(sManager);
*/
  return sManager->ErrorNumber;
}

void HDerror_dump(FILE* strm) {
  HDerror_manager_dump(sManager, strm);
}

void HDerror_get_info(int errNum, char* found, const char** pkgid,
    const char** code, const char** text) {
  HDerror_manager_find_map_entry(sManager, errNum, found, pkgid, code, text);
}

void HDerror_dump_silence(int silent) {
  HDerror_manager_dump_silence(sManager, silent);
}

int HDerror_dump_is_silent(void) {
  return HDerror_manager_dump_is_silent(sManager);
}

hd_error_manager* HDerror_manager_get_default(void) {
  return sManager;
}
/******************************************************************************/

#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log: headas_error.c,v $
 * Revision 1.4  2004/01/07 14:45:52  peachey
 * Add functions HDerror_dump_silence and HDerror_dump_is_silent which
 * set and get a silent flag. When this flag is set, HDerror_dump is a no-op.
 *
 * Revision 1.3  2002/10/09 16:41:49  peachey
 * Make sure that heaerr, heaout, heaprom are used consistently.
 *
 * Revision 1.2  2002/10/04 21:35:33  peachey
 * Replace HDerror_get_code with a (it is hoped) more useful function,
 * HDerror_get_info, which gets the package id and default text in addtion
 * to the code. Also some cosmetic reorganization of the code.
 *
 * Revision 1.1  2002/09/16 13:45:41  peachey
 * Add headas error handling facilities.
 *
 ******************************************************************************/
