/******************************************************************************
 *   File name: headas_error_internal.h                                       *
 *                                                                            *
 * Description: This is not a public header file. It should not be installed. *
 *     It contains private type definitions for internal use by HEAdas error  *
 *     manager code only.                                                     *
 *                                                                            *
 *      Author: James Peachey, LAC, for HEASARC/GSFC/NASA                     *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/
#ifndef HEADAS_ERROR_MANAGER_H
#define HEADAS_ERROR_MANAGER_H

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
/* Include forward declarations for types to be defined here. */
#include "headas_error.h"
/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Type definitions/declarations.                                             *
 ******************************************************************************/

/* Error manager type declarations. */
/******************************************************************************/
struct hd_error_msg_s;
typedef struct hd_error_msg_s hd_error_msg;

struct hd_error_manager_s;
typedef struct hd_error_manager_s hd_error_manager;

/* Definition of hd_error_msg_s structure. */
struct hd_error_msg_s {
  char* Text;
  const char* FileName;
  int LineNumber;
  int ErrorNumber;
};

/* Definition of hd_error_manager_s structure. */
struct hd_error_manager_s {
  hd_error_msg* MsgArrayBegin;
  hd_error_msg* MsgLast;
  hd_error_msg* MsgArrayEnd;
  char* TextArrayBegin;
  char* TextLast;
  char* TextArrayEnd;
  hd_error_map* Map;
  int ErrorNumber;
  int SilenceDump;
};
/******************************************************************************/

/******************************************************************************
 * API for creating and manipulating custom error managers.                   *
 ******************************************************************************/

/* Error manager functions and macros. */
/******************************************************************************/
int HDerror_manager_init(hd_error_manager* manager, hd_error_msg* msgArray,
    unsigned int msgArraySize, char* buf, hd_error_map* map,
    unsigned int bufSize, int status);

int HDerror_manager_throw(hd_error_manager* manager, const char* text,
    const char* fileName, int line, int errNum);

int HDerror_manager_hint(hd_error_manager* manager, const char* text,
    const char* fileName, int line, int errNum);

int HDerror_manager_reset(hd_error_manager* manager);

void HDerror_manager_pop_msg(hd_error_manager* manager, char* stack_empty,
    const char** text, const char** fileName, int* line);

int HDerror_manager_get_err(hd_error_manager* manager);

void HDerror_manager_dump(hd_error_manager* manager, FILE* strm);

void HDerror_manager_find_map_entry(hd_error_manager* manager, int errNum,
    char* found, const char** pkgid, const char** code, const char** text);

int HDerror_manager_get_map(hd_error_manager* manager, hd_error_map** errMap,
    int status);

void HDerror_manager_dump_silence(hd_error_manager* manager, int silent);

int HDerror_manager_dump_is_silent(hd_error_manager* manager);

hd_error_manager* HDerror_manager_get_default(void);
/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif

/******************************************************************************
 * $Log: headas_error_internal.h,v $
 * Revision 1.4  2004/01/07 14:45:52  peachey
 * Add functions HDerror_dump_silence and HDerror_dump_is_silent which
 * set and get a silent flag. When this flag is set, HDerror_dump is a no-op.
 *
 * Revision 1.3  2002/10/09 16:42:22  peachey
 * Make sure that heaerr, heaout, heaprom are used consistently.
 *
 * Revision 1.2  2002/10/04 21:44:39  peachey
 * Replace HDerror_manager_get* family with HDerror_manager_find_map_entry,
 * which gets all useful information from the manager's maps.
 * Add functions to dump the error stack, to get the default text, and to
 * find the map entry corresponding to a given error number.
 *
 * Revision 1.1  2002/09/16 13:45:41  peachey
 * Add headas error handling facilities.
 *
 ******************************************************************************/
