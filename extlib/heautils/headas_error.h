/******************************************************************************
 *   File name: headas_error.h                                                *
 *                                                                            *
 * Description: Public header file with declarations needed to use the        *
 *     HEAdas C-based error handling library functions.                       *
 *                                                                            *
 *      Author: James Peachey, LAC, for HEASARC/GSFC/NASA                     *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/
#ifndef HEADAS_ERROR_H
#define HEADAS_ERROR_H

#include <stdio.h>
#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Constants.                                                                 *
 ******************************************************************************/
/* Enumerated error codes for HEAdas standard errors. */
typedef enum hd_error_code {
  HD_OK = 0,
  HD_ERR_BASE = 1000,
  HD_ERR_DYN_ALLOC_FAIL,
  HD_ERR_PROGRAMMER_ERRORS = HD_ERR_BASE + 500,
  HD_ERR_MNGR_ERROR,
  HD_ERR_NULL_POINTER,
  HD_ERR_OUT_OF_BOUNDS
} hd_error_code;

/* Enumerated general purpose constants. */
enum hd_error_constants {
  HD_NO_LINE_NUMBER = -1
};
/******************************************************************************/

/******************************************************************************
 * API for basic error handling.                                              *
 ******************************************************************************/
int HDerror_init(int status);

int HDerror_throw(const char* msg, const char* fileName, int line, int errNum);

int HDerror_hint(const char* msg, const char* fileName, int line, int errNum);

int HDerror_reset(void);

int HDerror_get(void);

void HDerror_dump(FILE* strm);

void HDerror_get_info(int errNum, char* found, const char** pkgid,
    const char** code, const char** text);

void HDerror_dump_silence(int silent);

int HDerror_dump_is_silent(void);

#define HD_ERROR_GET HDerror_get

#define HD_ERROR_SET(STATUS) \
  HDerror_throw((void*) 0, __FILE__, __LINE__, STATUS)

#define HD_ERROR_THROW(MSG, STATUS) \
  HDerror_throw(MSG, __FILE__, __LINE__, STATUS)

#define HD_ERROR_HINT(MSG, STATUS) \
  HDerror_hint(MSG, __FILE__, __LINE__, STATUS)
/******************************************************************************/

/******************************************************************************
 * API for creating and manipulating error message maps.                      *
 ******************************************************************************/

/* Error message map type definitions. */
/******************************************************************************/
struct hd_error_map_entry_s;
typedef struct hd_error_map_entry_s hd_error_map_entry;

struct hd_error_map_s;
typedef struct hd_error_map_s hd_error_map;

struct hd_error_map_entry_s {
  int Number;
  const char* Code;
  const char* Text;
};

struct hd_error_map_s {
  hd_error_map_entry* EntryArray;
  hd_error_map_entry* Current;
  hd_error_map_entry* End;
  hd_error_map* NextMap;
  const char* PkgId;
  const char* (*MsgFunc)(hd_error_map*, int);
};
/******************************************************************************/

/* Error message map functions and macros. */
/******************************************************************************/
int HDerror_map_init(hd_error_map* map, hd_error_map_entry* entries,
    unsigned int mapSize, int status);

int HDerror_map_set_msg_func(hd_error_map* map,
    const char* (*func)(hd_error_map*, int), int status);

int HDerror_map_set_pkgid(hd_error_map* map, const char* pkgid, int status);

int HDerror_map_push(hd_error_map* currentMap, hd_error_map* nextMap,
    int status);

int HDerror_map_add_entry(hd_error_map* map, int errNum, const char* code,
    const char* text, int status);

void HDerror_map_get_entry(hd_error_map* map, int errNum, char* found,
    const char** pkgid, const char** code, const char** text);

hd_error_map* HDerror_map_next_map(hd_error_map* map);

#define HD_ERR_MAP_STATIC_ENTRY(A, B) { A, #A, B }
#define HD_ERR_MAP_INIT(A) { A, A + sizeof(A)/sizeof(hd_error_map_entry), \
    A + sizeof(A)/sizeof(hd_error_map_entry), (void*) 0, (void*) 0, (void*) 0}
/******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif

/******************************************************************************
 * $Log: headas_error.h,v $
 * Revision 1.4  2004/01/07 14:45:52  peachey
 * Add functions HDerror_dump_silence and HDerror_dump_is_silent which
 * set and get a silent flag. When this flag is set, HDerror_dump is a no-op.
 *
 * Revision 1.3  2002/10/09 16:41:58  peachey
 * Make sure that heaerr, heaout, heaprom are used consistently.
 *
 * Revision 1.2  2002/10/04 21:40:47  peachey
 * 1) Replace HDerror_get_code with a (it is hoped) more useful function,
 * HDerror_get_info, which gets the package id and default text in addtion
 * to the code. Also some cosmetic reorganization of the code.
 *
 * 2) Remove the clunky "LookUp" construct from hd_error_map_s. The
 * lookup is now handled in a more natural way by passing a pointer
 * to the map in the MsgFunc construct.
 *
 * 3) Replace HDerror_map_get* family with HDerror_map_get_entry, which
 * gets all the parts of a map entry one might want, while allowing
 * NULL for ignored return arguments.
 *
 * 4) Add PkgId and functions to set/get it from error maps.
 *
 * Revision 1.1  2002/09/16 13:45:41  peachey
 * Add headas error handling facilities.
 *
 ******************************************************************************/
