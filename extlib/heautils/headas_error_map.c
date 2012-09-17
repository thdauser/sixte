/******************************************************************************
 *   File name: headas_error_map.c                                            *
 *                                                                            *
 * Description: Implementation of the HEAdas error map, which maintains       *
 *     an association between error numbers, symbolic codes, and default      *
 *     error messages. Each map also contains a pointer to other maps,        *
 *     allowing lists of maps to be formed conveniently.                      *
 *                                                                            *
 *      Author: James Peachey, LAC, for HEASARC/GSFC/NASA                     *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include <stdio.h>
#include "headas_error.h"
#include "headas_error_internal.h"
/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Global function definitions.                                               *
 ******************************************************************************/
int HDerror_map_init(hd_error_map* map, hd_error_map_entry* entries,
    unsigned int mapSize, int status) {
  if(HD_OK != status) return status;

  do {
    if(NULL == map || NULL == entries) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    /* Initialize the map. */
    map->EntryArray = entries;
    map->Current = entries;
    map->End = entries + mapSize;
    map->NextMap = NULL;

  } while(0);

  return status;
}

int HDerror_map_set_msg_func(hd_error_map* map,
    const char* (*func)(hd_error_map*, int), int status) {
  if(HD_OK != status) return status;
  do {
    if(NULL == map) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    map->MsgFunc = func;

  } while(0);
  return status;
}

int HDerror_map_set_pkgid(hd_error_map* map, const char* pkgid, int status) {
  if(HD_OK != status) return status;
  do {
    if(NULL == map) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    map->PkgId = pkgid;

  } while(0);
  return status;
}

int HDerror_map_push(hd_error_map* currentMap, hd_error_map* nextMap,
    int status) {
  if(HD_OK != status) return status;

  do {
    if(NULL == currentMap) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    /* Go to the end of the chain. */
    while(NULL != currentMap->NextMap) currentMap = currentMap->NextMap;

    /* Insert the next map at the end. */
    if(NULL != nextMap) nextMap->NextMap = currentMap->NextMap;
    currentMap->NextMap = nextMap;

  } while(0);

  return status;
}

int HDerror_map_add_entry(hd_error_map* map, int errNum, const char* code,
    const char* text, int status) {
  hd_error_map_entry* entry_p;
  if(HD_OK != status) return status;

  do {
    if(NULL == map) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    /* Do not allow duplicates; go through the list looking for an
       entry with the given error number. */
    for(entry_p = map->EntryArray; entry_p < map->Current; ++entry_p) {
      if(errNum == entry_p->Number) break;
    }

    /* If in bounds, supercede this entry with the given information. */
    if(entry_p < map->End) {
      entry_p->Number = errNum;
      entry_p->Code = code;
      entry_p->Text = text;
    } else {
      status = HD_ERR_OUT_OF_BOUNDS;
      continue;
    }

    /* Update current pointer if appropriate. */
    if(entry_p > map->Current) map->Current = entry_p;

  } while(0);

  return status;
}

void HDerror_map_get_entry(hd_error_map* map, int errNum, char* found,
    const char** pkgid, const char** code, const char** text) {
  do {
    hd_error_map_entry* entry_p;

    if(NULL == found) continue;
    *found = 0;

    if(NULL == map) continue;

    /* Go through the list looking for an entry with the given error number. */
    for(entry_p = map->EntryArray; entry_p < map->Current; ++entry_p) {
      if(errNum == entry_p->Number) {
        *found = 1;
        break;
      }
    }

    /* If found in this map, return the requested information. */
    if(*found) {
      if(NULL != pkgid) *pkgid = map->PkgId;
      if(NULL != code) *code = entry_p->Code;
      if(NULL != text) {
        if(NULL != entry_p->Text) *text = entry_p->Text;
        else if(NULL != map->MsgFunc) *text = map->MsgFunc(map, errNum);
        else *text = NULL;
      }
    } else {
      if(NULL != pkgid) *pkgid = NULL;
      if(NULL != code) *code = NULL;
      if(NULL != text) *text = NULL;
    }

  } while(0);
}

hd_error_map* HDerror_map_next_map(hd_error_map* map) {
  /* This deliberately does not check status because it is used
     when reporting status. */
  do {
    if(NULL == map) continue;
    map = map->NextMap;
  } while(0);

  return map;
}
/******************************************************************************/

#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log: headas_error_map.c,v $
 * Revision 1.2  2002/10/04 21:48:51  peachey
 * To support improvements in the error manager, change function
 * HDerror_map_set_msg_func, and replace HDerror_map_get* family
 * with HDerror_map_get_entry, which gets all the information from
 * the map in one call.
 *
 * Revision 1.1  2002/09/16 13:45:41  peachey
 * Add headas error handling facilities.
 *
 ******************************************************************************/
