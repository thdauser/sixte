/** \file ape_list.h
    \brief Declaration of list facility.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_list_h
#define ape_ape_list_h

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ApeList;
typedef struct ApeList ApeList;

struct ApeListEntry;
typedef struct ApeListEntry ApeListEntry;

typedef ApeListEntry * ApeListIterator;

ApeList * ape_list_create(void);

void ape_list_destroy(ApeList * ape_list);

ApeListIterator ape_list_insert(ApeList * ape_list, ApeListIterator position, void * data);

ApeListIterator ape_list_append(ApeList * ape_list, void * data);

void ape_list_erase(ApeList * ape_list, void * data);

ApeListIterator ape_list_remove_entry(ApeList * ape_list, ApeListIterator position);

ApeListIterator ape_list_begin(ApeList * ape_list);

ApeListIterator ape_list_end(ApeList * ape_list);

ApeListIterator ape_list_next(ApeListIterator itor);

ApeListIterator ape_list_prev(ApeListIterator itor);

void * ape_list_get(ApeListIterator position);

void ape_list_set(ApeListIterator position, void * data);

size_t ape_list_get_size(ApeList * ape_list);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_list.h,v $
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
