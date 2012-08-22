/** \file ape_par_group.c
    \brief Implementation of parameter group facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_par.h"
#include "ape/ape_par_group.h"
#include "ape/ape_test.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ApeParGroup {
  ApeList * ape_list;
};

ApeParGroup * ape_par_group_create(void) {
  ApeParGroup * group = (ApeParGroup *) malloc(sizeof(ApeParGroup));
  if (0 != group) {
    group->ape_list = ape_list_create();
    if (0 == group->ape_list) {
      free(group);
      group = 0;
    }
  }
  return group;
}

void ape_par_group_destroy(ApeParGroup * ape_par_group) {
  if (0 != ape_par_group) {
    ApeList * ape_list = ape_par_group->ape_list;
    if (0 != ape_list) {
      ApeListIterator itor = 0;
      /* Destroy parameters from list in reverse order. */
      for (itor = ape_list_end(ape_list); itor != ape_list_begin(ape_list); ) {
        itor = ape_list_prev(itor);
        ape_par_destroy((ApePar *)(ape_list_get(itor)));
      }
      /* Destroy the list itself. */
      ape_list_destroy(ape_list);
    }
    /* Finally, destroy the group. */
    free(ape_par_group);
  }
}

ApeList * ape_par_group_get_list(ApeParGroup * group) {
  ApeList * ape_list = 0;
  if (0 != group) ape_list = group->ape_list;
  return ape_list;
}

ApePar * ape_par_group_get_par(ApeListIterator position) { return (ApePar *) ape_list_get(position); }

void ape_par_group_test(void) {
  /* Pointer to the list in the test group. */
  ApeList * ape_list = 0;

  /* Create a test group. */
  ApeParGroup * group = ape_par_group_create();
  if (0 == group) ape_test_failed("Could not create ApeParGroup.\n");

  /* Get list from group and confirm it has a list of parameters. */
  ape_list = ape_par_group_get_list(group);
  if (0 == ape_list) ape_test_failed("ApeParGroup's list is null.\n");

  /* Destroy the test group. */
  ape_par_group_destroy(group);
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_par_group.c,v $
 * Revision 1.2  2006/04/10 21:11:41  peachey
 * Fix a small memory leak.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
