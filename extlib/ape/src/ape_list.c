/** \file ape_list.h
    \brief Implementation of list facility.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_list.h"
#include "ape/ape_msg.h"
#include "ape/ape_test.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ApeList {
  ApeListEntry * end;
  size_t size;
};

struct ApeListEntry {
  void * data;
  ApeListEntry * next;
  ApeListEntry * prev;
};

ApeList * ape_list_create(void) {
  /* Create a list with no elements. */
  ApeList * ape_list = (ApeList *) malloc(sizeof(ApeList));
  if (0 != ape_list) {
    ApeListEntry * itor = 0;

    /* Initialize members. */
    ape_list->end = 0;
    ape_list->size = 0;

    /* Create terminal for list. */
    itor = (ApeListEntry *) malloc(sizeof(ApeListEntry));
    if (0 != itor) {
      /* Terminal has no data, and points to itself on either side. */
      itor->data = 0;
      itor->next = itor;
      itor->prev = itor;
      ape_list->end = itor;
    } else {
      ape_msg_debug("ape_list_create: %s.\n", "allocation failed for new ApeListEntry");
      /* Could not create terminal, so creating a list must fail in a clean way. */
      free(ape_list);
      ape_list = 0;
    }
  } else {
    ape_msg_debug("ape_list_create: %s.\n", "allocation failed for new ApeList");
  }
  return ape_list;
}

void ape_list_destroy(ApeList * ape_list) {
  if (0 != ape_list) {
    ApeListEntry * begin = ape_list_begin(ape_list);
    ApeListEntry * itor = 0;
    ApeListEntry * curr = 0;

    for (itor = ape_list_end(ape_list); itor != begin;) {
      /* Save current entry. */
      curr = itor;

      /* Follow link to previous entry. */
      itor = ape_list_prev(curr);

      /* Destroy current entry. */
      curr->prev = 0;
      curr->next = 0;
      curr->data = 0;
      free(curr);
    }

    /* Finally, destroy 0th entry. */
    begin->prev = 0;
    begin->next = 0;
    begin->data = 0;
    free(begin);

    /* To assist in spotting memory problems, zero-out the list. */
    ape_list->size = 0;
    ape_list->end = 0;

    /* Finally, free the list. */
    free(ape_list);
  }
}

ApeListIterator ape_list_insert(ApeList * ape_list, ApeListIterator position, void * data) {
  ApeListEntry * itor = 0;
  if (0 != ape_list) {
    /* Add at the end by default. */
    ApeListEntry * end = ape_list_end(ape_list);
    ApeListEntry * found = end;

    /* If caller gave a position hint as to where to add the new element, use it. */
    if (0 != position && position != found) {
      /* Find the given position in this list. */
      for (itor = ape_list_begin(ape_list); itor != end; itor = ape_list_next(itor)) {
        if (itor == position) {
          found = itor;
          break;
        }
      }
      if (end == found) {
        ape_msg_debug("ape_list_insert: %s.\n", "caller's position argument was not found in list");
      }
    }

    /* Create new entry at the found position. (Reuse itor for new entry). */
    itor = (ApeListEntry *) malloc(sizeof(ApeListEntry));
    if (0 != itor) {
      /* Place new entry just before the "found" position. */
      itor->data = data;
      itor->next = found;
      itor->prev = found->prev;

      /* Connect up next and previous elements. It is not necessary to check that they exist; the terminal
         guarantees they do. */
      itor->next->prev = itor;
      itor->prev->next = itor;

      /* Increment the number of items in the list. */
      ++ape_list->size;
    }
  } else {
    ape_msg_debug("ape_list_insert: %s.\n", "NULL passed for first argument \"ape_list\"");
  }
  return itor;
}

ApeListIterator ape_list_append(ApeList * ape_list, void * data) {
  return ape_list_insert(ape_list, ape_list_end(ape_list), data);
}

void ape_list_erase(ApeList * ape_list, void * data) {
  if (0 != ape_list) {
    ApeListEntry * itor = 0;
    for (itor = ape_list_begin(ape_list); itor != ape_list_end(ape_list); itor = ape_list_next(itor)) {
      /* Look for entries with matching data. */
      if (itor->data == data) {
        /* Store the previous position so that we can continue iterating after removing this entry.
           Note that the iterator is incremented at the end of the loop, which is why we store the previous position. */
        ApeListEntry * prev = ape_list_prev(itor);
        /* Remove the entry. */
        ape_list_remove_entry(ape_list, itor);
        /* Go back to previous position. This is OK even at the beginning of the list, because the first element
           connects back to the terminal. */
        itor = prev;
      }
    }
  } else {
    ape_msg_debug("ape_list_erase: %s.\n", "NULL passed for first argument \"ape_list\"");
  }
}

ApeListIterator ape_list_remove_entry(ApeList * ape_list, ApeListIterator position) {
  ApeListEntry * next = 0;
  if (0 != ape_list) {
    ApeListEntry * itor = 0;
    /* If position is not found, return end of the list. */
    ApeListEntry * end = ape_list_end(ape_list);
    next = end;

    if (0 != position) {
      /* Find the given position in the list, if it is present. */
      for (itor = ape_list_begin(ape_list); itor != end; itor = ape_list_next(itor)) {
        if (itor == position) {
          /* Found this iterator, so save the next iterator for return. */
          next = position->next;

          /* Connect up next and previous elements. It is not necessary to check that they exist; the terminal
             guarantees they do. */
          position->prev->next = position->next;
          position->next->prev = position->prev;

          /* To assist in spotting memory problems, zero-out the entry. */
          position->prev = 0;
          position->next = 0;
          position->data = 0;
          free(position);

          /* Decrement the number of items in the list. */
          --ape_list->size;
          break;
        }
      }
    } else {
      ape_msg_debug("ape_list_remove_entry: %s.\n", "NULL passed for second argument \"position\"");
    }
  } else {
    ape_msg_debug("ape_list_remove_entry: %s.\n", "NULL passed for first argument \"ape_list\"");
    if (0 == position)
      ape_msg_debug("ape_list_remove_entry: %s.\n", "NULL passed for second argument \"position\"");
  }
  return next;
}

ApeListIterator ape_list_begin(ApeList * ape_list) {
  ApeListEntry * begin = 0;
  if (0 != ape_list) {
    /* The beginning is the element after the terminal. */
    begin = ape_list->end->next;
  } else {
    ape_msg_debug("ape_list_begin: %s.\n", "NULL passed for argument \"ape_list\"");
  }
  return begin;
}

ApeListIterator ape_list_end(ApeList * ape_list) {
  ApeListEntry * end = 0;
  if (0 != ape_list) {
    /* The beginning is the element after the terminal. */
    end = ape_list->end;
    /* Just in case anyone assigned to the terminal, which is a no-no. */
    end->data = 0;
  } else {
    ape_msg_debug("ape_list_end: %s.\n", "NULL passed for argument \"ape_list\"");
  }
  return end;
}

ApeListIterator ape_list_next(ApeListIterator itor) {
  if (0 != itor) {
    /* Go on to the next element. */
    itor = itor->next;
  } else {
    ape_msg_debug("ape_list_next: %s.\n", "NULL passed for argument \"itor\"");
  }
  return itor;
}

ApeListIterator ape_list_prev(ApeListIterator itor) {
  if (0 != itor) {
    /* Go back to the previous element. */
    itor = itor->prev;
  } else {
    ape_msg_debug("ape_list_prev: %s.\n", "NULL passed for argument \"itor\"");
  }
  return itor;
}

void * ape_list_get(ApeListIterator position) {
  if (0 != position) {
    return position->data;
  } else {
    ape_msg_debug("ape_list_get: %s.\n", "NULL passed for argument \"position\"");
  }
  return 0;
}

void ape_list_set(ApeListIterator position, void * data) {
  if (0 != position) {
    position->data = data;
  } else {
    ape_msg_debug("ape_list_set: %s.\n", "NULL passed for first argument \"position\"");
  }
}

size_t ape_list_get_size(ApeList * ape_list) {
  size_t size = 0;
  if (0 != ape_list) {
    size = ape_list->size;
  } else {
    ape_msg_debug("ape_list_get_size: %s.\n", "NULL passed for argument \"ape_list\"");
  }
  return size;
}

void ape_list_test(void) {
  /* For use below, create a list pointer. */
  ApeList * ape_list = 0;

  /* For use below, create a list entry. */
  ApeListEntry ape_list_entry;

  /* For use below, create a list iterator. */
  ApeListIterator itor = 0;

  /* Fake data to store. */
  double data[] = { 0., 1., 2., 3., 4. };

  const size_t data_size = sizeof(data) / sizeof(double);

  void * data_ptr;

  /* Loop variable. */
  size_t ii;

  /* Variable holding size of list. */
  size_t size;

  ape_msg_debug("ape_list_test: %s.\n", "beginning test");

  /* Verify that operations on a NULL list are handled properly. */
  ape_list_destroy(ape_list);

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_insert(0, 0, data + 4).");
  itor = &ape_list_entry;
  itor = ape_list_insert(0, 0, data + 4);
  if (0 != itor) ape_test_failed("ape_list_insert(0, 0, data + 4) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_erase(0, data + 4).");
  ape_list_erase(0, data + 4);

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_remove_entry(0, 0)");
  itor = ape_list_remove_entry(0, 0);
  if (0 != itor) ape_test_failed("ape_list_remove_entry(0, 0) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_remove_entry(0, &ape_list_entry)");
  itor = ape_list_remove_entry(0, &ape_list_entry);
  if (0 != itor) ape_test_failed("ape_list_remove_entry(0, &ape_list_entry) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_begin(0)");
  itor = ape_list_begin(0);
  if (0 != itor) ape_test_failed("ape_list_begin(0) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_end(0)");
  itor = ape_list_end(0);
  if (0 != itor) ape_test_failed("ape_list_end(0) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_get_size(0)");
  size = ape_list_get_size(0);
  if (0 != size) {
    ape_test_failed("ape_list_get_size(0) returned a size of %u, not 0 as expected.\n", size);
  }

  /* Verify that operations using a NULL iterator are handled properly. */
  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_next(0)");
  itor = ape_list_next(0);
  if (0 != itor) ape_test_failed("ape_list_next(0) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_prev(0)");
  itor = ape_list_prev(0);
  if (0 != itor) ape_test_failed("ape_list_prev(0) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_get(0)");
  data_ptr = ape_list_get(0);
  if (0 != data_ptr) ape_test_failed("ape_list_get(0) did not return a NULL pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_set(0, data + 4)");
  ape_list_set(0, data + 4);

  /* Create a test list. */
  ape_list = ape_list_create();
  if (0 == ape_list) ape_test_failed("ape_list_create returned a NULL pointer.\n");

  /* Confirm it has 0 elements. */
  size = ape_list_get_size(ape_list);
  if (0 != size) {
    ape_test_failed("ape_list_get_size(ape_list) returned a size of %u, not 0, for an empty list.\n", size);
  }

  /* Verify that operations on an empty list are handled properly. */
  for (itor = ape_list_begin(ape_list); itor != ape_list_end(ape_list); itor = ape_list_next(itor)) {
    ape_test_failed("iterating over an empty list entered the for loop.\n");
    break;
  }

  ape_list_erase(ape_list, 0);

  ape_list_erase(ape_list, data + 4);

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n", "about to call ape_list_remove_entry(ape_list, 0) for an empty list");
  itor = ape_list_remove_entry(ape_list, 0);
  if (ape_list_end(ape_list) != itor)
    ape_test_failed("ape_list_remove_entry(ape_list, 0) for an empty list did not return end pointer.\n");

  itor = ape_list_remove_entry(ape_list, ape_list_begin(ape_list));
  if (ape_list_end(ape_list) != itor)
    ape_test_failed("ape_list_remove_entry(ape_list, ape_list_begin(ape_list)) for an empty list did not return end pointer.\n");

  itor = ape_list_remove_entry(ape_list, ape_list_end(ape_list));
  if (ape_list_end(ape_list) != itor)
    ape_test_failed("ape_list_remove_entry(ape_list, ape_list_end(ape_list)) for an empty list did not return end pointer.\n");

  ape_msg_debug("ape_list_test: %s.\n", "");
  ape_msg_debug("ape_list_test: %s.\n",
    "about to call ape_list_remove_entry(ape_list, &ape_list_entry) for an entry not in the list");
  itor = ape_list_remove_entry(ape_list, &ape_list_entry);
  if (ape_list_end(ape_list) != itor)
    ape_test_failed("ape_list_remove_entry(ape_list, &ape_list_entry) for an empty list did not return end pointer.\n");

  ape_list_destroy(ape_list);

  /* Create a test list. */
  ape_list = ape_list_create();
  if (0 == ape_list) {
    ape_test_failed("ape_list_create returned a NULL pointer.\n");
    ape_test_failed("ape_list_create is skipping some tests.\n");
  } else {
    /* Populate the list. */
    for (ii = 0; ii < data_size; ++ii) {
      itor = ape_list_insert(ape_list, 0, data + ii);
      if (0 == itor) {
        ape_test_failed("ape_list_insert(ape_list, 0, data + ii) returned a NULL pointer.\n");
      } else if (data + ii != ape_list_get(itor)) {
        ape_test_failed("ape_list_get(itor) did not return expected dataer after valid ape_list_insert.\n");
      } else if (ape_list_begin(ape_list) == ape_list_end(ape_list)) {
        ape_test_failed("after ape_list_insert(ape_list, 0, data + ii) loop, "
          "ape_list's begin and end iterators point to the same entry.\n");
      }
      if (0 != itor) {
        /* Confirm the size. */
        size = ape_list_get_size(ape_list);
        if (ii + 1 != size) {
          ape_test_failed("after ape_list_insert(ape_list, 0, data + ii), "
            "ape_list_get_size(ape_list) returned a size of %u, not %u.\n", size, ii + 1);
        }
      }
    }

    /* Add a couple duplicates of item 3, one at the end, one in front of the other. */
    ape_list_insert(ape_list, 0, data + 3);
    for (ii = 0, itor = ape_list_begin(ape_list); ii < data_size && itor != ape_list_end(ape_list);
      ++ii, itor = ape_list_next(itor)) {
      if (data + 3 == ape_list_get(itor)) ape_list_insert(ape_list, itor, data + 3);
    }

    /* Confirm the size. */
    size = ape_list_get_size(ape_list);
    if (data_size + 2 != size) {
      ape_test_failed("after ape_list_insert(ape_list, 0, data + 3), "
        "ape_list_get_size(ape_list) returned a size of %u, not %u.\n", size, data_size + 2);
    }

    /* Confirm that the contents are as expected. */
    {
      unsigned int expected[] = { 0, 1, 2, 3, 3, 4, 3 };
      size_t num_expected = sizeof(expected) / sizeof(expected[0]);
      for (ii = 0, itor = ape_list_begin(ape_list); ii < num_expected && itor != ape_list_end(ape_list);
        ++ii, itor = ape_list_next(itor)) {
        if (data + expected[ii] != ape_list_get(itor)) {
          ape_test_failed("after ape_list was populated, %dth element was %p, not %p as expected.\n", ii,
            ape_list_get(itor), (void *)(data + expected[ii]));
        }
      }
      if (ii != num_expected) ape_test_failed("after ape_list was populated, too few elements in list.\n");
      if (ape_list_end(ape_list) != itor) ape_test_failed("after ape_list was populated, too many elements in list.\n");
    }

    /* Remove all entries containing data + 3. */
    ape_list_erase(ape_list, data + 3);

    /* Confirm the size. */
    size = ape_list_get_size(ape_list);
    if (data_size - 1 != size) {
      ape_test_failed("after ape_list_erase(ape_list, data + 3), "
        "ape_list_get_size(ape_list) returned a size of %u, not %u.\n", size, data_size - 1);
    }

    /* Confirm that the contents are as expected. */
    {
      int expected[] = { 0, 1, 2, 4 };
      size_t num_expected = sizeof(expected) / sizeof(expected[0]);
      for (ii = 0, itor = ape_list_begin(ape_list); ii < num_expected && itor != ape_list_end(ape_list);
        ++ii, itor = ape_list_next(itor)) {
        if (data + expected[ii] != ape_list_get(itor)) {
          ape_test_failed("after entries with data + 3 were erased, %dth element was %p, not %p as expected.\n", ii,
            ape_list_get(itor), (void *) (data + expected[ii]));
        }
      }
      if (ii != num_expected) ape_test_failed("after entries with data + 3, too few elements in list.\n");
      if (ape_list_end(ape_list) != itor) ape_test_failed("after entries with data + 3, too many elements in list.\n");
    }

    /* Jump to the end of the list. */
    itor = ape_list_end(ape_list);

    /* Move around some. */
    itor = ape_list_prev(itor);
    itor = ape_list_prev(itor);
    itor = ape_list_next(itor);
    itor = ape_list_prev(itor);

    /* Remove the current entry, which should be the third out of four. */
    ape_list_remove_entry(ape_list, itor);

    /* Confirm the size. */
    size = ape_list_get_size(ape_list);
    if (data_size - 2 != size) {
      ape_test_failed("after ape_list_remove_entry(ape_list, itor), "
        "ape_list_get_size(ape_list) returned a size of %u, not %u as expected.\n", size, data_size - 2);
    }

    /* Confirm that the contents are as expected: the element which had data + 2 should have been removed. */
    {
      int expected[] = { 0, 1, 4 };
      size_t num_expected = sizeof(expected) / sizeof(expected[0]);
      for (ii = 0, itor = ape_list_begin(ape_list); ii < num_expected && itor != ape_list_end(ape_list);
        ++ii, itor = ape_list_next(itor)) {
        if (data + expected[ii] != ape_list_get(itor)) {
          ape_test_failed("after third entry was removed, %dth element was %p, not %p as expected.\n", ii,
            ape_list_get(itor), (void *) (data + expected[ii]));
        }
      }
      if (ii != num_expected) ape_test_failed("after third entry was removed, too few elements in list.\n");
      if (ape_list_end(ape_list) != itor) ape_test_failed("after third entry was removed, too many elements in list.\n");
    }

    /* Destroy the test list. */
    ape_list_destroy(ape_list);
  }
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_list.c,v $
 * Revision 1.3  2007/05/15 20:52:48  peachey
 * Additional speed optimizations:
 * o Reduce number of calls to ape_list_end by calling it once and storing
 *   the returned iterator for use in loop logic when it is safe to do so.
 * o Make ape_list_destroy more efficient by iterating through the list
 *   back to front once and freeing the elements instead of repeatedly
 *   calling ape_list_remove_entry.
 *
 * Revision 1.2  2006/06/23 19:00:03  peachey
 * Correct a bug in which sizeof(void*) was used instead of sizeof(array[0]).
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
