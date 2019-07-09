/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/



#include "namelist.h"

// Dump a list for debugging
static inline void __dump_list(name_list_t *list)
{
  printf("***Dumping list:\n");
  for (int ii=0; ii<list->num; ii++)
  {
    printf("%d: \"%s\"\n", ii, list->names[ii]);
  }
  printf("***END\n");
}

static inline int __find_first_substr(char *str, char tofind)
{
  int pos = 0;

  while ( str[pos] != '\0' && str[pos] != tofind )
  {
    pos++;
  }

  if ( str[pos] == '\0' )
  {
    return -1;
  } else {
    return pos;
  }
}

int get_num_files_in_list(char *list, char delim)
{
  int num = 0;
  while (*list != '\0')
  {
    if ( *(list++) == delim )
    {
      num++;
    }
  }
  return num+1;
}


name_list_t *split_list(const char *str, char delim)
{
  name_list_t *list =(name_list_t *) malloc(sizeof(name_list_t));
  CHECK_MALLOC_RET_NULL(list);
  list->num = get_num_files_in_list((char *) str, delim);
  list->names = (char **) malloc(list->num * sizeof(char *));
  CHECK_MALLOC_RET_NULL(list->names);

  char *search = (char *) str;
  for (int ii=0; ii<list->num; ii++)
  {
    int pos = __find_first_substr(search, delim);
    if (pos > 0)
    {
      list->names[ii] = (char *) malloc((pos+1) * sizeof(char));
      CHECK_MALLOC_RET_NULL(list->names[ii]);
      strncpy(list->names[ii], (const char *) search, pos);
      list->names[ii][pos] = '\0';
      search += pos+1;
    } else {
      list->names[ii] = (char *) malloc(strlen(search)+1);
      CHECK_MALLOC_RET_NULL(list->names[ii]);
      strcpy(list->names[ii], (const char *) search);
    }
  }

  return list;
}


void destroy_name_list(name_list_t *list)
{
  for (int ii=0; ii<list->num; ii++)
  {
    free(list->names[ii]);
  }
  free(list);
}


/*
 * Add two lists. The old input lists will be freed
 */
name_list_t *add_list_to_list(name_list_t *org, name_list_t *append)
{
  name_list_t *ret = (name_list_t *) malloc(sizeof(name_list_t));
  CHECK_MALLOC_RET_NULL(ret);
  ret->num = org->num + append->num;
  ret->names = (char **) malloc(ret->num * sizeof(char *));
  CHECK_MALLOC_RET_NULL(ret->names);

  for(int ii=0; ii<org->num; ii++)
  {
    ret->names[ii] = org->names[ii];
  }
  for (int ii=0; ii<append->num; ii++)
  {
    ret->names[org->num+ii] = append->names[ii];
  }

  free(org->names);
  free(org);
  free(append->names);
  free(append);

  return ret;
}
