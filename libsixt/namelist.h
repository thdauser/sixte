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


/*
 * Parameter parsing for more convient access to e.g. reading multiple simput files
 */

#include "sixt.h"


// Delimiter for list of multiple filenames
#define DELIM ','

typedef struct {
  char **names;
  int num;
} name_list_t;

/*
 * For a string containing a list of names delimited by delim get the number
 * of these names in the string
 */
int get_num_files_in_list(char *list, char delim);

/*
 * For a string containing a list of names delimited by delim split the string
 * apart and store the individual strings in a list and return it
 */
name_list_t *split_list(const char *str, char delim);

/*
 * Merge to lists into a new one
 * The old lists will be freed!
 */
name_list_t *add_list_to_list(name_list_t *org, name_list_t *append);
