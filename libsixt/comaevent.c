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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include "comaevent.h"

CoMaEvent* getCoMaEvent(int* status){
  CoMaEvent* ce=(CoMaEvent*)malloc(sizeof(CoMaEvent));
  if (NULL==ce) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for CoMaEvent!\n",
		   *status);
    return(ce);
  }

  //Initialization:
  ce->time=0.;
  ce->rawx=0;
  ce->rawy=0;
  ce->charge=0.;

  return(ce);
}

void freeCoMaEvent(CoMaEvent** const ce)
{
  if (NULL!=*ce) {
    free(*ce);
    *ce=NULL;
  }
}
