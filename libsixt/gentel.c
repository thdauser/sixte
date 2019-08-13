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

#include "gentel.h"


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


GenTel* newGenTel(int* const status)
{
  // Allocate memory.
  GenTel* tel=(GenTel*)malloc(sizeof(GenTel));
  if (NULL==tel) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenTel failed");
    return(tel);
  }

  // Initialize all pointers with NULL.
  tel->psf         =NULL;
  tel->vignetting  =NULL;
  tel->arf         =NULL;
  tel->arf_filename=NULL;

  // Set initial values.
  tel->focal_length=0.;
  tel->fov_diameter=0.;
  tel->num_imaged=0;

  return(tel);
}


void destroyGenTel(GenTel** const tel)
{
  if (NULL!=*tel) {
    if (NULL!=(*tel)->arf_filename) {
      free((*tel)->arf_filename);
    }
    freeARF((*tel)->arf);
    destroyPSF(&(*tel)->psf);
    destroyVignetting(&(*tel)->vignetting);

    free(*tel);
    *tel=NULL;
  }
}

void check_if_imaged(const GenTel* const tel) {
  headas_chat(5, "Telescope imaged %ld photons\n", tel->num_imaged);

  if (tel->num_imaged == 0) {
    SIXT_WARNING("No photons imaged by the telescope! Check your pointing or exposure time");
  }

};
