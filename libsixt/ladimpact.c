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

#include "ladimpact.h"


LADImpact* getLADImpact(int* const status)
{
  LADImpact* imp=(LADImpact*)malloc(sizeof(LADImpact));
  CHECK_NULL_RET(imp, *status,
		 "memory allocation for LADImpact failed", imp);

  // Initalize.
  imp->panel  =0;
  imp->module =0;
  imp->element=0;
  imp->position.x=0.;
  imp->position.y=0.;
  imp->energy =0.;
  imp->time   =0.;
  imp->ph_id  =0;
  imp->src_id =0;

  return(imp);
}


void freeLADImpact(LADImpact** const impact)
{
  if (NULL!=*impact) {
    free(*impact);
    *impact=NULL;
  }
}
