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


   Copyright 2014 Jelle de Plaa, SRON, Thorsten Brand, FAU
*/

#ifndef PULSEGEN_H
#define PULSEGEN_H 1

#include "sixt.h"
#include "tesproftemplates.h"

#define TOOLSUB pulsetemplgen_main
#include "headas_main.c"



struct Parameters {
  int nver;
  char filename[MAXMSG];
  char clobber;
  char history;
  
  double energy_high;
  double energy_low;
  int energy_steps;
  
  TESTemplateInput pinp;
};

void freePar(struct Parameters *par);
int getpar(struct Parameters* const par);

#endif /* PULSEGEN_H */
