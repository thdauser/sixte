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

#ifndef PULSEIMP_H
#define PULSEIMP_H 1

#include "sixt.h"
#include "tesproftemplates.h"

#define TOOLSUB pulsetemplimport_main
#include "headas_main.c"



struct Parameters {
  
  char filename[MAXMSG];
  char ascii[MAXMSG];
  char version[9];
  
  double energy;
  
};

int getpar(struct Parameters* const par);

void read_ascii_pulse(char *ascii, 
		      TESProfilesEntries *prof,
		      double energy,
		      int* const status);

#endif /* PULSEIMP_H */
