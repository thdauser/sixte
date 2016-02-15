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


   Copyright 2015 Thorsten Brand, FAU
*/

#ifndef XML2SVG_H
#define XML2SVG_H 1

#include "sixt.h"
#include "obj2d.h"
#include "detstruct2obj2d.h"
#include "sixtesvg.h"
#include "geninst.h"
#include "advdet.h"

#define TOOLSUB xml2svg_main

#include "headas_main.c"

struct Parameters {
  char *XMLFiles;
  char *SVGName;
  double svgwidth;
  double border;
  int drawn;
  char writeid;
  char usegcol;
};

int parse_xmlnames(char *xmlnames, int *nxmls, char ***xmls);

int xml2svg_getpar(struct Parameters* const par, 
		   int *nxmls, 
		   char ***xmls, 
		   Obj2D_instance ***obj);

#endif /* XML2SVG_H */