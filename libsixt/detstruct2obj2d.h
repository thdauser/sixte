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

#ifndef DET2OBJ2D_H
#define DET2OBJ2D_H 1

#include "sixt.h"
#include "obj2d.h"
#include "advdet.h"
#include "gendet.h"
#include "geninst.h"
#include "sixtesvg.h"

/** Function declarations */

/** Copy the information from an advdet-object into a new Obj2D object */
Obj2D_instance *getObj2DFromAdvdet(AdvDet *det, int* const status);

/** Copy the information from a gendet-object into a new Obj2D object */
Obj2D_instance *getObj2DFromGendet(GenDet *det, int* const status);

/** Function to get an Obj2D from any SIXTE XML */
Obj2D_instance *getObj2DFromXML(char *XMLName, int* const status);

/** Draw an Obj2D to a SVG file. */
void Obj2D_DrawObjectSVG(Obj2D *obj, 
			 SixteSVGObj *svg,
			 double linewidth,
			 char *linecolor,
			 char *fillcolor,
			 int fill,
			 int writeid,
			 double textsize,
			 int* const status);

/** Draw an Obj2D_instance recursively to a SVG file. */
void Obj2D_DrawInstanceSVG(Obj2D_instance *obj, 
			   SixteSVGObj *svg,
			   char **linecolor,
			   double *linewidth,
			   char **fillcolor,
			   int *fill,
			   int ndraw,
			   int writeid,
			   double *textsize,
			   int* const status);


#endif /* DET2OBJ2D_H */