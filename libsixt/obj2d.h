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
   Copyright 2016-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef OBJ2D_H
#define OBJ2D_H 1

#include "sixt.h"

typedef enum {
  OBJ2D_UNDEF	= 0,
  OBJ2D_RECT	= 1,
  OBJ2D_CIRC	= 2
} Obj2D_Type;

/** Type definitions */

/** Data structure that gives information about the geometrical
properties of an object in two dimensions as seen from the illuminated
side of the detector. All coordinates are absolute [m].*/
typedef struct{

  /** Type identifier. */
  Obj2D_Type type;

  /** Object ID. */
  int id;

  /** Width for rectangular object, radius for circle [m]. */
  double width;
  /** Height for rectangular object [m]. */
  double height;

  /** Minimum of bounding box in x-direction [m]. */
  double bbxmin;
  /** Maximum of bounding box in x-direction [m]. */
  double bbxmax;
  /** Minimum of bounding box in y-direction [m]. */
  double bbymin;
  /** Maximum of bounding box in y-direction [m]. */
  double bbymax;

  /** Absolute x-coordinate of central point [m]. */
  double cx;
  /** Absolute y-coordinate of central point [m]. */
  double cy;

  /** Rotation angle [deg] around central point. Counterclockwise
   seen from the illuminated side. */
  double rota;

  /** Number of vertices for rectangulas and other polygons. */
  int nvertices;

  /** x-coordinate array for rectangulas and other polygons [m]. */
  double *vert_x;
  /** y-coordinate array for rectangulas and other polygons [m]. */
  double *vert_y;

  /** group to which the pixel belongs */
  int group_id;
  /** attribute of the pixel */
  double attribute;

}Obj2D;

/** Data structure defining one geometrical object and the
reference to the sub-objects. This structure can be used to
define a detector, the sub-objects can be the individual pixels.
All cordinates contained within this structure and the sub-
structures must be absolute [m]. */
typedef struct Obj2D_instance{

  /** Geometry description of the instance. */
  Obj2D *geometry;

  /** Number of sub-objects of this instance. */
  int n_subobj;
  /** References to the sub-objects. */
  struct Obj2D_instance **subobj;

  /** References to the parent. */
  struct Obj2D_instance *parent;

}Obj2D_instance;

/** Function declarations */

/** Initialization function for Obj2D structure. */
Obj2D *getObj2D(int* const status);

/** Free function for Obj2D structure. */
void freeObj2D(Obj2D *obj);

/** Initialization function for Obj2D_instance structure. */
Obj2D_instance *getObj2D_instance(int* const status);

/** Free function for Obj2D_instance structure. */
void freeObj2D_instance(Obj2D_instance *obj);

/** Function to create a rectangular object as Obj2D. */
Obj2D *create_rectangular_Obj2D(int id,
				double cx,
				double cy,
				double w,
				double h,
				double rota,
				int* const status);

/** Function to create a circular object as Obj2D. */
Obj2D *create_circular_Obj2D(int id,
			     double cx,
			     double cy,
			     double r,
			     int* const status);

/** Function to find the bounding box limits of all geometries
    contained in the object and it's subobjects. Calls itself
    recursively for all sublevels. Needs pre-initialized values
    for x/y max/min. */
void Obj2D_inst_findBBLimits_recursively(Obj2D_instance *obj,
					 double *xmin,
					 double *xmax,
					 double *ymin,
					 double *ymax);

/** Function to find the absolute bounding box of all geometries
    contained in the object and all subobjects. Default values
    of x/y max/min are neglected. Uses
    Obj2D_inst_findBBLimits_recursively after finding start values.*/
void Obj2D_inst_findBBLimits(Obj2D_instance *obj,
			     double *xmin,
			     double *xmax,
			     double *ymin,
			     double *ymax);

/** Function to assign a group id and an attribute to an object. */
void Obj2D_assign_group_attribute(Obj2D_instance *obj,
				  int *group_id,
				  double *attribute);



#endif /* OBJ2D_H */
