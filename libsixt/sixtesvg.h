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

#ifndef SIXTESVG_H
#define SIXTESVG_H 1

#include "sixt.h"

/** Data structure that holds a SVG file and the
    neccessary information for the environment handling
    and sizes. */
typedef struct{
  
  /** File pointer to SVG file. */
  FILE* file;
  /** Width of the SVG in pixels. */
  double canvaswidth;
  /** Height of the SVG in pixels. */
  double canvasheight;
  /** Factor to convert from world units to pixels. */
  double scalefactor;
  /** World width of the SVG. */
  double width;
  /** World height of the SVG. */
  double height;
  /** Horizontal offset of world coordinates **/
  double x0;
  /** Vertical offset of world coordinates **/
  double y0;
  /** Init flag. */
  int Init;
  /** Header flag. */
  int Head;
  /** Path flag. */
  int path_open;
  /** Text flag. */
  int text_open;
  
}SixteSVGObj;

/** Function to get a SixteSVGObj and default values. */
SixteSVGObj* getSixteSVGObj(int* const status);

/** Function to initialize a SixteSVGObj and to open file. */
void SixteSVG_init(SixteSVGObj *svg, 
		   char *svgfilename,
		   double worldx0,
		   double worldy0,
		   double worldwidth, 
		   double worldheight, 
		   double svgwidth, 
		   int* const status);

/** Function to end the SVG file and to close an active SixteSVGObj. */
void SixteSVG_close(SixteSVGObj *svg,
		    int* const status);

/** Function to write the header into an initialized SixteSVGObj. */
void SixteSVG_makeHeader(SixteSVGObj *svg,
			 int* const status);

/** Function to draw a basic line between two points. */
void SixteSVG_draw_line(SixteSVGObj *svg, 
			double x1, 
			double y1, 
			double x2, 
			double y2, 
			double linewidth, 
			char* color, 
			int* const status);

/** Function to draw a circle. */
void SixteSVG_draw_circle(SixteSVGObj *svg, 
			  double cx, 
			  double cy, 
			  double R, 
			  double linewidth, 
			  char* linecolor, 
			  char* fillcolor, 
			  int fill, 
			  int* const status);

/** Function to open the path environment and set the starting point. */
void SixteSVG_begin_path(SixteSVGObj *svg, 
			 double sx, 
			 double sy, 
			 int* const status);

/** Function to close the path environment. */
void SixteSVG_close_path(SixteSVGObj *svg, 
			 char* linecolor, 
			 char* fillcolor, 
			 double linewidth, 
			 int fill, 
			 int* const status);

/** Function to add another point to a path in absolute coordinates. */
void SixteSVG_path_line(SixteSVGObj *svg, 
			double x, 
			double y, 
			int* const status);

/** Function to add another point to a path in coordinates relatively 
    to the last point. */
void SixteSVG_path_line_rel(SixteSVGObj *svg, 
			    double dx, 
			    double dy, 
			    int* const status);

/** Function to add another point to a path at a given distance and
    angle. */
void SixteSVG_path_line_angle(SixteSVGObj *svg, 
			      double r, 
			      double ang, 
			      int* const status);

/** Function to draw an arc in a path defined by its start and end
    angle and the radius. */
void SixteSVG_path_arc(SixteSVGObj *svg, 
		       double r, 
		       double ang0, 
		       double ang1, 
		       int* const status);

/** Function to draw a polygon. */
void SixteSVG_draw_polygon(SixteSVGObj *svg,
			   int npoints,
			   double *x, 
			   double *y, 
			   char* linecolor, 
			   char* fillcolor, 
			   double linewidth, 
			   int fill, 
			   int* const status);

/** Function to open the text environment. */
void SixteSVG_open_text(SixteSVGObj *svg, 
			double x, 
			double y, 
			double size, 
			int* const status);

/** Function to close the text environment. */
void SixteSVG_close_text(SixteSVGObj *svg, 
			 int* const status);

/** Function to add text inside the text environment. */
void SixteSVG_write_text(SixteSVGObj *svg, 
			 char* text, 
			 int* const status);

/** Function to write text centered to a given point */
void SixteSVG_write_centered_text(SixteSVGObj *svg, 
			     char* text, 
			     double cx,
			     double cy,
			     double textsize,
			     int* const status);

#endif /* SIXTESVG_H */