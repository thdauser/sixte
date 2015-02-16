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

#include "sixtesvg.h"

SixteSVGObj* getSixteSVGObj(int* const status){
  
  SixteSVGObj *svg=(SixteSVGObj*)malloc(sizeof(SixteSVGObj));
  if(svg==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for SixteSVGObj failed.");
    return(svg);
  }
  
  svg->canvaswidth=0.;
  svg->canvasheight=0.;
  svg->scalefactor=0.;
  svg->width=0.;
  svg->height=0.;
  svg->x0=0.;
  svg->y0=0.;
  svg->Init=0;
  svg->Head=0;
  svg->path_open=0;
  svg->text_open=0;
  svg->file=NULL;
  
  return svg;
}

void SixteSVG_init(SixteSVGObj *svg, 
		   char *svgfilename,
		   double worldx0,
		   double worldy0,
		   double worldwidth, 
		   double worldheight, 
		   double svgwidth, 
		   int* const status){

  if(svg->Init!=0){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG object already initialized");
    return;
  }
  svg->file=fopen(svgfilename, "w");
  svg->x0=worldx0;
  svg->y0=worldy0;
  svg->width=worldwidth;
  svg->height=worldheight;
  svg->scalefactor=svgwidth/worldwidth;
  svg->canvaswidth=svg->width*svg->scalefactor;
  svg->canvasheight=svg->height*svg->scalefactor;
  svg->Init=1;

}

void SixteSVG_close(SixteSVGObj *svg,
		    int* const status){

  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG can not be closed at this point.");
    return;
  }
  fprintf(svg->file, "\n</svg>");
  fclose(svg->file);
  
  svg->Init=0;
  svg->Head=0;
  svg->path_open=0;

}

void SixteSVG_makeHeader(SixteSVGObj *svg,
			 int* const status){

  if(svg->Init==0||svg->Head==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG Header can not be written to uninitialized file.");
    return;
  }
  fprintf(svg->file, "<svg version=\"1.1\"\n"); 
  fprintf(svg->file, "     baseProfile=\"full\"\n");
  fprintf(svg->file, "     width=\"%.0lfpx\" height=\"%.0lfpx\"\n", svg->canvaswidth, svg->canvasheight); 
  fprintf(svg->file, "     xmlns=\"http://www.w3.org/2000/svg\">\n");
  svg->Head=1;

}


void SixteSVG_draw_line(SixteSVGObj *svg, 
			double x1, 
			double y1, 
			double x2, 
			double y2, 
			double linewidth, 
			char* color, 
			int* const status){

  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG line can not be written to uninitialized file.");
    return;
  }
  double r1[2], r2[2];
  r1[0]=svg->scalefactor*(x1-svg->x0);
  r1[1]=svg->canvasheight-svg->scalefactor*(y1-svg->y0);
  r2[0]=svg->scalefactor*(x2-svg->x0);
  r2[1]=svg->canvasheight-svg->scalefactor*(y2-svg->y0);

  fprintf(svg->file, "<line x1=\"%.3lf\" y1=\"%.3lf\" x2=\"%.3lf\" y2=\"%.3lf\"\n  stroke=\"%s\" stroke-width=\"%.2lfpx\"/>\n", 
	  r1[0], r1[1], r2[0], r2[1], color, linewidth);

}

void SixteSVG_draw_circle(SixteSVGObj *svg, 
			  double cx, 
			  double cy, 
			  double R, 
			  double linewidth, 
			  char* linecolor, 
			  char* fillcolor, 
			  int fill, 
			  int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG circle can not be written to uninitialized file.");
    return;
  }
  double r, x, y;
  r=svg->scalefactor*R;
  x=svg->scalefactor*(cx-svg->x0);
  y=svg->canvasheight-svg->scalefactor*(cy-svg->y0);   
     
  fprintf(svg->file, "<circle cx=\"%.3lf\" cy=\"%.3lf\" r=\"%.3lf\"\n        ", x, y, r);
  if(fill!=0){
    fprintf(svg->file, "style=\"fill:%s ;", fillcolor);
  }else{
    fprintf(svg->file, " fill:none; ");
  }
  fprintf(svg->file, "stroke:%s; stroke-width:%.2lfpx \"/>\n", linecolor, linewidth);
     
}

void SixteSVG_begin_path(SixteSVGObj *svg, 
			 double sx, 
			 double sy, 
			 int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG path can not be created in uninitialized file or if another path is already opened.");
    return;
  }
  double x, y;
  x=svg->scalefactor*(sx-svg->x0);
  y=svg->canvasheight-svg->scalefactor*(sy-svg->y0);
  fprintf(svg->file, "<path d=\"M %.3lf %.3lf ", x, y );
  svg->path_open=1;   
     
}

void SixteSVG_close_path(SixteSVGObj *svg, 
			 char* linecolor, 
			 char* fillcolor, 
			 double linewidth, 
			 int fill, 
			 int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==0||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG path can not be closed in uninitialized file or if another path is already opened.");
    return;
  }
  fprintf(svg->file, "z\"\n      style=\"stroke:%s; stroke-width:%.2lfpx;", linecolor, linewidth);
  if(fill!=0){
    fprintf(svg->file, " fill:%s", fillcolor);
  }else{
    fprintf(svg->file, " fill:none");
  }
  fprintf(svg->file, "\" />\n");
  
  svg->path_open=0;
     
}

void SixteSVG_path_line(SixteSVGObj *svg, 
			double x, 
			double y, 
			int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==0||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG path can not be drawn in uninitialized file or if no path is open.");
    return;
  }
  double rx, ry;
  rx=svg->scalefactor*(x-svg->x0);
  ry=svg->canvasheight-svg->scalefactor*(y-svg->y0);   
  
  fprintf(svg->file, "\n      L %.3lf %.3lf ", rx, ry);
     
}

void SixteSVG_path_line_rel(SixteSVGObj *svg, 
			    double dx, 
			    double dy, 
			    int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==0||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG path can not be drawn in uninitialized file or if no path is open.");
    return;
  }
  double rx, ry;
  rx=svg->scalefactor*dx;
  ry=svg->scalefactor*dy;   
  
  fprintf(svg->file, "\n      l %.3lf %.3lf ", rx, ry);
     
}

void SixteSVG_path_line_angle(SixteSVGObj *svg, 
			      double r, 
			      double ang, 
			      int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==0||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG path can not be drawn in uninitialized file or if no path is open.");
    return;
  }
  double rx, ry;
  rx=svg->scalefactor*r*cos(2*M_PI*ang/360.);
  ry=-svg->scalefactor*r*sin(2*M_PI*ang/360.);   
  
  fprintf(svg->file, "\n      l %.3lf %.3lf ", rx, ry);
     
}

void SixteSVG_path_arc(SixteSVGObj *svg, 
		       double r, 
		       double ang0, 
		       double ang1, 
		       int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==0||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG path can not be drawn in uninitialized file or if no path is open.");
    return;
  }
  double R, s[2], e[2], a0, a1, ha;
  int sf, laf, reverse;
  if(ang1<ang0){
    reverse=1;
    ha=ang1;
    ang1=ang0;
    ang0=ha;
  }else{
    reverse=0;
  }
  R=svg->scalefactor*r;
  a0=ang0*2.*M_PI/360.;
  a1=ang1*2.*M_PI/360.;
  if(a1-a0>M_PI){
    laf=1;
  }else{ 
    laf=0;
  }
  if(reverse==1){
    sf=1;
  }else{
    sf=0;
  }
  
  s[0]=R*cos(a0);
  s[1]=-R*sin(a0);
  e[0]=(R*cos(a1)-s[0]);
  e[1]=(-R*sin(a1)-s[1]);
  if(reverse==1){
    e[0]*=-1.;
    e[1]*=-1.;             
  }
  
  fprintf(svg->file, "\n      a %.3lf,%.3lf 0 %d,%d %.3lf,%.3lf ", R, R, laf, sf, e[0], e[1]);
     
}

void SixteSVG_draw_polygon(SixteSVGObj *svg,
			   int npoints,
			   double *x, 
			   double *y, 
			   char* linecolor, 
			   char* fillcolor, 
			   double linewidth, 
			   int fill, 
			   int* const status){
			     
  
  int ii;
  if(*status!=EXIT_SUCCESS){
    return;
  }
  
  SixteSVG_begin_path(svg, x[0], y[0], status);
  if(*status!=EXIT_SUCCESS){
    SIXT_ERROR("SixteSVG failed to draw polygon.");
    return;
  }
  
  for(ii=1; ii<npoints; ii++){
    SixteSVG_path_line(svg, x[ii], y[ii], status);
    if(*status!=EXIT_SUCCESS){
      SIXT_ERROR("SixteSVG failed to draw polygon.");
      return;
    }
  }
  
  SixteSVG_close_path(svg, 
		      linecolor, 
		      fillcolor, 
		      linewidth, 
		      fill, 
		      status);
  
  if(*status!=EXIT_SUCCESS){
    SIXT_ERROR("SixteSVG failed to draw polygon (close path).");
    return;
  }
  
}

void SixteSVG_open_text(SixteSVGObj *svg, 
			double x, 
			double y, 
			double size, 
			int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==1){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG text can not be written in uninitialized file or if another instance is open.");
    return;
  }  
     
  double r[2];
  
  r[0]=svg->scalefactor*(x-svg->x0);
  r[1]=svg->canvasheight-svg->scalefactor*(y-svg->y0);

  fprintf(svg->file, "\n<text x=\"%.3lf\" y=\"%.3lf\" font-family=\"arial\" font-size=\"%.1lfpx\">\n", r[0], r[1], size*svg->scalefactor);
  svg->text_open=1;
}

void SixteSVG_close_text(SixteSVGObj *svg, 
			 int* const status){
  
  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==0){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG text can not be closed at this point.");
    return;
  } 
     
  fprintf(svg->file, "\n</text>\n");
  svg->text_open=0;
     
}

void SixteSVG_write_text(SixteSVGObj *svg, 
			 char* text, 
			 int* const status){
     
  if(svg->Init==0||svg->Head==0||svg->path_open==1||svg->text_open==0){
    *status=EXIT_FAILURE;
    SIXT_ERROR("SixteSVG text can not be written in uninitialized file or if another instance is open.");
    return;
  }
  
  fprintf(svg->file, "\n%s", text);  
     
}
