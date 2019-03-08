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

#include "obj2d.h"

Obj2D *getObj2D(int* const status){

  Obj2D *obj=NULL;

  obj=(Obj2D*)malloc(sizeof(Obj2D));
  if(obj==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for Obj2D failed.");
    return(obj);
  }

  obj->id=0;
  obj->width=0.;
  obj->height=0.;
  obj->bbxmin=0.;
  obj->bbxmax=0.;
  obj->bbymin=0.;
  obj->bbymax=0.;
  obj->cx=0.;
  obj->cy=0.;
  obj->rota=0.;
  obj->nvertices=0;
  obj->vert_x=NULL;
  obj->vert_y=NULL;
  obj->type=OBJ2D_UNDEF;
  obj->group_id=0;
  obj->attribute=0.;

  return obj;

}

void freeObj2D(Obj2D *obj){

  if(obj->vert_x!=NULL){
    free(obj->vert_x);
    obj->vert_x=NULL;
  }
  if(obj->vert_y!=NULL){
    free(obj->vert_y);
    obj->vert_y=NULL;
  }

  obj->id=0;
  obj->width=0.;
  obj->height=0.;
  obj->bbxmin=0.;
  obj->bbxmax=0.;
  obj->bbymin=0.;
  obj->bbymax=0.;
  obj->cx=0.;
  obj->cy=0.;
  obj->rota=0.;
  obj->nvertices=0;
  obj->type=OBJ2D_UNDEF;
  obj->group_id=0;
  obj->attribute=0.;
}

Obj2D_instance *getObj2D_instance(int* const status){

  Obj2D_instance *obj=NULL;

  obj=(Obj2D_instance*)malloc(sizeof(Obj2D_instance));
  if(obj==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for Obj2D_instance failed.");
    return(obj);
  }

  obj->geometry=NULL;

  obj->n_subobj=0;
  obj->subobj=NULL;

  obj->parent=NULL;

  return obj;
}

void freeObj2D_instance(Obj2D_instance *obj){

  int ii;

  if(obj!=NULL){
    if(obj->subobj!=NULL && obj->n_subobj>0){
      for(ii=0; ii<(obj->n_subobj); ii++){
    	  if(obj->subobj[ii]!=NULL){
    		  freeObj2D_instance(obj->subobj[ii]);
    		  free(obj->subobj[ii]);
    	  }
      }
      free(obj->subobj);
      obj->subobj=NULL;
    }

    obj->n_subobj=0;
    if(obj->geometry!=NULL){
      freeObj2D(obj->geometry);
      free(obj->geometry);
      obj->geometry=NULL;
    }
  }
}

Obj2D *create_rectangular_Obj2D(int id,
				double cx,
				double cy,
				double w,
				double h,
				double rota,
				int* const status){

  int ii;

  Obj2D *obj=getObj2D(status);
  CHECK_STATUS_RET(*status, NULL);

  obj->id=id;
  obj->type=OBJ2D_RECT;
  obj->width=w;
  obj->height=h;
  obj->cx=cx;
  obj->cy=cy;
  obj->rota=rota;

  obj->vert_x=(double*)malloc(4*sizeof(double));
  if(obj->vert_x==NULL){
    freeObj2D(obj);
    free(obj);
    obj=NULL;
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for rectangular Obj2D vertex array failed.");
    return(obj);
  }

  obj->vert_y=(double*)malloc(4*sizeof(double));
  if(obj->vert_y==NULL){
    freeObj2D(obj);
    free(obj);
    obj=NULL;
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for rectangular Obj2D vertex array failed.");
    return(obj);
  }

  obj->nvertices=4;

  // Set coordinates according to pixel width and height
  double px[4], py[4];
  px[0]=0.5*w;
  px[1]=-0.5*w;
  px[2]=-0.5*w;
  px[3]=0.5*w;

  py[0]=0.5*h;
  py[1]=0.5*h;
  py[2]=-0.5*h;
  py[3]=-0.5*h;

  // Rotate them around the center and add center position
  double sinr, cosr;
  sinr=sin(rota*M_PI/180.);
  cosr=cos(rota*M_PI/180.);
  for(ii=0; ii<4; ii++){
    obj->vert_x[ii] = cx + px[ii]*cosr - py[ii]*sinr;
    obj->vert_y[ii] = cy + px[ii]*sinr + py[ii]*cosr;
  }

  // Find bounding box borders
  obj->bbxmin=obj->vert_x[0];
  obj->bbxmax=obj->vert_x[0];
  obj->bbymin=obj->vert_y[0];
  obj->bbymax=obj->vert_y[0];
  for(ii=1; ii<4; ii++){
    obj->bbxmin=MIN(obj->bbxmin, obj->vert_x[ii]);
    obj->bbxmax=MAX(obj->bbxmax, obj->vert_x[ii]);
    obj->bbymin=MIN(obj->bbymin, obj->vert_y[ii]);
    obj->bbymax=MAX(obj->bbymax, obj->vert_y[ii]);
  }

  return obj;
}

Obj2D *create_circular_Obj2D(int id,
			     double cx,
			     double cy,
			     double r,
			     int* const status){

  Obj2D *obj=getObj2D(status);
  CHECK_STATUS_RET(*status, NULL);

  obj->id=id;
  obj->type=OBJ2D_CIRC;
  obj->width=r;
  obj->height=r;
  obj->cx=cx;
  obj->cy=cy;

  obj->bbxmin=cx-r;
  obj->bbxmax=cx+r;
  obj->bbymin=cy-r;
  obj->bbymax=cy+r;

  return obj;
}

void Obj2D_inst_findBBLimits_recursively(Obj2D_instance *obj,
					 double *xmin,
					 double *xmax,
					 double *ymin,
					 double *ymax){

  int ii;
  int n=obj->n_subobj;

  if(obj->geometry!=NULL){
    *xmin=MIN(*xmin, obj->geometry->bbxmin);
    *xmax=MAX(*xmax, obj->geometry->bbxmax);
    *ymin=MIN(*ymin, obj->geometry->bbymin);
    *ymax=MAX(*ymax, obj->geometry->bbymax);
  }

  if(obj->subobj!=NULL){
    for(ii=0; ii<n; ii++){
      Obj2D_inst_findBBLimits_recursively(obj->subobj[ii],
					xmin, xmax, ymin, ymax);
    }
  }

}

void Obj2D_inst_findBBLimits(Obj2D_instance *obj,
			     double *xmin,
			     double *xmax,
			     double *ymin,
			     double *ymax){


  int n=obj->n_subobj;
  int ii;

  if(obj->geometry!=NULL){
    *xmin=obj->geometry->bbxmin;
    *xmax=obj->geometry->bbxmax;
    *ymin=obj->geometry->bbymin;
    *ymax=obj->geometry->bbymax;
  }else if(n>0){
    *xmin=obj->subobj[0]->geometry->bbxmin;
    *xmax=obj->subobj[0]->geometry->bbxmax;
    *ymin=obj->subobj[0]->geometry->bbymin;
    *ymax=obj->subobj[0]->geometry->bbymax;
  }else{
    return ;
  }

  if(obj->subobj!=NULL){
    for(ii=0; ii<n; ii++){
      Obj2D_inst_findBBLimits_recursively(obj->subobj[ii],
					xmin, xmax, ymin, ymax);
    }
  }

}

void Obj2D_assign_group_attribute(Obj2D_instance *obj,
				  int *group_id,
				  double *attribute){

  if(obj!=NULL){
    if(group_id!=NULL){
      obj->geometry->group_id=(*group_id);
    }
    if(attribute!=NULL){
      obj->geometry->attribute=(*attribute);
    }
  }
}
