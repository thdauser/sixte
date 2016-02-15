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

#include "detstruct2obj2d.h"

Obj2D_instance *getObj2DFromAdvdet(AdvDet *det, int* const status){
  
  Obj2D_instance *obj=NULL;
  
  obj=getObj2D_instance(status);
  CHECK_STATUS_RET(*status, NULL);
  
  int n, ii;
  double sx, sy;
  
  do{
  
    n=det->npix;
    sx=det->sx;
    sy=det->sy;
  
    obj->subobj=(Obj2D_instance**)malloc(n*sizeof(Obj2D_instance*));
    CHECK_NULL_BREAK(obj->subobj,*status,"Malloc of Obj2D subobjects failed.")
    
    obj->n_subobj=n;
  
    for(ii=0; ii<n; ii++){
      obj->subobj[ii]=getObj2D_instance(status);
      CHECK_STATUS_BREAK(*status);
      obj->subobj[ii]->geometry=create_rectangular_Obj2D(det->pix[ii].pindex, 
							sx+det->pix[ii].sx,
							sy+det->pix[ii].sy,
							det->pix[ii].width,
							det->pix[ii].height,
							0.,
							status);
      CHECK_STATUS_BREAK(*status);
      if(det->pix[ii].channel!=NULL){
	Obj2D_assign_group_attribute(obj->subobj[ii],
				    &(det->pix[ii].channel->channel_id),
				    &(det->pix[ii].freq));
      }
      obj->subobj[ii]->parent=obj;
    }
    CHECK_STATUS_BREAK(*status);
    
    break;
  }while(0);
  
  if(*status!=EXIT_SUCCESS){
    freeObj2D_instance(obj);
    free(obj);
    obj=NULL;
    SIXT_ERROR("Converting advdet to Obj2D failed.");
  }
  
  return(obj);  
}


Obj2D_instance *getObj2DFromGendet(GenDet *det, int* const status){
  
  Obj2D_instance *obj=NULL;
  
  obj=getObj2D_instance(status);
  CHECK_STATUS_RET(*status, NULL);
  
  int nx, ny, n, ii, jj, ll;
  double dx, dy, w, h, rota, xrpix, yrpix, xrval, 
	 yrval, posx, posy, x, y, sinr, cosr;
  
  nx=det->pixgrid->xwidth;
  ny=det->pixgrid->ywidth;
  
  n=nx*ny;
  
  dx=det->pixgrid->xdelt;
  dy=det->pixgrid->ydelt;
  
  xrpix=det->pixgrid->xrpix;
  yrpix=det->pixgrid->yrpix;
  
  xrval=det->pixgrid->xrval;
  yrval=det->pixgrid->yrval;
  
  w=det->pixgrid->xdelt-2.*det->pixgrid->xborder;
  h=det->pixgrid->ydelt-2.*det->pixgrid->yborder;
  
  rota=det->pixgrid->rota;
	
  sinr=sin(rota);
  cosr=cos(rota);
  
  x=(0.5*(nx+1)-xrpix)*dx;
  y=(0.5*(ny+1)-yrpix)*dy;
  
  posx=xrval+x*cosr-y*sinr;
  posy=yrval+x*sinr+y*cosr;
  
  obj->geometry=create_rectangular_Obj2D(n,
					 posx,
					 posy,
					 nx*dx,
					 ny*dy,
					 rota*180./M_PI,
					 status);
  
  do{   
    obj->subobj=(Obj2D_instance**)malloc(n*sizeof(Obj2D_instance*));
    CHECK_NULL_BREAK(obj->subobj,*status,"Malloc of Obj2D subobjects failed.");
    CHECK_STATUS_BREAK(*status);
    obj->n_subobj=n; 
    ll=0;
    for(jj=0; jj<ny; jj++){
      for(ii=0; ii<nx; ii++){
	x=(1.*(ii+1)-xrpix)*dx;
	y=(1.*(jj+1)-yrpix)*dy;
	
	posx=xrval+x*cosr-y*sinr;
	posy=yrval+x*sinr+y*cosr;
	
	obj->subobj[ll]=NULL;
	obj->subobj[ll]=getObj2D_instance(status);
	CHECK_STATUS_BREAK(*status);
	obj->subobj[ll]->geometry=create_rectangular_Obj2D(ll,
							   posx,
							   posy,
							   w,
							   h,
							   rota*180./M_PI,
							   status);
	CHECK_STATUS_BREAK(*status);	
	obj->subobj[ll]->parent=obj;
	
	ll++;
      }
      CHECK_STATUS_BREAK(*status);
    }
    CHECK_STATUS_BREAK(*status);
    
    break;    
  }while(0);
  
  if(*status!=EXIT_SUCCESS){
    if(obj!=NULL){
      freeObj2D_instance(obj);
      free(obj);
      obj=NULL;
    }
    SIXT_ERROR("Converting gendet to Obj2D failed.");
  }
  return(obj);   
}

Obj2D_instance *getObj2DFromXML(char *XMLName, int* const status){
  
  int detstatus=EXIT_SUCCESS;
  
  Obj2D_instance *obj=NULL;
  printf("Load %s\n", XMLName);
  // First try to interprete the XML as GenDetXML
  GenInst *inst=NULL;
  do{
    inst=loadGenInst(XMLName, 0, &detstatus); 
    if(detstatus==EXIT_SUCCESS && inst->det->pixgrid->ywidth>0){
      obj=getObj2DFromGendet(inst->det, status);
      if(*status!=EXIT_SUCCESS){
	SIXT_ERROR("Converting GenDet XML file failed.");
	break;
      }else{
	puts("GenDet XML file successfully converted to Obj2D.");
      }
    // If it was not an GenDetXML, interprete it as AdvDetXML
    }else{
      detstatus=EXIT_SUCCESS;
      AdvDet *adet=loadAdvDet(XMLName, &detstatus);

      if(detstatus==EXIT_SUCCESS){
	if(adet->channel_file!=NULL){
	  puts("Found crosstalk information in AdvDet XML file.");
	  init_crosstalk(adet, status);
	}
	obj=getObj2DFromAdvdet(adet, status);
	if(*status!=EXIT_SUCCESS){
	  SIXT_ERROR("Converting AdvDet XML file failed.");
	  break;
	}else{
	  puts("AdvDet XML file successfully converted to Obj2D.");
	}
	destroyAdvDet(&adet);
      }else{
	SIXT_ERROR("XML format is not known.");
	*status=detstatus;
      }
    }
  }while(0);
  destroyGenInst(&inst, status);
  if(*status!=EXIT_SUCCESS && obj!=NULL){
    freeObj2D_instance(obj);
    free(obj);
    obj=NULL;
  }
  CHECK_STATUS_RET(*status,obj);
  printf("getObj2DFromXML returns object with %d subobjects.\n", obj->n_subobj);
  return obj;
}

void Obj2D_DrawObjectSVG(Obj2D *obj, 
			 SixteSVGObj *svg,
			 double linewidth,
			 char *linecolor,
			 char *fillcolor,
			 int fill,
			 int writeid,
			 double textsize,
			 int* const status){
				   
  if(obj->type==OBJ2D_CIRC){
    SixteSVG_draw_circle(svg, 
			 obj->cx, 
			 obj->cy, 
			 obj->width, 
			 linewidth, 
			 linecolor, 
			 fillcolor, 
			 fill, 
			 status);
    if(*status!=EXIT_SUCCESS){
      SIXT_ERROR("Drawing Circle Obj2D to SVG failed.");
      return ;
    }
  }else if(obj->nvertices>1){
    SixteSVG_draw_polygon(svg,
			  obj->nvertices,
			  obj->vert_x, 
			  obj->vert_y, 
			  linecolor, 
			  fillcolor, 
			  linewidth, 
			  fill, 
			  status);
    if(*status!=EXIT_SUCCESS){
      SIXT_ERROR("Drawing Polygon Obj2D to SVG failed.");
      return ;
    }
  }

  if(writeid){
    char label[10];
    sprintf(label, "%d", obj->id);
    SixteSVG_write_centered_text(svg, 
				 label, 
				 obj->cx,
				 obj->cy,
				 textsize,
				 status);
    if(*status!=EXIT_SUCCESS){
      SIXT_ERROR("Writing text to SVG failed.");
      return ;
    }
  }
  
}

void Obj2D_DrawInstanceSVG(Obj2D_instance *obj, 
			   SixteSVGObj *svg,
			   char **linecolor,
			   double *linewidth,
			   char **fillcolor,
			   int *fill,
			   int ndraw,
			   int writeid,
			   double *textsize,
			   int usegcol,
			   int* const status){
			     
  int ii, nmax;
  if(obj->geometry!=NULL){
    char *fc=fillcolor[0];
    if(usegcol!=0){
      fc=fillcolor[(obj->geometry->group_id % usegcol)];
    }
    Obj2D_DrawObjectSVG(obj->geometry, 
			svg,
			linewidth[0],
			linecolor[0],
			fc,
			fill[0],
			writeid,
			textsize[0],
			status);
    if(*status!=EXIT_SUCCESS){
      return ;
    }
  }
  
  if(ndraw>=0){
    nmax=ndraw;
    if(nmax>obj->n_subobj){
      nmax=obj->n_subobj;
    }
  }else{
    nmax=obj->n_subobj;
  }
  char **nextfill=&(fillcolor[1]);
  if(usegcol!=0){
    nextfill=fillcolor;
  }
  if(obj->n_subobj>0){
    for(ii=0; ii<nmax; ii++){
      Obj2D_DrawInstanceSVG(obj->subobj[ii], 
			    svg,
			    &(linecolor[1]),
			    &(linewidth[1]),
			    nextfill,
			    &(fill[1]),
			    ndraw,
			    writeid,
			    &(textsize[1]),
			    usegcol,
			    status);
      if(*status!=EXIT_SUCCESS){
	return ;
      }
    }
  }
  
}
