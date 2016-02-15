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

#include "xml2svg.h"

#define OBJ2D_NCOL 30
/** static color list */
static char *obj2d_fc[]={"#ff7979",
  			 "#72bfd2",
  			 "#d2bf72",
  			 "#d272c6",
  			 "#78d272",
  			 "#b1d272",
  			 "#72d2ad",
  			 "#edb2bb",
  			 "#ede7b2",
  			 "#edceb2",
  			 "#b2edcb",
  			 "#b2e3ed",
  			 "#b9b2ed",
  			 "#e5b2ed",
  			 "#d9d9d9",
  			 "#828282",
  			 "#ff34c3",
  			 "#34ff47",
  			 "#fff834",
  			 "#ff4234",
  			 "#5a34ff",
  			 "#34b5ff",
  			 "#0a9f3f",
  			 "#c13d1a",
  			 "#c17c1a",
  			 "#9cc11a",
  			 "#1ac17c",
  			 "#291ac1",
  			 "#861bae",
  			 "#875681",
			 "#d2143c"
};

/** Main */

int xml2svg_main() {
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  
  // Error status.
  int status=EXIT_SUCCESS;
  
  int nxmls=0;
  char **xmls=NULL;
  
  double xmin, xmax, ymin, ymax, xxmin, xxmax, yymin, yymax;
  double worldxmin, worldxmax, worldymin, worldymax, worldborder;
  
  char *linecolor[]={"black", "red"};
  char *fillcolor[]={"#ffeeaa", "white"};
  double linewidth[2]={2.0, 1.0};
  double textsize[2]={5.,4.};
  int fill[2]={1, 0};
  
  Obj2D_instance **obj=NULL;
  SixteSVGObj *svg=NULL;
  
  int ii;
  
  // Register HEATOOL:
  set_toolname("xml2svg");
  set_toolversion("0.05");
  
  do{
    
    // ---- Initialization ----   
    // Read the parameters using PIL.
    status=xml2svg_getpar(&par, &nxmls, &xmls, &obj);
    CHECK_STATUS_BREAK(status);    
    for(ii=0; ii<nxmls; ii++){
      obj[ii]=getObj2DFromXML(xmls[ii], &status);
      CHECK_STATUS_BREAK(status);
    }    
    CHECK_STATUS_BREAK(status);
    int writeid=0;
    if(par.writeid!=0){
      writeid=1;
    }
    int usegcol=0;
    char **fillc=fillcolor;
    if(par.usegcol!=0){
      usegcol=OBJ2D_NCOL;
      fillc=obj2d_fc;
      fill[1]=1;
    }
    // Find min, max values in all objects.
    Obj2D_inst_findBBLimits(obj[0], &xmin, &xmax, &ymin, &ymax);
    if(nxmls>1){
      for(ii=1; ii<nxmls; ii++){
	Obj2D_inst_findBBLimits(obj[ii], &xxmin, &xxmax, &yymin, &yymax);
	xmin=MIN(xmin, xxmin);
	xmax=MAX(xmax, xxmax);
	ymin=MIN(ymin, yymin);
	ymax=MAX(ymax, yymax);
      }
    }    
    // Open the SVG objects
    worldborder=par.border*(xmax-xmin)/par.svgwidth;
    worldxmin=xmin-worldborder;
    worldxmax=xmax+worldborder;
    worldymin=ymin-worldborder;
    worldymax=ymax+worldborder;
    
    svg=getSixteSVGObj(&status);
    SixteSVG_init(svg, 
		  par.SVGName, 
		  worldxmin, 
		  worldymin, 
		  (worldxmax-worldxmin),
		  (worldymax-worldymin),
		  par.svgwidth,
		  &status);
    CHECK_STATUS_BREAK(status);
    SixteSVG_makeHeader(svg, &status);   
    // Draw all Objects    
    for(ii=0; ii<nxmls; ii++){
      Obj2D_DrawInstanceSVG(obj[ii], 
			    svg, 
			    linecolor, 
			    linewidth, 
			    fillc, 
			    fill,
			    par.drawn,
			    writeid,
			    textsize,
			    usegcol,
			    &status);
    CHECK_STATUS_BREAK(status);
    }
    // Draw coordinate system axes
    SixteSVG_draw_line(svg, 
		       0., worldymin, 
		       0., worldymax, 
		       1., "#999999", 
		       &status);
    SixteSVG_draw_line(svg, 
		       worldxmin, 0., 
		       worldxmax, 0.,
		       1., "#999999", 
		       &status);
    CHECK_STATUS_BREAK(status);
    // Close SVG Objects
    SixteSVG_close(svg, &status);
    CHECK_STATUS_BREAK(status);
    
  }while(0); // END of the error handling loop.
  
  if(obj!=NULL){
    for(ii=0; ii<nxmls; ii++){
      freeObj2D_instance(obj[ii]);
      obj[ii]=0;
    }
    free(obj);
    obj=NULL;
  }  
  if(xmls!=NULL){
    for(ii=0; ii<nxmls; ii++){
      free(xmls[ii]);
      xmls[ii]=0;
    }
    free(xmls);
    xmls=NULL;
  }
  
  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}

int parse_xmlnames(char *xmlnames, int *nxmls, char ***xmls){
  
  char *ptr=xmlnames;
  char buffer[MAXFILENAME];
  int ll;
  
  *nxmls=0;
  
  while(1){
    ll=0;
    while(*ptr!=' ' && *ptr!=';' && *ptr!=',' && *ptr!='\'' && *ptr!='\"' && ll<MAXFILENAME-2 && *ptr!='\0'){
      buffer[ll]=*ptr;
      ptr++;
      ll++;
    }
    buffer[ll]='\0';
    *xmls=(char**)realloc(*xmls, (*nxmls+1)*sizeof(char*));
    if(*xmls==NULL){
      SIXT_ERROR("failed allocating memory for XML file names.");
      return EXIT_FAILURE;
    }
    (*xmls)[*nxmls]=(char*)malloc(((strlen(buffer)+1)*sizeof(char)));
    if((*xmls)[*nxmls]==NULL){
      SIXT_ERROR("failed allocating memory for XML file names.");
      return EXIT_FAILURE;
    }
    strcpy((*xmls)[*nxmls], buffer);
    (*nxmls)++;
    if(*ptr==' ' || *ptr==',' || *ptr==';'){
      ptr=ptr+1;
    }else{
      break;
    }
  }
  
  return EXIT_SUCCESS;
}

int xml2svg_getpar(struct Parameters* const par, 
		   int *nxmls, 
		   char ***xmls, 
		   Obj2D_instance ***obj){

  // Error status.
  int status=EXIT_SUCCESS;
  
  int ii;
  
  status=ape_trad_query_string("XMLFiles", &(par->XMLFiles));
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the XML file(s).");
    return(status);
  }
  
  status=ape_trad_query_string("SVGName", &(par->SVGName));
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the output SVG file name.");
    return(status);
  }
  
  status=ape_trad_query_double("SVGWidth", &par->svgwidth);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading SVGWidth");
    return(status);
  }
  
  status=ape_trad_query_double("Border", &par->border);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading Border");
    return(status);
  }
  
  status=ape_trad_query_int("DrawN", &par->drawn);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading DrawN");
    return(status);
  }
  
  status=parse_xmlnames(par->XMLFiles, nxmls, xmls);
  
  if(status==EXIT_SUCCESS){
    printf("%d XML files were found:\n", *nxmls);
    for(ii=0; ii<*nxmls; ii++){
      printf("#%d: %s\n", ii+1, (*xmls)[ii]);
    }
  }else{
    SIXT_ERROR("Failed to parse XML files.");
    return status;
  }
  
  (*obj)=(Obj2D_instance**)malloc(*nxmls*sizeof(Obj2D_instance*));
  if(*obj==NULL){
    status=EXIT_FAILURE;
    SIXT_ERROR("Failed to allocate memory for Obj2D instances.");
    return status;
  }
  for(ii=0; ii<(*nxmls); ii++){
    (*obj)[ii]=NULL;
  }
  
  status=ape_trad_query_bool("WriteID", &par->writeid);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the WriteID parameter");
    return(status);
  }
  
  status=ape_trad_query_bool("UseGCol", &par->usegcol);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the UseGCol parameter");
    return(status);
  }
  
  return status;
}
