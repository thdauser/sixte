#include "projectedmask.h"

ProjectedMask* newProjectedMask(int* const status)
{
  ProjectedMask* proj=(ProjectedMask*)malloc(sizeof(ProjectedMask));
  if (NULL==proj) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for ProjectedMask!\n",
		   *status);
    return(proj);
  }

  //Initialization:
  proj->map=NULL;

  proj->pixelwidth1 = 0.;
  proj->pixelwidth2 = 0.;
  proj->naxis1 = 0;
  proj->naxis2 = 0;

  return(proj);
}


ProjectedMask* getProjectedMask(const CodedMask* const mask, SquarePixels* det_pix,
				const double proj_dist, int* const status)
{
  ProjectedMask* proj=NULL;
  int x,y;                   //counts
  int xcount, ycount;        //counts 


  //Get empty ProjectedMask-object
  proj=newProjectedMask(status);
  if (EXIT_SUCCESS!=*status) return(proj);

  //size of axes in pixels: 2*amount_of_mask_pixels+1
  proj->naxis1=2*mask->naxis1+1;
  proj->naxis2=2*mask->naxis2+1;

  //memory-allocation for the map
  //only amount of pixels important,not their actual alternating sizes
  //memory-allocation for Rmap
  proj->map=(double**)malloc(proj->naxis1*sizeof(double*));
  if(NULL!=proj->map){
    for(x=0; x < proj->naxis1; x++){
      proj->map[x]=(double*)malloc(proj->naxis2*sizeof(double));
	if(NULL==proj->map[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "ProjectedMask!\n", *status);
	  return(proj);
	}
	//Clear the pixels
	for(y=0; y < proj->naxis2; y++){
	  proj->map[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(proj);
  } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "ReconstructionArray!\n", *status);
      return(proj);
  }//end of memory-allocation

  //get the different pixelsizes in meters
  //NOTE: only for square pixels!
  proj->pixelwidth1=proj_dist*(mask->cdelt1-det_pix->xpixelwidth);
  proj->pixelwidth2=proj_dist*det_pix->xpixelwidth;

  //scanning over whole mask and fill in the ProjectedMask-pixels
    //equivalent to former mask pixels
  for(y=0; y<mask->naxis2; y++){
    if(y==0){ycount=1;}else{ycount=ycount+2;}
    for(x=0; x<mask->naxis1; x++){
      if(x==0){xcount=1;}else{xcount=xcount+2;}
      proj->map[xcount][ycount]=(double)mask->map[x][y];
     }
   }
  
  
  //scanning over all ProjectedMask-pixels and fill in the new intermediate ones

  for(y=1; y<(proj->naxis2-1); y+=2){ //all middle odd rows -> contains former mask-pix (each odd col in this row)

	  //1st pix
	  if(proj->map[1][y] == 0){proj->map[0][y]=0.;}
	  else{proj->map[0][y]=0.5;}

	  for(x=2; x<(proj->naxis1-1); x+=2){ //all middle odd cols without former mask pixels (already set)
	    proj->map[x][y]=(proj->map[x-1][y]+proj->map[x+1][y])/2;
	  }

	  //last pix
	  if(proj->map[proj->naxis2-2][y] == 0){proj->map[proj->naxis2-1][y]=0.;}
	  else{proj->map[proj->naxis2-1][y]=0.5;}
  }
 

  for(y=2; y<(proj->naxis2-1); y+=2){ //all middle even rows -> only new intermediate pixels

	  //1st pix
	  proj->map[0][y]=(proj->map[0][y-1]+proj->map[0][y+1])/2;

	  for(x=1; x<(proj->naxis1-1); x+=2){ //all middle odd cols -> derived from former already set mask pixels
	    proj->map[x][y]=(proj->map[x][y-1]+proj->map[x][y+1])/2;
	  }

	  for(x=2; x<(proj->naxis1-1); x+=2){ //all middle even cols -> derived from just set odd cols
	    proj->map[x][y]=(proj->map[x-1][y]+proj->map[x+1][y])/2;
	  }

	  //last pix
	  proj->map[proj->naxis2-1][y]=(proj->map[proj->naxis2-1][y-1]+proj->map[proj->naxis2-1][y+1])/2; 
  }

  
 //special case: 1st row (bottom)
   //1st pix
  proj->map[0][0]=(proj->map[1][0]+proj->map[0][1])/2;

  for(x=1; x<(proj->naxis1-1); x+=2){ //all middle odd cols -> look at former mask-pixels
      if(proj->map[x][1] == 0.){
	proj->map[x][0]=0.;
      }else{
	proj->map[x][0]=0.5;
      }
  }

  for(x=2; x<(proj->naxis1-1); x+=2){ //even cols -> intermediate values
      proj->map[x][0]=(proj->map[x-1][0]+proj->map[x+1][0])/2;
  }

   //last pix
  proj->map[proj->naxis1-1][0]=(proj->map[proj->naxis1-2][0]+proj->map[proj->naxis1-1][1])/2;

  
  //special case: last row (top)
   //1st pix
  proj->map[0][proj->naxis1-1]=(proj->map[0][proj->naxis1-2]+proj->map[1][proj->naxis1-1])/2;

  for(x=1; x<(proj->naxis1-1); x+=2){ //all middle odd cols -> look at former mask-pixels
      if(proj->map[x][proj->naxis1-2] == 0.){
	proj->map[x][proj->naxis1-1]=0.;
      }else{
	proj->map[x][proj->naxis1-1]=0.5;
      }
    }

   for(x=2; x<(proj->naxis1-1); x+=2){ //even cols -> intermediate values
      proj->map[x][proj->naxis1-1]=(proj->map[x-1][proj->naxis1-1]+proj->map[x+1][proj->naxis1-1])/2;
    }
 

  //last pix
  proj->map[proj->naxis1-1][proj->naxis1-1]=(proj->map[proj->naxis1-2][proj->naxis1-1]+
  proj->map[proj->naxis1-1][proj->naxis1-2])/2;
  

  return(proj);
}
