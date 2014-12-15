#include "reconstruction.h"

ReconArray* newReconArray(int* const status)
{
  ReconArray* recon=(ReconArray*)malloc(sizeof(ReconArray));
  if (NULL==recon) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Reconstruction Array!\n",
		   *status);
    return(recon);
  }

  //Initialization:
  recon->Rmap=NULL;
  recon->open_fraction=0.;

  recon->naxis1 = 0;
  recon->naxis2 = 0;

  return(recon);
}


ReconArray* getReconArray(CodedMask* mask,int type, SquarePixels* detector_pixels, int* const status)
{
  ReconArray* recon=NULL;
  int x,y;                          //count for memory allocation
  int xcount, ycount;               //count for getting Rmap in case1:same pixel size

  //Get empty reconstruction array-object
  recon=newReconArray(status);
  if (EXIT_SUCCESS!=*status) return(recon);

 //in order to get same pixelsize as detector:
 //[naxis1 in maskpixels]*[pixelsize of mask]->absolut size of mask in meters
 //divided by [detector pixelsize]
    recon->naxis1=2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth);
    recon->naxis2=2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth);

  //memory-allocation for Rmap
  recon->Rmap=(double**)malloc(recon->naxis1*sizeof(double*));
  if(NULL!=recon->Rmap){
    for(x=0; x < recon->naxis1; x++){
      recon->Rmap[x]=(double*)malloc(recon->naxis2*sizeof(double));
	if(NULL==recon->Rmap[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "ReconstructionArray!\n", *status);
	  return(recon);
	}
	//Clear the pixels
	for(y=0; y < recon->naxis2; y++){
	  recon->Rmap[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(recon);
  } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "ReconstructionArray!\n", *status);
      return(recon);
  }//end of memory-allocation

   recon->open_fraction = mask->transparency;

  // decide whether pixelsize is different for det and mask
   if(mask->cdelt1 == detector_pixels->xpixelwidth){

  //Scanning over all mask-elements
  for(xcount=0; xcount < mask->naxis1; xcount++){
    for(ycount=0; ycount < mask->naxis2; ycount++){
      if(mask->map[xcount][ycount]==1){
	recon->Rmap[xcount][ycount]=1.;
      }else if(mask->map[xcount][ycount]==0){

	if(type == 1){
	  recon->Rmap[xcount][ycount]=
	    (recon->open_fraction)/(recon->open_fraction - 1.);
	}else if(type == 2){
	  recon->Rmap[xcount][ycount]=-1;
	}

      }else{
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error while scanning mask-elements!\n", *status);
	}
    }
  }//End of scanning over mask-elements
   }//End equal pixelsize

   else{//begin diff pixelsize
     //initialize all to max neg value (depending on OF)
     double MinVal;
     if(type == 1){
       MinVal=(recon->open_fraction)/(recon->open_fraction - 1.);
     }else if(type == 2){
       MinVal=-1.;
     }

     for(xcount=0; xcount < (recon->naxis1/2); xcount++){
       for(ycount=0; ycount < (recon->naxis2/2); ycount++){
	 recon->Rmap[xcount][ycount]=MinVal;
	}
     }

     repixWithReminder(mask,recon,5,mask->naxis1,mask->naxis2,mask->cdelt1,mask->cdelt1,detector_pixels->xpixelwidth,MinVal);

   }//end diff pixelsize
   return(recon);
}


double* SaveReconArray1d(ReconArray* recon, int* status)
{
 double* ReconArray1d=NULL;

 //Memory-Allocation for 1d-image
 ReconArray1d = (double*)malloc(recon->naxis1*recon->naxis2*sizeof(double));
 if (!ReconArray1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-recon-array!\n", *status);
    return(ReconArray1d);
 }

    //Create the 1D-image from ReconArray 
  int x, y;
  for (x=0; x<recon->naxis1; x++) {
    for (y=0; y<recon->naxis2; y++) {
	ReconArray1d[(x+ recon->naxis1*y)] = recon->Rmap[x][y];
   }
  }
 return(ReconArray1d);  
}


void FreeReconArray1d(double* ReconArray1d)
{
  if (ReconArray1d!=NULL) free(ReconArray1d);
}


void FreeReconArray(ReconArray** const recon)
{
  if ((*recon)!=NULL) {
    if (((*recon)->naxis1>0)&&(NULL!=(*recon)->Rmap)) {
      int count;
      for(count=0; count<(*recon)->naxis1; count++) {
	if (NULL!=(*recon)->Rmap[count]) {
	  free((*recon)->Rmap[count]);
	}
      }
      free((*recon)->Rmap);
    }
    free(*recon);
    *recon=NULL;
  }
}



 
