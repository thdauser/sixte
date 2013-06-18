#include "reconstruction.h"

ReconArrayFits* newReconArrayForFits(int* const status)
{
  ReconArrayFits* recon=(ReconArrayFits*)malloc(sizeof(ReconArrayFits));
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
  recon->cdelt1 = 0.;
  recon->cdelt2 = 0.;
  recon->crpix1 = 0.;
  recon->crpix2 = 0.;
  recon->crval1 = 0.;
  recon->crval2 = 0.;

  return(recon);
}

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

ReconArrayFits* getReconArrayForFits(const CodedMask* const mask, int* const status)
{

  ReconArrayFits* recon=NULL;
  int x,y;
  int xcount, ycount;

  //Get empty reconstruction array-object
  recon=newReconArrayForFits(status);
  if (EXIT_SUCCESS!=*status) return(recon);

  //properties are equal to those of the mask
  recon->naxis1=mask->naxis1;
  recon->naxis2=mask->naxis2;
  recon->cdelt1=mask->cdelt1;
  recon->cdelt2=mask->cdelt2;
  recon->crpix1=mask->crpix1;
  recon->crpix2=mask->crpix2;
  recon->crval1=mask->crval1;
  recon->crval2=mask->crval2;

  //memory-allocation
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

  //Scanning over all mask-elements
  for(xcount=0; xcount < mask->naxis1; xcount++){
    for(ycount=0; ycount < mask->naxis2; ycount++){
      if(mask->map[xcount][ycount]==1){
	recon->Rmap[xcount][ycount]=1.;
      }else if(mask->map[xcount][ycount]==0){
	recon->Rmap[xcount][ycount]=(recon->open_fraction)/(recon->open_fraction - 1.);
      }else{
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error while scanning mask-elements!\n", *status);
      }
      
    }
  }//End of scanning over mask-elements

  return(recon);
}


ReconArray* getReconArray(const CodedMask* const mask, SquarePixels* detector_pixels, int* const status)
{

  ReconArray* recon=NULL;
  int x,y;                          //count for memory allocation
  int xcount, ycount;               //count for getting Rmap in case1:same pixel size
  int xpixelcount=0, ypixelcount=0; //count for getting Rmap in case2:diff pixel size
  double left=0.,leftbig=0.;        //left border of small(detector) and big(mask) pixel
  double top=0.,topbig=0.;          //top border of small(detector) and big(mask) pixel
  double ry, rx;
  double w, h;                      //width and height of part of small pixel inside current big pixel
  

  //Get empty reconstruction array-object
  recon=newReconArray(status);
  if (EXIT_SUCCESS!=*status) return(recon);

 //in order to get same pixelsize as detector:
 //[naxis1 in maskpixels]*[pixelsize of mask]->absolut size of mask in meters
 //divided by [detector pixelsize]
  
  /*recon->naxis1=mask->naxis1;
    recon->naxis2=mask->naxis2;*/

  
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


  //for testing: Loop to decide whether pixelsize is different for det and mask
   if(mask->cdelt1 == detector_pixels->xpixelwidth){



  //Scanning over all mask-elements
  for(xcount=0; xcount < mask->naxis1; xcount++){
    for(ycount=0; ycount < mask->naxis2; ycount++){
      if(mask->map[xcount][ycount]==1){


		recon->Rmap[xcount][ycount]=1.;
      }else if(mask->map[xcount][ycount]==0){
	//Recon Array1
	recon->Rmap[xcount][ycount]=
	(recon->open_fraction)/(recon->open_fraction - 1.);

	//padding:move Rmap to middle 
	/*	recon->Rmap[xcount+(recon->naxis1)/4-1][ycount+(recon->naxis2)/4-1]=1.;
      }else if(mask->map[xcount][ycount]==0){
	//Recon Array1
	recon->Rmap[xcount+(recon->naxis1)/4-1][ycount+(recon->naxis2)/4-1]=
	(recon->open_fraction)/(recon->open_fraction - 1.);*/



	//padding:move Rmap to upper right corner 
	/*	recon->Rmap[xcount+(recon->naxis1)/2-1][ycount+(recon->naxis2)/2-1]=1.;
      }else if(mask->map[xcount][ycount]==0){
	//Recon Array1
	recon->Rmap[xcount+(recon->naxis1)/2-1][ycount+(recon->naxis2)/2-1]=
	(recon->open_fraction)/(recon->open_fraction - 1.);*/
      }else{
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error while scanning mask-elements!\n", *status);
	}
      
    }
}//End of scanning over mask-elements
   }//End equal pixelsize

   else{//begin diff pixelsize
  //Scanning over all mask-elements to get Rmap with same pixel-size as detector
     for(ycount=0; ycount<mask->naxis2;ycount++){
       for(xcount=0; xcount<mask->naxis1;xcount++){

	 topbig=ycount*mask->cdelt2;   //top of current big pixel
	 top=topbig;                   //top of current small pixel starts at big pixel
	 ypixelcount=topbig/detector_pixels->ypixelwidth;  //count for small pixel(Rmap)
	 while(top< (topbig+mask->cdelt2)){
	   ry=0.;
	   if((top+detector_pixels->ypixelwidth)>(topbig+mask->cdelt2)){
	     ry=top+detector_pixels->ypixelwidth-topbig+mask->cdelt2;
	   }//end (last) small pixel partially inside big pixel
	   
	   h=detector_pixels->ypixelwidth - (top-ypixelcount*detector_pixels->ypixelwidth)-ry;
	   
	   leftbig=xcount*mask->cdelt1;
	   left=leftbig;
	   xpixelcount=leftbig/detector_pixels->xpixelwidth;
	   while(left< (leftbig+mask->cdelt1)){
	     rx=0.;
	     if((left+detector_pixels->xpixelwidth)>(leftbig+mask->cdelt1)){
	     rx=left+detector_pixels->xpixelwidth-leftbig+mask->cdelt1;
	     }//end (last) small pixel partially inside big pixel
	     
	     w=detector_pixels->xpixelwidth - (left-xpixelcount*detector_pixels->xpixelwidth)-rx;

	     if(mask->map[xcount][ycount]==1){
	       recon->Rmap[xpixelcount][ypixelcount]+=
		h*w/(detector_pixels->xpixelwidth*detector_pixels->ypixelwidth);
	     }//end transparent pixel
	     else{
	       recon->Rmap[xpixelcount][ypixelcount]+=
		 h*w/(detector_pixels->xpixelwidth*detector_pixels->ypixelwidth)*
		 (recon->open_fraction)/(recon->open_fraction - 1.);
	     }//end opaque pixel

	     xpixelcount++;
	     left+=left+w;
	   }//end sanning x small pixel

	   ypixelcount++;
	   top+=top+h;
	 }//end scanning y small pixel

       }//end x bigpixel
     }//end y bigpixel
   }
  return(recon);
}


void SaveReconArrayFits(ReconArrayFits* recon, char* filename, int* status)
{
  fitsfile *fptr=NULL;
  double *image1d=NULL;

  // Print information to STDOUT.
  char msg[MAXMSG];
  sprintf(msg, "Store ReconArray in file '%s' ...\n", filename);
  headas_chat(5, msg);
 
  do { // ERROR handling loop


    // If the specified file already exists, remove the old version.
    remove(filename);
 
    // Create a new FITS-file:
    if (fits_create_file(&fptr, filename, status)) break;

    // Allocate memory for the 1-dimensional image buffer (required for
    // output to FITS file).
    image1d = (double*)malloc(recon->naxis1*recon->naxis2*sizeof(double));
    if (!image1d) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error allocating memory!\n", *status);
      break;
    }

    // Store the ReconArray in the 1-dimensional buffer to handle it 
    // to the FITS routine.
    int x, y;
    for (x=0; x<recon->naxis1; x++) {
      for (y=0; y<recon->naxis2; y++) {
	image1d[(x+ recon->naxis1*y)] = recon->Rmap[x][y];
      }
    }

    
    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = {(long)(recon->naxis1), (long)(recon->naxis2)};
    if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status)) break;
    //                                   |-> naxis
    //    int hdutype;
    if (fits_movabs_hdu(fptr, 1, NULL, status)) break;

    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1)
    // Upper right corner.
    long lpixel[2] = {recon->naxis1, recon->naxis2}; 
    fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, image1d, status);

  } while (0); // END of ERROR handling loop

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  if (NULL!=image1d) free(image1d);
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

 //for testing: multiply mask and ReconArray to see whether it equals the unity matrix
double* MultiplyMaskRecon(ReconArray* recon, CodedMask* mask, int* status)
{
  double* multiplyMR1d=NULL;
  double** Multiply=NULL;
  double sum=0.;
  int i,j,k;
  
//Memory-Allocation for 1d-image
  multiplyMR1d = (double*)malloc((mask->naxis1)*(mask->naxis2)*sizeof(double));
 if (!multiplyMR1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-recon-array!\n", *status);
    return(multiplyMR1d);
 }


//memory-allocation for Multiply
 Multiply=(double**)malloc(mask->naxis1*sizeof(double*));
  if(NULL!=Multiply){
    for(i=0; i < mask->naxis1; i++){
      Multiply[i]=(double*)malloc(mask->naxis2*sizeof(double));
	if(NULL==Multiply[i]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "MultiplyArray!\n", *status);
	  return(Multiply);
	}
	//Clear the pixels
	for(j=0; j < mask->naxis2; j++){
	 Multiply[i][j]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(Multiply);
  } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "MultiplyArray!\n", *status);
      return(Multiply);
  }//end of memory-allocation

 for(i=0; i<mask->naxis1; i++)
   {
     for(j=0; j<mask->naxis2; j++)
    {
      for(k=0; k<mask->naxis2; k++)
	{
	  sum+=mask->map[i][k] * recon->Rmap[k + mask->naxis1][j + mask->naxis2];
	  Multiply[i][j]=sum;
	}
      sum=0.;
    }
   }


  //Create the 1D-image from Multiply 
  for (i=0; i<mask->naxis1; i++) {
    for (j=0; j<mask->naxis2; j++) {
	multiplyMR1d[(i+ mask->naxis1*j)] = Multiply[i][j];
   }
  }


  free(Multiply);
  return(multiplyMR1d);
}


void FreeReconArray1d(double* ReconArray1d)
{
  if (ReconArray1d!=NULL) free(ReconArray1d);
}

void FreeReconArrayFits(ReconArrayFits* recon)
{
  if (recon!=NULL) {
    if ((recon->naxis1>0)&&(NULL!=recon->Rmap)) {
      int count;
      for(count=0; count< recon->naxis1; count++) {
	if (NULL!=recon->Rmap[count]) {
	  free(recon->Rmap[count]);
	}
      }
      free(recon->Rmap);
    }
    free(recon);
  }
}

void FreeReconArray(ReconArray* recon)
{
  if (recon!=NULL) {
    if ((recon->naxis1>0)&&(NULL!=recon->Rmap)) {
      int count;
      for(count=0; count< recon->naxis1; count++) {
	if (NULL!=recon->Rmap[count]) {
	  free(recon->Rmap[count]);
	}
      }
      free(recon->Rmap);
    }
    free(recon);
  }
}
