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

  return(recon);
}

ReconArrayFits* getReconArrayForFits(const CodedMask* const mask, int* const status)
{

  ReconArrayFits* recon=NULL;
  int x,y;
  int xcount, ycount;

  //Get empty reconstruction array-object
  recon=newReconArray(status);
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
	//Recon Array1
	recon->Rmap[xcount][ycount]=(recon->open_fraction)/(recon->open_fraction - 1.);
	//ReconArray2
	//double total_elements = (double)(mask->naxis1*mask->naxis2);

	//recon->Rmap[xcount][ycount]=
	// ((total_elements*recon->open_fraction-1.)/(total_elements*recon->open_fraction))*
	// ((recon->open_fraction)/(recon->open_fraction - 1.));

      }else{
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error while scanning mask-elements!\n", *status);
      }
      
    }
  }//End of scanning over mask-elements

  return(recon);
}


ReconArray* getReconArray(const CodedMask* const mask, int* const status)
{

  ReconArray* recon=NULL;
  int x,y;
  int xcount, ycount;

  //Get empty reconstruction array-object
  recon=newReconArray(status);
  if (EXIT_SUCCESS!=*status) return(recon);

  //memory-allocation
  recon->Rmap=(double**)malloc(mask->naxis1*sizeof(double*));
  if(NULL!=recon->Rmap){
    for(x=0; x < recon->naxis1; x++){
      recon->Rmap[x]=(double*)malloc(mask->naxis2*sizeof(double));
	if(NULL==recon->Rmap[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "ReconstructionArray!\n", *status);
	 
	  return(recon);
	}
	//Clear the pixels
	for(y=0; y < mask->naxis2; y++){
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
	//Recon Array1
	recon->Rmap[xcount][ycount]=(recon->open_fraction)/(recon->open_fraction - 1.);
      }else{
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error while scanning mask-elements!\n", *status);
      }
      
    }
  }//End of scanning over mask-elements

  return(recon);
}


void SaveReconArray(ReconArray* recon, char* filename, int* status)
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
