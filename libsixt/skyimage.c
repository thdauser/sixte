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


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include "skyimage.h"


SkyImage* get_SkyImage()
{
  SkyImage* si=NULL;

  // Allocate memory:
  si=(SkyImage*)malloc(sizeof(SkyImage));
  if(si!=NULL) {
    si->pixel=NULL;
    si->ra=NULL;
    si->dec=NULL;

    si->naxis1=0;
    si->naxis2=0;

    si->cdelt1 = 0.;
    si->cdelt2 = 0.;
    si->crpix1 = 0.;
    si->crpix2 = 0.;
    si->crval1 = 0.;
    si->crval2 = 0.;

    si->minra  = 0.;
    si->maxra  = 0.;
    si->mindec = 0.;
    si->maxdec = 0.;
  }

  return(si);
}


SkyImage* getEmptySkyImage(struct SkyImageParameters* sip, int* status)
{
  SkyImage* si=NULL;

  // Check if the requested dimensions are reasonable.
  if ((sip->naxis1<0) || (sip->naxis2<0)) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: array dimensions for SkyImage pixel array cannot "
		   "be negative!\n", *status);
    return(si);
  }

  // Obtain a bare SkyImage object from the standard constructor.
  si=get_SkyImage();
  if (NULL==si) return(si);

  si->naxis1 = sip->naxis1;
  si->naxis2 = sip->naxis2;

  si->cdelt1 = sip->cdelt1;
  si->cdelt2 = sip->cdelt2;
  si->crval1 = sip->crval1;
  si->crval2 = sip->crval2;
  si->crpix1 = sip->crpix1;
  si->crpix2 = sip->crpix2;

  si->minra  = sip->minra;
  si->maxra  = sip->maxra;
  si->mindec = sip->mindec;
  si->maxdec = sip->maxdec;

  // Allocate memory for pixel-array
  si->pixel=(double**)malloc(sip->naxis1*sizeof(double*));
  if (NULL==si->pixel) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory to store "
		   "the SkyImage!\n", *status);
    return(si);
  }
  int xcount, ycount;
  for(xcount=0; xcount<si->naxis1; xcount++) {
    si->pixel[xcount]=(double*)malloc(si->naxis2*sizeof(double));
    if (NULL==si->pixel[xcount]) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SkyImage!\n", *status);
      return(si);
    }
    // Clear the pixels.
    for(ycount=0; ycount<si->naxis2; ycount++) {
      si->pixel[xcount][ycount] = 0.;
    }

    // Allocate memory for ra-array
    si->ra=(double**)malloc(sip->naxis1*sizeof(double*));
    if (NULL==si->ra) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SkyImage RA-array!\n", *status);
      return(si);
    }
  }

  for(xcount=0; xcount<si->naxis1; xcount++) {
    si->ra[xcount]=(double*)malloc(si->naxis2*sizeof(double));
    if (NULL==si->ra[xcount]) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SkyImage RA-array!\n", *status);
      return(si);
    }
    // Clear the pixels.
    for(ycount=0; ycount<si->naxis2; ycount++) {
      si->ra[xcount][ycount] = 0.;
    }

    // Allocate memory for dec-array
    si->dec=(double**)malloc(sip->naxis1*sizeof(double*));
    if (NULL==si->dec) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SkyImage DEC-array!\n", *status);
      return(si);
    }
  }

  for(xcount=0; xcount<si->naxis1; xcount++) {
    si->dec[xcount]=(double*)malloc(si->naxis2*sizeof(double));
    if (NULL==si->dec[xcount]) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SkyImage DEC-array!\n", *status);
      return(si);
    }
    // Clear the pixels.
    for(ycount=0; ycount<si->naxis2; ycount++) {
      si->dec[xcount][ycount] = 0.;
    }
  }

  return(si);
}

void fillRaDecArrays(SkyImage* si)
{
  int xcount, ycount;
  float cdelt_ra,cdelt_dec;

  cdelt_ra=0;
  for(xcount=0; xcount<si->naxis1; xcount++){ //one ra
    cdelt_dec=0.;
    for(ycount=0; ycount<si->naxis2; ycount++){ //all dec-positions for this ra
      si->ra[xcount][ycount]=(double)(si->minra+cdelt_ra);
      si->dec[xcount][ycount]=(double)(si->mindec+cdelt_dec);
      //next dec:
      cdelt_dec+=si->cdelt2;
    }
    //next ra:
    cdelt_ra+=si->cdelt1;
  }

}

void free_SkyImage(SkyImage* si)
{
  if(si != NULL) {
    if((si->naxis1 > 0)&&(NULL!=si->pixel)) {
      int count;
      for(count=0; count<si->naxis1; count++) {
	if(si->pixel[count] != NULL) free(si->pixel[count]);
      }
      free(si->pixel);
    }

    if((si->naxis1 > 0)&&(NULL!=si->ra)) {
      int count;
      for(count=0; count<si->naxis1; count++) {
	if(si->ra[count] != NULL) free(si->ra[count]);
      }
      free(si->ra);
    }

    if((si->naxis1 > 0)&&(NULL!=si->dec)) {
      int count;
      for(count=0; count<si->naxis1; count++) {
	if(si->dec[count] != NULL) free(si->dec[count]);
      }
      free(si->dec);
    }

    free(si);
  }
}


void saveSkyImage(SkyImage* si, char* filename, int* status)
{
  fitsfile *fptr=NULL;
  double *image1d=NULL;

  // Print information to STDOUT.
  char msg[MAXMSG];
  sprintf(msg, "Store SkyImage in file '%s' ...\n", filename);
  headas_chat(5, msg);

  do { // ERROR handling loop


    // If the specified file already exists, remove the old version.
    remove(filename);

    // Create a new FITS-file:
    if (fits_create_file(&fptr, filename, status)) break;

    // Allocate memory for the 1-dimensional image buffer (required for
    // output to FITS file).
    image1d = (double*)malloc(si->naxis1*si->naxis2*sizeof(double));
    if (!image1d) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error allocating memory!\n", *status);
      break;
    }

    // Store the source image in the 1-dimensional buffer to handle it
    // to the FITS routine.
    int x, y;
    for (x=0; x<si->naxis1; x++) {
      for (y=0; y<si->naxis2; y++) {
	image1d[(x+ si->naxis1*y)] = si->pixel[x][y];
      }
    }

    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = {(long)(si->naxis1), (long)(si->naxis2)};
    if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status)) break;
    //                                   |-> naxis
    //    int hdutype;
    if (fits_movabs_hdu(fptr, 1, NULL, status)) break;


    // Write the header keywords for the SourceImage.
    if (fits_write_key(fptr, TSTRING, "HDUCLAS1", "Image", "", status)) break;
    if (fits_write_key(fptr, TSTRING, "CTYPE1", "RA---TAN",
		       "sky coordinate system", status)) break;
    if (fits_write_key(fptr, TSTRING, "CTYPE2", "DEC--TAN",
		       "sky coordinate system", status)) break;

    if (fits_write_key(fptr, TSTRING, "CUNIT1", "degree", "", status)) break;
    if (fits_write_key(fptr, TSTRING, "CUNIT2", "degree", "", status)) break;
    if (fits_write_key(fptr, TFLOAT, "CRPIX1", &si->crpix1,
                       "X axis reference pixel", status)) break;
    if (fits_write_key(fptr, TFLOAT, "CRPIX2", &si->crpix2,
    		       "Y axis reference pixel", status)) break;

    float fbuffer;
    fbuffer = si->crval1;
    if (fits_write_key(fptr, TFLOAT, "CRVAL1", &fbuffer,
		       "coord of X ref pixel", status)) break;
    fbuffer = si->crval2;
    if (fits_write_key(fptr, TFLOAT, "CRVAL2", &fbuffer,
    		       "coord of Y ref pixel", status)) break;
    fbuffer = si->cdelt1;
    if (fits_write_key(fptr, TFLOAT, "CDELT1", &fbuffer,
		       "X axis increment", status)) break;
    fbuffer = si->cdelt2;
    if (fits_write_key(fptr, TFLOAT, "CDELT2", &fbuffer,
		       "Y axis increment", status)) break;

    HDpar_stamp(fptr, 1, status);
    if (EXIT_SUCCESS!=*status) break;
    // END of writing header keywords.


    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1)
    // Upper right corner.
    long lpixel[2] = {si->naxis1, si->naxis2};
    fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, image1d, status);

  } while (0); // END of ERROR handling loop

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  if (NULL!=image1d) free(image1d);
}
