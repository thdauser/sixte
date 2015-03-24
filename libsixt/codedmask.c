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


   Copyright 2007-2014 Christian Schmid, Mirjam Oertel, FAU
*/

#include "codedmask.h"


CodedMask* getCodedMask(int* const status)
{
  CodedMask* mask=(CodedMask*)malloc(sizeof(CodedMask));
  if (NULL==mask) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Coded Mask!\n",
		   *status);
    return(mask);
  }

  // Set initial default values.
  mask->map=NULL;
  mask->transparent_pixels=NULL;
  mask->n_transparent_pixels=0;
  mask->transparency=0.;
  mask->naxis1=0;
  mask->naxis2=0;
  mask->cdelt1=0.;
  mask->cdelt2=0.;
  mask->crpix1=0.;
  mask->crpix2=0.;
  mask->crval1=0.;
  mask->crval2=0.;

  return(mask);
}


void destroyCodedMask(CodedMask** const mask)
{
  if (NULL!=(*mask)) {
    if (NULL!=(*mask)->map) {
      int count;
      for(count=0; count<(*mask)->naxis1; count++) {
	if (NULL!=(*mask)->map[count]) {
	  free((*mask)->map[count]);
	}
      }
      free((*mask)->map);
    }
    if (NULL!=(*mask)->transparent_pixels) {
      int count;
      for(count=0; count<(*mask)->n_transparent_pixels; count++) {
	if (NULL!=(*mask)->transparent_pixels[count]) {
	  free((*mask)->transparent_pixels[count]);
	}
      }
      free((*mask)->transparent_pixels);
    }
    free(*mask);
    *mask=NULL;
  }
}


CodedMask* getCodedMaskFromFile(const char* const filename, int* const status)
{
  // Obtain a new (empty) CodedMask object using the basic constructor.
  CodedMask* mask = getCodedMask(status);
  if (EXIT_SUCCESS!=*status) return(mask);

  fitsfile* fptr=NULL;
  int* input_buffer=NULL; // Input buffer for FITS image.
  int x, y;

  do { // Beginning of Error Handling loop. 

    // Open the FITS file containing the coded mask.
    headas_chat(5, "open Coded Mask FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, status)) break;

    // Determine the width of the image.
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, status)) break;
    assert(naxes[0]==naxes[1]);
    mask->naxis1 = (int)naxes[0];
    mask->naxis2 = (int)naxes[1];
    
    // Determine the width of one image pixel.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &mask->cdelt1, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &mask->cdelt2, comment, status))
      break;
    // Determine the WCS keywords of the image.
    if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &mask->crpix1, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &mask->crpix2, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &mask->crval1, comment, status))
      break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &mask->crval2, comment, status))
      break;

    // Allocate memory for the pixels of the image:
    mask->map = (int**)malloc(mask->naxis1*sizeof(int*));
    if (NULL!=mask->map) {
      for(x=0; x<mask->naxis1; x++) {
	mask->map[x] = (int*)malloc(mask->naxis2*sizeof(int));
	if(NULL==mask->map[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "CodedMask!\n", *status);
	  break;
	}
      }
      if (EXIT_SUCCESS!=*status) break;
    } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "CodedMask!\n", *status);
      break;
    }     
    // Allocate memory for input buffer (1D array):
    input_buffer=(int*)malloc(mask->naxis1*mask->naxis2*sizeof(int));
    if(NULL==input_buffer) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "CodedMask!\n", *status);
      break;
    }
    // END of memory allocation

    
    // READ the FITS image:
    int anynul=0;
    int null_value=0.;
    long fpixel[2] = {1, 1}; // lower left corner
    //                |--|--> FITS coordinates start at (1,1)
    long lpixel[2] = {mask->naxis1, mask->naxis2}; // upper right corner
    long inc[2] = {1, 1};
    if (fits_read_subset(fptr, TINT, fpixel, lpixel, inc, &null_value, 
			 input_buffer, &anynul, status)) break;

    // Copy the image from the input buffer to the map array.
    mask->n_transparent_pixels=0;
    for(x=0; x<mask->naxis1; x++) {
      for(y=0; y<mask->naxis2; y++) {
	mask->map[x][y] = input_buffer[x+ mask->naxis1*y];
	
	// Check if the pixel is transparent or opaque.
	if (TRANSPARENT==mask->map[x][y]) {
	  mask->n_transparent_pixels++;
	}
      }
    }


    // Determine the masks transparency.
    mask->transparency = 
      ((double)mask->n_transparent_pixels)/((double)(mask->naxis1*mask->naxis2));
    
    // Allocate memory for the list of transparent pixels.
    mask->transparent_pixels = (int**)malloc(mask->n_transparent_pixels*sizeof(int*));
    if (NULL==mask->transparent_pixels) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "CodedMask!\n", *status);
      break;
    }
    for(x=0; x<mask->n_transparent_pixels; x++) {
      mask->transparent_pixels[x] = (int*)malloc(2*sizeof(int));
      if (NULL==mask->transparent_pixels[x]) {
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: could not allocate memory to store the "
		       "CodedMask!\n", *status);
	break;
      }
    }
    if (EXIT_SUCCESS!=*status) break;

    int count=0;
    for (x=0; x<mask->naxis1; x++) {
      for (y=0; y<mask->naxis2; y++) {
	if (TRANSPARENT==mask->map[x][y]) {
	  mask->transparent_pixels[count][0] = x;
	  mask->transparent_pixels[count][1] = y;
	  count++;
	}
      }
    }  

  } while (0); // END of Error Handling loop.

  // Release memory from image input buffer.
  if (NULL!=input_buffer) free(input_buffer);

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  return(mask);
}


//impact position via projection downwards,'per hand'-> for comparison
int getImpactPos (struct Point2d* const position,
		  const Vector* const phodir,
		  const CodedMask* const mask, 
		  const struct Telescope* const telescope,
		  const Vector* const nz,
		  const float distance,
		  const float x_det,
		  const float y_det,
		  int* const status)
{// Check if a CodedMask is specified. If not, break.
  if (NULL==mask) return(0);

  // Get a random number.
  double rand=sixt_get_random_number(status);
    CHECK_STATUS_RET(*status, 0); 

   //Photon passes through transparent pixel.
   //Determine its impact position on the detection plane.

   //First:Determine the pixel(random)the photon passes through (mask plane).
  int pixel = (int)(rand * mask->n_transparent_pixels);
   //Pixel points to an arbitrary pixel out of all transparent ones.

   //Second: get the impact positon in the mask plane in meters.
   //(pixel_pos - reference_pixel  +0.5[since crpix is the middle of the ref. pixel,
   //so to start at the left pixel border, add 0.5] + random number[to have a random
   //position within the whole pixel])*width of one pixel+value at reference_pixel
   position->x =
     ((double)(mask->transparent_pixels[pixel][0])-mask->crpix1+0.5+sixt_get_random_number(status))
     *mask->cdelt1+mask->crval1;
   CHECK_STATUS_RET(*status, 0);
   position->y =
     ((double)(mask->transparent_pixels[pixel][1])-mask->crpix2+0.5+sixt_get_random_number(status))
     *mask->cdelt2+mask->crval2;
   CHECK_STATUS_RET(*status, 0);

   
   //third:if off-axis photon
   if(phodir->x!=nz->x || phodir->y!=nz->y || phodir->z!=nz->z){
   //Determine the components of the photon direction with respect to the
   //detector coordinate axes nx, ny (in mask plane).
   double x_comp = scalar_product(phodir, &telescope->nx);
   double y_comp = scalar_product(phodir, &telescope->ny);

   //Determine the component of the photon direction within the mask plane.
   double radius = sqrt(pow(x_comp,2.)+pow(y_comp,2.));
   //And the azimuthal angle (with respect to the nx-axis)
   double alpha=atan2(y_comp, x_comp);

   //Finally get the impact position in the detection plane.
   //Shift the above according to the off-axis position and
   //the distance mask-detector.
   position->x -= cos(alpha) * radius * distance;
   position->y -= sin(alpha) * radius * distance;
   }

   //shift origin from center to back left corner (if pos x-dir goes from back to front
   // and pos y-dir from left to right)
   position->x += x_det/2.;
   position->y += y_det/2.;

   //ensure that photon does not hit the walls
   if (position->x > x_det || position->x < 0. || position->y > y_det || position->y < 0.){
     return(0);
   }

return(1);
}


//impact position via wcs-projection
int getImpactPos_wcs (struct wcsprm* wcs, struct Point2d* const position, const CodedMask* const mask,
		   double const photon_ra, double const photon_dec, float const det_pixelwidth,
		   const float x_det, const float y_det, int* const status)
{
  // Check if a CodedMask is specified. If not, break.
  if (NULL==mask) return(0);

  double pixcrd[2];
  double imgcrd[2];
  double world[2]={photon_ra, photon_dec};
  double phi, theta;
  int stat=0;

  wcss2p(wcs,1,2,world,&phi,&theta,imgcrd,pixcrd,&stat);
  if(0!=stat){
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("wcs projection failed!\n", *status);
  }

  // Get a random number.
  double rand=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0); 

  //Photon passes through transparent pixel.
  //Determine its impact position on the detection plane.

  //Determine the pixel(random)the photon passes through (mask plane).
   int pixel =(int)(rand * mask->n_transparent_pixels);
  //Pixel points to an arbitrary pixel out of all transparent ones.
  //use current mask pixel - mask->crpix -> set coordinate system to center
  //above*mask->cdelt -> in meters
   position->x=((double)(mask->transparent_pixels[pixel][0])-mask->crpix1+0.5+sixt_get_random_number(status))
    *mask->cdelt1+mask->crval1+(double)((pixcrd[0]-wcs->crpix[0]+0.5)*det_pixelwidth);
   position->y=((double)(mask->transparent_pixels[pixel][1])-mask->crpix2+0.5+sixt_get_random_number(status))
    *mask->cdelt2+mask->crval2+(double)((pixcrd[1]-wcs->crpix[1]+0.5)*det_pixelwidth);

  position->x+=x_det/2;
  position->y+=y_det/2;

  //ensure that photon does not hit the walls
  if (position->x > x_det || position->x < 0. || position->y > y_det || position->y < 0.){
    return(0);
  }

  return(1);
}


int getImpactPos_protoMirax(struct wcsprm* wcs, struct wcsprm* wcs2, struct Point2d* const position,
			   const CodedMask* const mask, double const photon_ra, double const photon_dec,
			   float const det_pixelwidth, const float det_width, const float x_det,
			   const float y_det, const float wall, int* const status)
{
  // Check if a CodedMask is specified. If not, break.
  if (NULL==mask) return(0);

  double pixcrd[2];
  double imgcrd[2];
  double world[2]={photon_ra, photon_dec};
  double phi, theta;
  int stat=0;

  wcss2p(wcs,1,2,world,&phi,&theta,imgcrd,pixcrd,&stat);
  if(0!=stat){
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("wcs projection failed!\n", *status);
  }

  double pixcrd2[2];
  double imgcrd2[2];
  double world2[2]={photon_ra, photon_dec};
  double phi2, theta2;
  int stat2=0;

  wcss2p(wcs2,1,2,world2,&phi2,&theta2,imgcrd2,pixcrd2,&stat2);
  if(0!=stat2){
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("wcs projection failed!\n", *status);
  }

  // Get a random number.
  double rand=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0);

  //Photon passes through transparent pixel.
  //Determine its impact position on the detection plane.

  //Determine the pixel(random)the photon passes through (mask plane).
  int pixel =(int)(rand * mask->n_transparent_pixels);
  //Pixel points to an arbitrary pixel out of all transparent ones.
  
  //impact position at collimator distance (between mask and det-plane) in meters
  double x_coll=((double)(mask->transparent_pixels[pixel][0])-mask->crpix1+0.5+sixt_get_random_number(status))
    *mask->cdelt1+mask->crval1+(double)((pixcrd2[0]-wcs2->crpix[0])*det_pixelwidth);
  double y_coll=((double)(mask->transparent_pixels[pixel][1])-mask->crpix2+0.5+sixt_get_random_number(status))
    *mask->cdelt2+mask->crval2+(double)((pixcrd2[1]-wcs2->crpix[1])*det_pixelwidth);

  //shift origin from center to corner
  x_coll+=x_det/2;
  y_coll+=y_det/2;

  //determine detector-element that is hit at collimator distance
  int x_coll_det=(int)(x_coll/det_width);
  int y_coll_det=(int)(y_coll/det_width);

  //impact position at detection plane in meters
  position->x=((double)(mask->transparent_pixels[pixel][0])-mask->crpix1+0.5+sixt_get_random_number(status))
    *mask->cdelt1+mask->crval1+(double)((pixcrd[0]-wcs->crpix[0]+0.5)*det_pixelwidth);
  position->y=((double)(mask->transparent_pixels[pixel][1])-mask->crpix2+0.5+sixt_get_random_number(status))
    *mask->cdelt2+mask->crval2+(double)((pixcrd[1]-wcs->crpix[1]+0.5)*det_pixelwidth);

  //shift origin from center to corner
  position->x+=x_det/2;
  position->y+=y_det/2;

  //determine detector that is hit at detection plane
  int x_detector=(int)(position->x/det_width);
  int y_detector=(int)(position->y/det_width);

  //ensure that photon does not hit the outer walls
  if(position->x > x_det || position->x < 0. || position->y > y_det || position->y < 0.){
    return(0);
  }

  //ensure that photon does not travel through collimator walls and hit another detector
  if(x_coll_det != x_detector || y_coll_det != y_detector){
    return(0);
  }

  //switch to coordinate sys of one detector (that is surrounded by the collimator)
  position->x=position->x-(x_detector*det_width);  //position [m] [0:det_width]
  position->y=position->y-(y_detector*det_width);

  if(position->x <= wall || position->x > (det_width-wall) || position->y <= wall || position->y > (det_width-wall)){
    return(0);
    }

  return(1);
}
