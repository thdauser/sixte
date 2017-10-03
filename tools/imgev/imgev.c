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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "teseventlist.h"

#define TOOLSUB imgev_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EvtFile[MAXFILENAME];
  char Image[MAXFILENAME];

  int coordinatesystem;
  char projection[MAXMSG];
  long naxis1, naxis2;
  char cunit1[MAXMSG], cunit2[MAXMSG];
  float crval1, crval2;
  float crpix1, crpix2;
  float cdelt1, cdelt2;
  
  char clobber;
};


static int imgev_getpar(struct Parameters *par);


////////////////////////////////////
/** Main procedure. */
int imgev_main() {
  // Program parameters.
  struct Parameters par; 

  // Input event file.
  EventFile* elf=NULL;

  // Input TES event file.
  TesEventFile* tes_elf=NULL;

  // Input file (either TES or normal event file).
  fitsfile* input_fptr=NULL;

  // Output image.
  long** img=NULL;
  struct wcsprm wcs={ .flag=-1 };

  // FITS file access.
  long* img1d=NULL;
  char* headerstr=NULL;
  fitsfile* imgfptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS;   


  // Register HEATOOL:
  set_toolname("imgev");
  set_toolversion("0.03");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=imgev_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Check if a valid coordinate system has been selected.
    if ((par.coordinatesystem<0)||(par.coordinatesystem>1)) {
      status=EXIT_FAILURE;
      SIXT_ERROR("invalid selection for coordinate system");
      break;
    }
    
    headas_chat(3, "initialize ...\n");

    // Set the event file.
    long nrows = 0;
    elf=openEventFile(par.EvtFile, READWRITE, &status);
    if(status==COL_NOT_FOUND){
      headas_chat(3, "Given file is not a standard Event File, trying to read it as TES Event File...\n");
      status=EXIT_SUCCESS;
      freeEventFile(&elf, &status);
      CHECK_STATUS_BREAK(status);
      tes_elf=openTesEventFile(par.EvtFile,READWRITE,&status);
      nrows = tes_elf->nrows;
      input_fptr = tes_elf->fptr;
    } else {
      nrows = elf->nrows;
      input_fptr = elf->fptr;
    }
    CHECK_STATUS_BREAK(status);

    // Allocate memory for the output image.
    img=(long**)malloc(par.naxis1*sizeof(long*));
    CHECK_NULL_BREAK(img, status, "memory allocation failed");
    long ii; 
    for (ii=0; ii<par.naxis1; ii++) {
      img[ii]=(long*)malloc(par.naxis2*sizeof(long));
      CHECK_NULL_BREAK(img[ii], status, "memory allocation failed");
      long jj;
      for (jj=0; jj<par.naxis2; jj++) {
	img[ii][jj]=0;
      }
    }

    // Determine the projection type.
    char ctype1[MAXMSG], ctype2[MAXMSG];
    if (0==par.coordinatesystem) {
      strcpy(ctype1, "RA---");
      strcpy(ctype2, "DEC--");
    } else if (1==par.coordinatesystem) {
      strcpy(ctype1, "GLON-");
      strcpy(ctype2, "GLAT-");
    }
    strcat(ctype1, par.projection);
    strcat(ctype2, par.projection);

    if (strlen(ctype1)!=8) {
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "invalid projection type: CTYPE1='%s'", ctype1);
      SIXT_ERROR(msg);
      break;
    }
    if (strlen(ctype2)!=8) {
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "invalid projection type: CTYPE2='%s'", ctype2);
      SIXT_ERROR(msg);
      break;
    }

    // Set up the WCS data structure.
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis=2;
    wcs.crpix[0]=par.crpix1;
    wcs.crpix[1]=par.crpix2;
    wcs.crval[0]=par.crval1;
    wcs.crval[1]=par.crval2;
    wcs.cdelt[0]=par.cdelt1;
    wcs.cdelt[1]=par.cdelt2;
    strcpy(wcs.cunit[0], par.cunit1);
    strcpy(wcs.cunit[1], par.cunit2);
    strcpy(wcs.ctype[0], ctype1);
    strcpy(wcs.ctype[1], ctype2);

    // --- END of Initialization ---


    // --- Beginning Image Binning ---

    headas_chat(5, "image binning ...\n");

    // LOOP over all events in the FITS table.
    long row;
    double ra,dec;
    int anynul=0;
    double dnull=0.;
    for (row=0; row<nrows; row++) {
      
      if (elf!=NULL){
	// Read the next event from the file.
	Event event;
	getEventFromFile(elf, row+1, &event, &status);
	CHECK_STATUS_BREAK(status);
	ra = event.ra;
	dec = event.dec;
      } else {
	fits_read_col(tes_elf->fptr, TDOUBLE, tes_elf->raCol, row+1, 1, 1,
		      &dnull, &ra, &anynul, &status);
	ra*=M_PI/180.;
	fits_read_col(tes_elf->fptr, TDOUBLE, tes_elf->decCol, row+1, 1, 1,
		      &dnull, &dec, &anynul, &status);
	dec*=M_PI/180.;
      }
      // Convert the coordinates to the desired coordinate system.
      double lon, lat; // [rad].
      if (0==par.coordinatesystem) {
    	  // Equatorial coordinates.
    	  lon=ra*180./M_PI;
    	  lat=dec*180./M_PI;
      } else {
    	  // Galactic coordinates.
    	  const double l_ncp=2.145566759798267518;
    	  const double ra_ngp=3.366033268750003918;
    	  const double cos_d_ngp=0.8899880874849542;
    	  const double sin_d_ngp=0.4559837761750669;
    	  double cos_d=cos(dec);
    	  double sin_d=sin(dec);
    	  lon=(l_ncp-atan2(cos_d*sin(ra-ra_ngp),
    			  cos_d_ngp*sin_d-sin_d_ngp*cos_d*cos(ra-ra_ngp)))
			  *180./M_PI;
    	  lat=asin(sin_d_ngp*sin_d + cos_d_ngp*cos_d*cos(ra-ra_ngp))*180./M_PI;
      }

      // Determine the image coordinates corresponding to the event.
      double pixcrd[2];
      double imgcrd[2];
      double world[2]={lon, lat};
      double phi, theta;
      int status2=0;
      wcss2p(&wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &status2);
      if (3==status2) { 
	// Pixel does not correspond to valid world coordinates.
	continue;
      } else if (0!=status2) {
	SIXT_ERROR("projection failed");
	status=EXIT_FAILURE;
	break;
      }
      
      // Increase the image value at the event position.
      long xx=((long)(pixcrd[0]+0.5))-1;
      long yy=((long)(pixcrd[1]+0.5))-1;
      if ((xx>=0)&&(xx<par.naxis1) && (yy>=0)&&(yy<par.naxis2)) {
	img[xx][yy]++;
      }
    }
    CHECK_STATUS_BREAK(status);

    // Store the image in the output file.
    // Convert it to a 1d-array to store it in the FITS image.
    img1d=(long*)malloc(par.naxis1*par.naxis2*sizeof(long));
    CHECK_NULL_BREAK(img1d, status, "memory allocation for 1d image failed");
    for (ii=0; ii<par.naxis1; ii++) {
      long jj;
      for (jj=0; jj<par.naxis2; jj++) {
	img1d[ii + jj*par.naxis1]=img[ii][jj];
      }
    }

    // Create a new FITS-file (remove existing one before):
    remove(par.Image);
    fits_create_file(&imgfptr, par.Image, &status);
    CHECK_STATUS_BREAK(status);

    // Create an image in the FITS-file (primary HDU):
    long naxes[2]={ par.naxis1, par.naxis2 };
    fits_create_img(imgfptr, LONG_IMG, 2, naxes, &status);
    //                                 |-> naxis
    CHECK_STATUS_BREAK(status);

    // Copy the mission header keywords.
    char comment[MAXMSG], telescop[MAXMSG]={""}, instrume[MAXMSG]={""};
    fits_read_key(input_fptr, TSTRING, "TELESCOP", &telescop, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(imgfptr, TSTRING, "TELESCOP", telescop, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_read_key(input_fptr, TSTRING, "INSTRUME", &instrume, comment, &status);
    CHECK_STATUS_BREAK(status);
    fits_update_key(imgfptr, TSTRING, "INSTRUME", instrume, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Write WCS header keywords.
    int nkeyrec;
    if (0!=wcshdo(0, &wcs, &nkeyrec, &headerstr)) {
      SIXT_ERROR("construction of WCS header failed");
      status=EXIT_FAILURE;
      break;
    }
    char* strptr=headerstr;
    while (strlen(strptr)>0) {
      char strbuffer[81];
      strncpy(strbuffer, strptr, 80);
      strbuffer[80]='\0';
      fits_write_record(imgfptr, strbuffer, &status);
      CHECK_STATUS_BREAK(status);
      strptr+=80;
    }
    CHECK_STATUS_BREAK(status);

    // Write the image to the file.
    long fpixel[2]={1, 1}; // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2]={par.naxis1, par.naxis2}; 
    fits_write_subset(imgfptr, TLONG, fpixel, lpixel, img1d, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Close the files.
  freeEventFile(&elf, &status);
  freeTesEventFile(tes_elf,&status);

  // Close image file.
  if (NULL!=imgfptr) fits_close_file(imgfptr, &status);

  // Release buffer memory.
  if (NULL!=headerstr) {
    free(headerstr);
  }

  // Free the image.
  if (NULL!=img) {
    long ii;
    for (ii=0; ii<par.naxis1; ii++) {
      if (NULL!=img[ii]) {
	free(img[ii]);
      }
    }
    free(img);
  }
  if (NULL!=img1d) free(img1d);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}


static int imgev_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("EvtFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  } 
  strcpy(par->EvtFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Image", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output image file");
    return(status);
  } 
  strcpy(par->Image, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("CoordinateSystem", &par->coordinatesystem);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading coordinate system");
    return(status);
  } 

  status=ape_trad_query_string("Projection", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading projection type");
    return(status);
  } 
  strcpy(par->projection, sbuffer);
  free(sbuffer);

  status=ape_trad_query_long("NAXIS1", &par->naxis1);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading NAXIS1");
    return(status);
  } 

  status=ape_trad_query_long("NAXIS2", &par->naxis2);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading NAXIS2");
    return(status);
  } 

  status=ape_trad_query_string("CUNIT1", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CUNIT1");
    return(status);
  } 
  strcpy(par->cunit1, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("CUNIT2", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CUNIT2");
    return(status);
  } 
  strcpy(par->cunit2, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("CRVAL1", &par->crval1);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CRVAL1");
    return(status);
  } 

  status=ape_trad_query_float("CRVAL2", &par->crval2);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CRVAL2");
    return(status);
  } 

  status=ape_trad_query_float("CRPIX1", &par->crpix1);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CRPIX1");
    return(status);
  } 

  status=ape_trad_query_float("CRPIX2", &par->crpix2);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CRPIX2");
    return(status);
  } 

  status=ape_trad_query_float("CDELT1", &par->cdelt1);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CDELT1");
    return(status);
  } 

  status=ape_trad_query_float("CDELT2", &par->cdelt2);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CDELT2");
    return(status);
  } 

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}



