#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "patternfile.h"
#include "pattern.h"
#include "wcslib.h"

#define TOOLSUB imgev_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char PatternList[MAXFILENAME];
  char Image[MAXFILENAME];

  long naxis1, naxis2;
  char ctype1[MAXMSG], ctype2[MAXMSG];
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

  // Input pattern list file.
  PatternFile* plf=NULL;

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
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=imgev_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Set the input pattern file.
    plf=openPatternFile(par.PatternList, READWRITE, &status);
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

    // Set up the WCS data structure.
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis = 2;
    wcs.crpix[0] = par.crpix1;
    wcs.crpix[1] = par.crpix2;
    wcs.crval[0] = par.crval1;
    wcs.crval[1] = par.crval2;
    wcs.cdelt[0] = par.cdelt1;
    wcs.cdelt[1] = par.cdelt2;
    strcpy(wcs.cunit[0], par.cunit1);
    strcpy(wcs.cunit[1], par.cunit2);
    strcpy(wcs.ctype[0], par.ctype1);
    strcpy(wcs.ctype[1], par.ctype2);


    // --- END of Initialization ---


    // --- Beginning Image Binning ---
    headas_chat(5, "image binning ...\n");

    // LOOP over all patterns in the FITS table.
    long row;
    for (row=0; row<plf->nrows; row++) {
      
      // Read the next pattern from the file.
      Pattern pattern;
      getPatternFromFile(plf, row+1, &pattern, &status);
      CHECK_STATUS_BREAK(status);
      
      // Determine the image coordinates of the pattern.
      double pixcrd[2];
      double imgcrd[2];
      double world[2] = {pattern.ra*180./M_PI, pattern.dec*180./M_PI};
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
      
      // Increase the image value at the pattern position.
      long xx = ((long)(pixcrd[0]+0.5))-1;
      long yy = ((long)(pixcrd[1]+0.5))-1;
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
	img1d[ii + jj*par.naxis1] = img[ii][jj];
      }
    }

    // Create a new FITS-file (remove existing one before):
    remove(par.Image);
    fits_create_file(&imgfptr, par.Image, &status);
    CHECK_STATUS_BREAK(status);

    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = { par.naxis1, par.naxis2 };
    fits_create_img(imgfptr, LONG_IMG, 2, naxes, &status);
    //                                 |-> naxis
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
      strbuffer[80] = '\0';
      fits_write_record(imgfptr, strbuffer, &status);
      CHECK_STATUS_BREAK(status);
      strptr += 80;
    }
    CHECK_STATUS_BREAK(status);

    // Write the image to the file.
    long fpixel[2] = {1, 1}; // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1), NOT (0,0).
    // Upper right corner.
    long lpixel[2] = {par.naxis1, par.naxis2}; 
    fits_write_subset(imgfptr, TLONG, fpixel, lpixel, img1d, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Close the files.
  destroyPatternFile(&plf, &status);

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

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}


static int imgev_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pattern list file");
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("Image", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output image file");
    return(status);
  } 
  strcpy(par->Image, sbuffer);
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

  status=ape_trad_query_string("CTYPE1", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CTYPE1");
    return(status);
  } 
  strcpy(par->ctype1, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("CTYPE2", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading CTYPE2");
    return(status);
  } 
  strcpy(par->ctype2, sbuffer);
  free(sbuffer);

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



