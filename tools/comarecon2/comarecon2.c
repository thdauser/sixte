#include "comarecon2.h"

/////////////////////////////////////////////////////////////
//simulates reconstruction process.The detected photons are//
//deconvolved via the ReconArray and FFT(TODO)             //
//Input:event-list(t,charge,RAWX,RAWY),mask-file,          //
//      det-width(in pixel), pixelwidth, distance          //
//Output:sky-image                                         //
/////////////////////////////////////////////////////////////


////////////////////////////////////
/** Main procedure. */
int comarecon2_main() {
  struct Parameters par;
  
  CoMaEventFile* eventfile=NULL;
  SquarePixels* detector_pixels=NULL;
  CodedMask* mask=NULL;
  SourceImage* sky_pixels=NULL;
  ReconArray* recon=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comarecon2");
  set_toolversion("0.02");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using the PIL library.
    if ((status=comarecon2_getpar(&par))) break;
    
    // Open the event file.
    eventfile=openCoMaEventFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Load the coded mask from the file.
    mask=getCodedMaskFromFile(par.Mask, &status);
    CHECK_STATUS_BREAK(status);
    
    // DETECTOR setup.
    struct SquarePixelsParameters spp = {
      .xwidth = par.width,
      .ywidth = par.width,
      .xpixelwidth = par.pixelwidth,
      .ypixelwidth = par.pixelwidth 
    };
    detector_pixels=newSquarePixels(&spp, &status);
    CHECK_STATUS_BREAK(status);
    // END of DETECTOR CONFIGURATION SETUP    

    // SKY IMAGE setup.
    float delta = atan(par.pixelwidth/par.MaskDistance);
    struct SourceImageParameters sip = {
      .naxis1 = 2*par.width -1,
      .naxis2 = 2*par.width -1,
      .cdelt1 = delta,
      .cdelt2 = delta,
      .crval1 = 0.,
      .crval2 = 0.,
      .crpix1 = par.width*1.,
      .crpix2 = par.width*1.
    };
    sky_pixels=getEmptySourceImage(&sip, &status);
    CHECK_STATUS_BREAK(status);
    // END of SKY IMAGE CONFIGURATION SETUP    
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start image reconstruction process ...\n");

    // Loop over all events in the FITS file.
    while (0==EventFileEOF(&eventfile->generic)) {

      CoMaEvent event;
      status=CoMaEventFile_getNextRow(eventfile, &event);
      CHECK_STATUS_BREAK(status);

      // Add the event to the SquarePixels array.
      detector_pixels->array[event.rawx][event.rawy].charge+=1.0;

    } // END of scanning the impact list.
    CHECK_STATUS_BREAK(status);

    //TODO: Write imaging algorithm

    //Get the reconstruction array:
    recon=getReconArray(mask, &status);
    SaveReconArray(recon, par.ReconArray, &status);

    int s_x, s_y, ii, jj;
    //determine the convolution source_image = R*D
    //two loops over the source-image (2*size_det-1)
    for(s_x=0; s_x < sky_pixels->naxis1; s_x++){
      for (s_y=0; s_y < sky_pixels->naxis2; s_y++){
	//two loops over the detection-area. Each positioning of the mask
	//over the detector is considered.Therefore the shift-parameter is needed.
	int xshift=s_x - sky_pixels->naxis1/2;
	int yshift=s_y - sky_pixels->naxis2/2;
	for (ii=MAX(0,-xshift); ii < detector_pixels->xwidth-1-MAX(0,xshift); ii++){
	  for(jj=MAX(0,-yshift); jj < detector_pixels->ywidth-1-MAX(0,yshift); jj++){
	    sky_pixels->pixel[s_x][s_y]+=
	     recon->Rmap[ii][jj] * detector_pixels->array[ii+xshift][jj+yshift].charge;
	  }
	}
      }
    }
    // Write the reconstructed source function to the output FITS file.
    saveSourceImage(sky_pixels, par.Image, &status);
    CHECK_STATUS_BREAK(status);

  } while(0);  // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the detector and sky image pixels.
  destroySquarePixels(&detector_pixels);
  free_SourceImage(sky_pixels);
  FreeReconArray(recon);

  // Close the FITS files.
  status=closeCoMaEventFile(eventfile);

  if (EXIT_SUCCESS==status) headas_chat(5, "finished successfully!\n\n");
  return(status);
}


int comarecon2_getpar(struct Parameters* par)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the filename of the mask reconstruction file (FITS input file).
  if ((status=PILGetFname("Mask", par->Mask))) {
    SIXT_ERROR("failed reading the filename of the mask reconstruction image");
  }
 
  // Get the filename of the event list file (FITS input file).
  else if ((status=PILGetFname("EventList", par->EventList))) {
    SIXT_ERROR("failed reading the filename of the event list");
  }

  //Get the filename of the image file (FITS output file).
   else if ((status=PILGetFname("Image", par->Image))) {
  SIXT_ERROR("failed reading the filename of the output image");
  }

 //Get the filename of the ReconArray file (FITS output file).
   else if ((status=PILGetFname("ReconArray", par->ReconArray))) {
  SIXT_ERROR("failed reading the filename of the ReconArray");
  }

  // Read the width of the detector in [pixel].
  else if ((status=PILGetInt("Width", &par->width))) {
    SIXT_ERROR("failed reading the detector width");
  }

  // Read the width of one detector pixel in [m].
  else if ((status=PILGetReal("Pixelwidth", &par->pixelwidth))) {
    SIXT_ERROR("failed reading the width of the detector pixels");
  }

  // Read the distance between the coded mask and the detector plane [m].
  else if ((status=PILGetReal("MaskDistance", &par->MaskDistance))) {
    SIXT_ERROR("failed reading the distance between the mask and the detector");
  }
  CHECK_STATUS_RET(status, status);

  return(status);
}


//function test

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
