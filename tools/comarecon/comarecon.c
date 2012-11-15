#include "comarecon.h"


////////////////////////////////////
/** Main procedure. */
int comarecon_main() {
  struct Parameters par;
  
  CoMaEventFile* eventfile=NULL;
  SquarePixels* detector_pixels=NULL;
  CodedMask* mask=NULL;
  SourceImage* sky_pixels=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comarecon");
  set_toolversion("0.02");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using the PIL library.
    if ((status=comarecon_getpar(&par))) break;
    
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

    // Perform the image reconstruction algorithm.
    // We calculate the correlation S = D * A, instead of S = D * G as
    // proposed by Gottesman and Fenimore 1989, because I currently don't see
    // the reason to use the proposed decoding function instead of the mask
    // function itself.
    int i,j, ishift, jshift, k,l; // Counters.
    for (i=0; i<sky_pixels->naxis1; i++) {
      ishift=i-sky_pixels->naxis1/2;
      for (j=0; j<sky_pixels->naxis2; j++) {
	jshift=j-sky_pixels->naxis2/2;

	for (k=MAX(0,-ishift); k<detector_pixels->xwidth-1-MAX(0,ishift); k++) {
	  for (l=MAX(0,-jshift); l<detector_pixels->ywidth-1-MAX(0,jshift); l++) {
	    sky_pixels->pixel[i][j] += 
	      detector_pixels->array[ishift+k][jshift+l].charge * mask->map[k][l];
	    // TODO Maybe the indices are wrong and we have to write
	    // 	  detector_pixels->array[k][l].charge * mask.map[ishift+k][jshift+l];
	    // instead.
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

  // Close the FITS files.
  status=closeCoMaEventFile(eventfile);

  if (EXIT_SUCCESS==status) headas_chat(5, "finished successfully!\n\n");
  return(status);
}


int comarecon_getpar(struct Parameters* par)
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

  // Get the filename of the image file (FITS output file).
  else if ((status=PILGetFname("Image", par->Image))) {
    SIXT_ERROR("failed reading the filename of the output image");
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



