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
 
 
  ReadEvent* ea=NULL;
  double* ReconArray1d=NULL; 
  double* ReconImage1d=NULL;
  double* EventArray1d=NULL; 
  double* EventImage1d=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comarecon2");
  set_toolversion("0.03");


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
      .naxis1 = 2*par.width, //-1,
      .naxis2 = 2*par.width, //-1,
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
    if (par.ReconType ==1)
      //Begin of reconstruction process type1 (direct deconvolution)
      {
      ReconArrayFits* recon_fits=NULL;

      // Loop over all events in the FITS file.
      while (0==EventFileEOF(&eventfile->generic)) {

      CoMaEvent event;
      status=CoMaEventFile_getNextRow(eventfile, &event);
      CHECK_STATUS_BREAK(status);

      // Add the event to the SquarePixels array.
      detector_pixels->array[event.rawx][event.rawy].charge+=1.0;

    } // END of scanning the event list.
    CHECK_STATUS_BREAK(status);

    //Get the reconstruction array:
    recon_fits=getReconArrayForFits(mask, &status);

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
	     recon_fits->Rmap[ii][jj] * detector_pixels->array[ii+xshift][jj+yshift].charge;
	  }
	}
      }
    }
    // Write the reconstructed source function to the output FITS file.
    saveSourceImage(sky_pixels, par.Image, &status);
    SaveReconArrayFits(recon_fits, par.ReconArray, &status);
    CHECK_STATUS_BREAK(status);

  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the detector and sky image pixels.
  destroySquarePixels(&detector_pixels);
  free_SourceImage(sky_pixels);
  FreeReconArrayFits(recon_fits);
  }//END of reconstruction process type1 (direct deconvolution)

  else if(par.ReconType==2)
    //Begin of reconstruction process type2 (FFT)
    {
     ReconArray* recon=NULL;

    //Get empty event array object (type: ReadEvent).
    ea=getEventArray(detector_pixels, &status);

    // Loop over all events in the FITS file.
    while (0==EventFileEOF(&eventfile->generic)) {

      status=readEventList_nextRow(eventfile, ea);
      CHECK_STATUS_BREAK(status);

      //Get the 2d-EventArray
      ea->EventArray[ea->rawx][ea->rawy]+=ea->charge;
      
      
    } // END of scanning the event list.
    CHECK_STATUS_BREAK(status);

    //Get the reconstruction array:
    recon=getReconArray(mask, &status);
    int Size1 = recon->naxis1;
    int Size2 = recon->naxis2;

    int ii, jj, count;

    //Get the 1d image of the reconstruction array:
    ReconImage1d=SaveReconArray1d(recon, &status);
    //testFitsImage1d(ReconImage1d, "reconTest.fits", Size1, Size2);

    //Get the 1d image of the event array:
    EventImage1d=SaveEventArray1d(ea, &status);
    //testFitsImage1d(EventImage1d, "eventTest.fits", Size1, Size2);

    //Check whether the ReconArray and the EventArray have the same size
    if ((recon->naxis1 != ea->naxis1) || (recon->naxis2 != ea->naxis2)){
      printf ("Error: ReconArrray and EventArray must have the same size!\n");
      break;
    }
   
    //reconstruct sky-image via FFT
    //perform a fft with the ReconArray
   
    fftw_complex* fftReconArray=NULL;
    fftReconArray=FFTOfArray_1d(ReconImage1d, Size1, Size2);

    //perform a fft with the EventArray
    fftw_complex* fftEventArray=NULL;
    fftEventArray=FFTOfArray_1d(EventImage1d, Size1, Size2 );

    //multiply the komplex conjugate of fftEventArray with fftReconArray
    fftw_complex* Multiply = NULL;
    Multiply=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Size1*Size2/*+Size1*/));
    for(count=0; count<(Size1+Size1*Size2);count++){
      Multiply[count][0]=fftEventArray[count][0]*fftReconArray[count][0]
	-fftEventArray[count][1]*fftReconArray[count][1];
      Multiply[count][1]=fftEventArray[count][0]*fftReconArray[count][1]
	+fftEventArray[count][1]*fftReconArray[count][0];
    }
    //Inverse FFT of Multilpy which already is of type fftw_complex
    fftw_complex* fftInvMultiply=NULL;
    fftInvMultiply=FFTOfArrayInverse_fftwcomplex(Multiply,Size1,Size2);
       
    for(ii=0; ii<Size1; ii++){
      for(jj=0; jj<Size2; jj++){
	//sky_pixels->pixel[ii][jj]=fftInvMultiply[ii+Size1*jj][0];
	recon->RImage[ii][jj]=fftInvMultiply[ii+Size1*jj][0];
	}
    }

    //Resize the image
     for(ii=0; ii<=Size1/2; ii++){
      for(jj=0; jj<=Size2/2; jj++){
	sky_pixels->pixel[Size1/2+ii][Size2/2+jj]=recon->RImage[ii][jj];
      }
      }

     for(ii=Size1/2; ii<=Size1; ii++){
      for(jj=Size2/2; jj<=Size2; jj++){
	sky_pixels->pixel[(int)fabs(Size1/2-ii)][(int)fabs(Size2/2-jj)]=recon->RImage[ii][jj];
      }
      }

      for(ii=Size1/2+1; ii<Size1; ii++){
      for(jj=0; jj<Size2/2; jj++){
	sky_pixels->pixel[(int)fabs(Size1/2-ii)-1][Size2/2+jj+1]=recon->RImage[ii][jj];
      }
      }

     for(ii=0; ii<Size1/2; ii++){
       for(jj=Size2/2+1; jj<Size2; jj++){
	 sky_pixels->pixel[ii+Size1/2+1][(int)fabs(Size2/2-jj)-1]=recon->RImage[ii][jj];
      }
      }

    // Write the reconstructed source function to the output FITS file.
    saveSourceImage(sky_pixels, par.Image, &status);
    CHECK_STATUS_BREAK(status);
    

  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the detector and sky image pixels.
  destroySquarePixels(&detector_pixels);
  free_SourceImage(sky_pixels);
  FreeReconArray(recon);
  FreeReconImage(recon);
  FreeReconArray1d(ReconArray1d);
  FreeEventArray(ea);
  FreeEventArray1d(EventArray1d);

    }//END of reconstruction process type2 (FFT)
  else {
    //Error: wrong reconstruction type
  }} while(0);  // END of the error handling loop.

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

  // Get reconstruction type.
  else if ((status=PILGetInt("ReconType", &par->ReconType))) {
    SIXT_ERROR("failed reading the reconstruction type");
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
