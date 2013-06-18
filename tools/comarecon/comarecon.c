#include "comarecon.h"

/////////////////////////////////////////////////////////////
//simulates reconstruction process.The detected photons are//
//deconvolved via the ReconArray and FFT                   //
//Input:event-list(t,charge,RAWX,RAWY),mask-file,          //
//      det-width(in pixel), pixelwidth, distance          //
//Output:sky-image                                         //
/////////////////////////////////////////////////////////////


////////////////////////////////////
/** Main procedure. */
int comarecon_main() {
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
  double* BalancingArray1d=NULL;
  double* BalancingImage1d=NULL;

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("comarecon");
  set_toolversion("0.03");


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

    
    // END of SKY IMAGE CONFIGURATION SETUP    
    
    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    // Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start image reconstruction process ...\n");
    if (par.ReconType ==1)
      //Begin of reconstruction process type1 (direct deconvolution)
      {
	struct SourceImageParameters sip = {
	.naxis1 = 2*mask->naxis1-1,
	.naxis2 = 2*mask->naxis2-1,
	.crpix1 = mask->naxis1,
	.crpix2 = mask->naxis2,

	.cdelt1 = delta,
	.cdelt2 = delta,
	.crval1 = 0.,
	.crval2 = 0.
	 };
	 sky_pixels=getEmptySourceImage(&sip, &status);
	 CHECK_STATUS_BREAK(status);



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
      struct SourceImageParameters sip = {
	    .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth),
	    .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth),
	    .crpix1 = (float)(mask->naxis1*mask->cdelt1)/detector_pixels->xpixelwidth,
	    .crpix2 = (float)(mask->naxis2*mask->cdelt2)/detector_pixels->ypixelwidth,

	    .cdelt1 = delta,
	    .cdelt2 = delta,
	    .crval1 = 0.,
	    .crval2 = 0.
	  };

       sky_pixels=getEmptySourceImage(&sip, &status);
       CHECK_STATUS_BREAK(status);


     ReconArray* recon=NULL;
     BalancingArray* balance=NULL;

    //Get empty event array object (type: ReadEvent)
     ea=getEventArray(mask, detector_pixels, &status);

    //mask size >= detector size
    //if mask size = det size -> shift is zero
    //for testing:same pixelsize:
     int xshift, yshift;
     if(mask->cdelt1 == detector_pixels->xpixelwidth){
       xshift=(mask->naxis1-detector_pixels->xwidth)/2;
       yshift=(mask->naxis2-detector_pixels->ywidth)/2;
     }else{ //diff pixelsize
       xshift=((ea->naxis1)/2-detector_pixels->xwidth)/2;
       yshift=((ea->naxis2)/2-detector_pixels->ywidth)/2;
     }

    // Loop over all events in the FITS file.
    while (0==EventFileEOF(&eventfile->generic)) {

      status=readEventList_nextRow(eventfile, ea);
      CHECK_STATUS_BREAK(status);

      //Get the 2d-EventArray

      //for testing:
      /* ea->EventArray[ea->rawx+xshift][ea->rawy+yshift]+=ea->charge;*/
      /* ea->EventArray[ea->rawx+(ea->naxis1)/4+xshift][ea->rawy+(ea->naxis2)/4+yshift]+=ea->charge;*/

      //minus 1:because if axis is of lenght N then pixels: 0 to N-1
      ea->EventArray[ea->rawx+(ea->naxis1)/2-1+xshift][ea->rawy+(ea->naxis2)/2-1+yshift]+=ea->charge;
      
    } // END of scanning the event list.
    CHECK_STATUS_BREAK(status);

    //Get the reconstruction array:
    recon=getReconArray(mask, detector_pixels, &status);
    int Size1 = recon->naxis1;
    int Size2 = recon->naxis2;

    int ii, jj;


    //for testing
    // double* multiplyPointer=NULL;
    //multiplyPointer=MultiplyMaskRecon(recon, mask, &status);
    // testFitsImage1d(multiplyPointer, "multiplyMR.fits", Size1/2, Size2/2, &status);

    //Get the 1d image of the reconstruction array:
    ReconImage1d=SaveReconArray1d(recon, &status);
    //testFitsImage1d(ReconImage1d, "R1dTest.fits", Size1, Size2, &status);

    //Get the 1d image of the event array:
    EventImage1d=SaveEventArray1d(ea, &status);
    //testFitsImage1d(EventImage1d, "E1dTest.fits", Size1, Size2, &status);

    //Check whether the ReconArray and the EventArray have the same size
    if ((recon->naxis1 != ea->naxis1) || (recon->naxis2 != ea->naxis2)){
      printf ("Error: ReconArrray and EventArray must have the same size!\n");
      break;
    }
   
    //reconstruct sky-image via FFT
    //perform a fft with the ReconArray
   
    fftw_complex* fftReconArray=NULL;
    fftReconArray=FFTOfArray_1d(ReconImage1d, Size1, Size2, -1);
    testFitsImagefft_real(fftReconArray, "fftReconArray_real.fits", Size1, Size2);
    testFitsImagefft_img(fftReconArray, "fftReconArray_img.fits", Size1, Size2);

    //testing: inverse FFT of ReconArray
    /* fftw_complex* InvfftReconArray=NULL;
    fftw_complex* fftReconBuffer=NULL;
    InvfftReconArray=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Size1*Size2));
    fftReconBuffer=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Size1*Size2));*/

    /*for(ii=Size1/2; ii<=Size1; ii++){
      for(jj=0; jj<=Size2/2; jj++){
	fftReconBuffer[(Size1/2-ii)+Size1*(jj+Size2/2)][0]=fftReconArray[ii+Size1*jj][0];
	fftReconBuffer[(Size1/2-ii)+Size1*(jj+Size2/2)][1]=fftReconArray[ii+Size1*jj][1];
      }
      }

     for(ii=0; ii<Size1/2; ii++){
      for(jj=0; jj<=Size2/2; jj++){
	fftReconBuffer[(ii+Size1/2)+Size1*(jj+Size2/2)][0]=fftReconArray[ii+Size1*jj][0];
	fftReconBuffer[(ii+Size1/2)+Size1*(jj+Size2/2)][1]=fftReconArray[ii+Size1*jj][1];
      }
      }*/
     //InvfftReconArray=FFTOfArray(fftReconBuffer, Size1, Size2, 1);
    /* testFitsImagefft_real(fftReconBuffer, "fftReconBuffer_real.fits", Size1, Size2);
    testFitsImagefft_img(fftReconBuffer, "fftReconBuffer_img.fits", Size1, Size2);*/
    //end testing: inverse FFT of ReconArray


    //perform a ft with the EventArray
    fftw_complex* fftEventArray=NULL;
    fftEventArray=FFTOfArray_1d(EventImage1d, Size1, Size2, -1);
    //testFitsImagefft_real(fftEventArray, "fftEventArray_real.fits", Size1, Size2);
    //testFitsImagefft_img(fftEventArray, "fftEventArray_img.fits", Size1, Size2);

    //multiply the komplex conjugate of fftEventArray with fftReconArray
    fftw_complex* Multiply = NULL;
    Multiply=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Size1*Size2));
    for(ii=0; ii<Size1; ii++){
      for(jj=0; jj<Size2; jj++){
	Multiply[ii+Size1*jj][0]=fftEventArray[ii+Size1*jj][0]*fftReconArray[ii+Size1*jj][0]
	  -fftEventArray[ii+Size1*jj][1]*((-1)*fftReconArray[ii+Size1*jj][1]);
	Multiply[ii+Size1*jj][1]=fftEventArray[ii+Size1*jj][0]*((-1)*fftReconArray[ii+Size1*jj][1])
	+fftEventArray[ii+Size1*jj][1]*fftReconArray[ii+Size1*jj][0];
      }
    }

    //test Multiply
    //testFitsImagefft_real(Multiply, "fftMultiply_real.fits", Size1, Size2);
    //testFitsImagefft_img(Multiply, "fftMultiply_img.fits", Size1, Size2);


    //Inverse FFT of Multilpy which already is of type fftw_complex
    fftw_complex* fftInvMultiply=NULL;
    fftInvMultiply=FFTOfArray(Multiply, Size1, Size2, +1);
    //testFitsImagefft_real(fftInvMultiply, "InvMultiply_real.fits", Size1, Size2);
    //testFitsImagefft_img(fftInvMultiply, "InvMultiply_img.fits", Size1, Size2);

    /*
    //get the balancing array
    balance=getBalancingArray(recon, detector_pixels, &eventfile->generic, &status);
    //testing:Get the 1d image of the balancing array:
     BalancingImage1d=SaveBalancingArray1d(balance, &status);
     testFitsImage1d(BalancingImage1d, "B1dTest.fits", Size1, Size2, &status);*/



     for(ii=0; ii<Size1; ii++){
      for(jj=0; jj<Size2; jj++){
	sky_pixels->pixel[ii][jj]=fftInvMultiply[ii+Size1*jj][0]/(Size1*Size2);//-balance->Bmap[ii][jj];
	}
	}

    // Write the reconstructed source function to the output FITS file.
    saveSourceImage(sky_pixels, par.Image, &status);
    CHECK_STATUS_BREAK(status);
    

  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the detector and sky image pixels.
  destroySquarePixels(&detector_pixels);
  FreeReconArray(recon);
  FreeReconArray1d(ReconArray1d);
  FreeEventArray(ea);
  FreeEventArray1d(EventArray1d);
  FreeBalancingArray(balance);
  fftw_free(fftReconArray);
  fftw_free(fftEventArray);
  fftw_free(fftInvMultiply);


  free_SourceImage(sky_pixels);

  

    }//END of reconstruction process type2 (FFT)
  else {
    //Error: wrong reconstruction type
  }} while(0);  // END of the error handling loop.

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
