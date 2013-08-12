#include "comarecon.h"

/////////////////////////////////////////////////////////////
//RECONSTRUCTION: The detected photons are                 //
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
  ReadEvent* ear=NULL;
  double* ReconArray1d=NULL; 
  double* ReconImage1d=NULL;
  double* EventArray1d=NULL; 
  double* EventImage1d=NULL;
  /*double* BalancingArray1d=NULL;
    double* BalancingImage1d=NULL;*/

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

      if(par.RePixSize==0.){
	struct SourceImageParameters sip = {
	    .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth)-1,
	    .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth)-1,
	    .crpix1 = (float)(mask->naxis1*mask->cdelt1)/detector_pixels->xpixelwidth,
	    .crpix2 = (float)(mask->naxis2*mask->cdelt2)/detector_pixels->ypixelwidth,

	    .cdelt1 = delta,
	    .cdelt2 = delta,
	    .crval1 = 0.,
	    .crval2 = 0.
	  };

       sky_pixels=getEmptySourceImage(&sip, &status);
       CHECK_STATUS_BREAK(status);
      }

       int ii, jj; //counts

     ReconArray* recon=NULL;
     //thought to contain contribution from non-sensitive detector area; see mail Yuri
     //TODO: check, if needed
     /*BalancingArray* balance=NULL;*/

    //Get empty event array object (type: ReadEvent)
     ea=getEventArray(mask, detector_pixels, &status);

    //detector size <= mask size
    //if mask size = det size -> shift is zero
      int xdiff, ydiff; //difference in size of mask and det plane
      xdiff=((ea->naxis1)/2-detector_pixels->xwidth)/2;
      ydiff=((ea->naxis2)/2-detector_pixels->ywidth)/2;


    // Loop over all events in the FITS file.
    while (0==EventFileEOF(&eventfile->generic)) {

      status=readEventList_nextRow(eventfile, ea);
      CHECK_STATUS_BREAK(status);

      //Get the 2d-EventArray, padded to the upper right corner

      //minus 1:because if axis is of lenght N then pixels: 0 to N-1
      ea->EventArray[ea->rawx+ea->naxis1/2-1+xdiff][ea->rawy+ea->naxis2/2-1+ydiff]+=ea->charge;

    } // END of scanning the event list.
    CHECK_STATUS_BREAK(status);

    if(par.RePixSize!=0.){//only works, if new smaller pixel fit without remainder in former big ones!
      //EventArray has been built with original pix-size, in order to distribute the events correctly
      //from now on the newe smaller size is used -> also for the later called ReconArray
      double pixelwidth_big=detector_pixels->xpixelwidth;//width of the former EventArray pixels (the real det-pix-size)
      detector_pixels->xpixelwidth=par.RePixSize;
      detector_pixels->ypixelwidth=par.RePixSize;

      //get source image with correct axes, since xpixelwidth has changed
      float delta = atan(par.RePixSize/par.MaskDistance);
      struct SourceImageParameters sip = {
	    .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth)-1,
	    .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth)-1,
	    //due to repix the axes get even numbers -> have to be shifted half a pixel
	    //since one former pixel (for MIRAX) corresponds to 6 now -> shift: 4.0
	    //TODO: determine crpix value automatically, depending on even/odd axes
	    .crpix1 = (float)(mask->naxis1*mask->cdelt1)/detector_pixels->xpixelwidth-4.,
	    .crpix2 = (float)(mask->naxis2*mask->cdelt2)/detector_pixels->ypixelwidth-4.,

	    .cdelt1 = delta,
	    .cdelt2 = delta,
	    .crval1 = 0.,
	    .crval2 = 0.
	  };

       sky_pixels=getEmptySourceImage(&sip, &status);
       CHECK_STATUS_BREAK(status);

      //new pointer for EventArray with more entries, since smaller pix-size (EventArrayRepix)
      ear=getEventArray(mask, detector_pixels, &status);

      int xcount, ycount;               //bigcount
      int xpixelcount=0, ypixelcount=0; //smallcount
      double leftsmall=0.,leftbig=0.;        //left border of small and big pixel
      double topsmall=0.,topbig=0.;          //top border of small and big pixel

      //Scanning over all EventArray-elements to get EventArray with smaller pixel-size
     for(ycount=0; ycount<ea->naxis2;ycount++){
       for(xcount=0; xcount<ea->naxis1;xcount++){

	 topbig=ycount*pixelwidth_big;   //top of current big pixel
	 ypixelcount=ceil(topbig/detector_pixels->ypixelwidth);  //count for small pixel (new pixels in EventArray)
	 //current y-pix: top border of big pix/width of one small pix->determines 1st small in current big

	 do{//as long as in current big pixel in y-direction
	   topsmall=ypixelcount*detector_pixels->ypixelwidth; //top border of small pix: current small pix*width of one
	   
	   leftbig=xcount*pixelwidth_big;
	   xpixelcount=ceil(leftbig/detector_pixels->xpixelwidth);
	   do{//as long as in current big pixel in x-direction
	     leftsmall=xpixelcount*detector_pixels->xpixelwidth;
    
	     ear->EventArray[xpixelcount][ypixelcount]=ea->EventArray[xcount][ycount];	     

	     xpixelcount++;
	   }while(leftsmall+detector_pixels->xpixelwidth < (leftbig+pixelwidth_big));
	   //end current big pixel x-direction
	   ypixelcount++;
	 }while(topsmall+detector_pixels->ypixelwidth < (topbig+pixelwidth_big));
	 //end current big pixel y-direction


       }
     }
      
    }//end re-pixel EventArray to smaller size given by RePixSize


    //Get the reconstruction array:
    recon=getReconArray(mask, detector_pixels, &status);
    int Size1 = recon->naxis1;
    int Size2 = recon->naxis2;

    //Get the 1d image of the reconstruction array -> needed by FFTW
    ReconImage1d=SaveReconArray1d(recon, &status);
    //create FITS image of 1d ReconArray for testing
    testFitsImage1d(ReconImage1d, "R1dTest.fits", Size1, Size2, &status);

    if(par.RePixSize!=0.){
    //Get the 1d image of the event array -> needed by FFTW
    EventImage1d=SaveEventArray1d(ear, &status);
    testFitsImage1d(EventImage1d, "E1dTest.fits", Size1, Size2, &status);
    //Check whether the ReconArray and the EventArray have the same size
    if ((recon->naxis1 != ear->naxis1) || (recon->naxis2 != ear->naxis2)){
      printf ("Error: ReconArrray and EventArray must have the same size!\n");
      break;
    }
		  
    }else{
    //Get the 1d image of the event array -> needed by FFTW
    EventImage1d=SaveEventArray1d(ea, &status);
    //create FITS image of 1d EventArray for testing
    testFitsImage1d(EventImage1d, "E1dTest.fits", Size1, Size2, &status);
    //Check whether the ReconArray and the EventArray have the same size
    if ((recon->naxis1 != ea->naxis1) || (recon->naxis2 != ea->naxis2)){
      printf ("Error: ReconArrray and EventArray must have the same size!\n");
      break;
    }
   
    }
    
    //reconstruct sky-image via FFT

    //perform a fft with the ReconArray
    fftw_complex* fftReconArray=NULL;
    fftReconArray=FFTOfArray_1d(ReconImage1d, Size1, Size2, -1);

    //perform a fft with the EventArray
    fftw_complex* fftEventArray=NULL;
    fftEventArray=FFTOfArray_1d(EventImage1d, Size1, Size2, -1);
    
    //multiply fftEventArray with komplex conjugate of fftReconArray
    //Re-part: E(Re)*R(Re)+E(Im)*R(Im); Im-part: E(Re)*R(Im)-E(Im)*R(Re)
    fftw_complex* Multiply = NULL;
    Multiply=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Size1*Size2));
    for(ii=0; ii<Size1; ii++){
      for(jj=0; jj<Size2; jj++){
	Multiply[ii+Size1*jj][0]=fftEventArray[ii+Size1*jj][0]*fftReconArray[ii+Size1*jj][0]
	  +fftEventArray[ii+Size1*jj][1]*(fftReconArray[ii+Size1*jj][1]);
	Multiply[ii+Size1*jj][1]=-fftEventArray[ii+Size1*jj][0]*(fftReconArray[ii+Size1*jj][1])
	+fftEventArray[ii+Size1*jj][1]*fftReconArray[ii+Size1*jj][0];
      }
    }

    //Inverse FFT of Multilpy which already is of type fftw_complex
    fftw_complex* fftInvMultiply=NULL;
    fftInvMultiply=FFTOfArray(Multiply, Size1, Size2, +1);

    /*
    //get the balancing array
    balance=getBalancingArray(recon, detector_pixels, &eventfile->generic, &status);
    //testing:Get the 1d image of the balancing array:
     BalancingImage1d=SaveBalancingArray1d(balance, &status);
     testFitsImage1d(BalancingImage1d, "B1dTest.fits", Size1, Size2, &status);*/

    //save real part of inverse fft in sky image
     for(ii=0; ii<Size1-1; ii++){
      for(jj=0; jj<Size2-1; jj++){
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

  if(par.RePixSize!=0.){
    FreeEventArray(ear);}
  /*FreeBalancingArray(balance);*/
  fftw_free(fftReconArray);
  fftw_free(fftEventArray);
  fftw_free(fftInvMultiply);


  free_SourceImage(sky_pixels);

  } while(0);  // END of the error handling loop.

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

  // Read the width of the detector in [pixel].
  else if ((status=PILGetInt("Width", &par->width))) {
    SIXT_ERROR("failed reading the detector width");
  }

  // Read the width of one detector pixel in [m].
  else if ((status=PILGetReal("Pixelwidth", &par->pixelwidth))) {
    SIXT_ERROR("failed reading the width of the detector pixels");
  }

  // Read the width of one re-pixeled detector pixel in [m].
  else if ((status=PILGetReal("RePixSize", &par->RePixSize))) {
    SIXT_ERROR("failed reading the width of the re-pixeled detector pixel");
  }

  // Read the distance between the coded mask and the detector plane [m].
  else if ((status=PILGetReal("MaskDistance", &par->MaskDistance))) {
    SIXT_ERROR("failed reading the distance between the mask and the detector");
  }
  CHECK_STATUS_RET(status, status);

  return(status);
}
