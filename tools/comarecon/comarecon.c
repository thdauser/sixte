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

#include "comarecon.h"

//////////////////////////////////////////////////////////////////////////////////////////
//RECONSTRUCTION: The detected photons are deconvolved via the ReconArray and EventArray//
//                using FFT like: SkyImg=FFTInv(FFT(EA)*FFT(RAcomplexconjugate));       //
//                IROS-algorithm for identifying sources                                //
//Input:event-list(t,charge,RAWX,RAWY),mask-file,detector-setup(width of det in pixels, //
//      pixelwidth,distance,pointing,gaps,RePixSize),sigma(threshold for sources)       //
//Output:sky-image, PositionList                                                        //
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////
/** Main procedure. */
int comarecon_main() {
  struct Parameters par;
  
  CoMaEventFile* eventfile=NULL;
  SquarePixels* detector_pixels=NULL;
  CodedMask* mask=NULL; 
  SourceImage* sky_pixels=NULL;
  ReconArray* recon=NULL;
  MaskShadow* mask_shadow=NULL;
  PixPositionList* position_list=NULL;
  double* median_list=NULL; //temp array of all background pix for determination of median 
  ReadEvent* ea=NULL;
  ReadEvent* ear=NULL;
  double* ReconImage1d=NULL;
  double* EventImage1d=NULL;
  fftw_complex* fftReconArray=NULL;
  fftw_complex* fftEventArray=NULL;
  fftw_complex* Multiply = NULL;
  fftw_complex* fftInvMultiply=NULL;

  int status=EXIT_SUCCESS; // Error status.

  int ii, jj; //counts
  int xdiff, ydiff; //difference in size of mask and det plane
  int PixAmount;

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
    double ra=par.RA;
    double dec=par.DEC;
    double distance=par.MaskDistance;

    double xdetsize_beforeRepix=par.width; //in the case of Repix the detector_pixels are overwritten, 
    double ydetsize_beforeRepix=par.width; //but some functions need the old smaller size (of the det in pixels)
    double xpixelsize_beforeRepix=par.pixelwidth; //original detector-pixelsize
    double ypixelsize_beforeRepix=par.pixelwidth;

    struct SquarePixelsParameters spp = {
      .xwidth = par.width,
      .ywidth = par.width,
      .xpixelwidth = par.pixelwidth,
      .ypixelwidth = par.pixelwidth,
      .DCU_length=par.DCU_length,
      .DCU_gap=par.DCU_gap,
      .DCA_gap=par.DCA_gap
    };
    detector_pixels=newSquarePixels(&spp, &status);
    CHECK_STATUS_BREAK(status);
    // END of DETECTOR CONFIGURATION SETUP    

    // SKY IMAGE setup.
       if(par.RePixSize==0.){
	 PixAmount=4;

	 float delta = atan(par.pixelwidth/distance);
	 struct SourceImageParameters sip = {
	   .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth),
	   .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth),
	   .crpix1 = (float)(sip.naxis1/2.+1), //even axes:ends to .5;odd axes (not possible,*2):ends to .0 (therefore: axis/2+0.5)
	   .crpix2 = (float)(sip.naxis2/2.+1), //further +0.5,since left boarder of 1st pix in sky image is 0.5(FITS)
	   .cdelt1 = delta,
	   .cdelt2 = delta,
	   .crval1 = ra*M_PI/180.,
	   .crval2 = dec*M_PI/180.
	};
       sky_pixels=getEmptySourceImage(&sip, &status);
       CHECK_STATUS_BREAK(status);
       }else{//only works, if new smaller pixel fit without remainder in former big ones!
	 //EventArray has been built with original pix-size, in order to distribute the events correctly
	 //from now on, the new smaller size is used -> also for the later called ReconArray
	 
	 detector_pixels->xwidth=(detector_pixels->xwidth*detector_pixels->xpixelwidth)/par.RePixSize;
	 detector_pixels->ywidth=(detector_pixels->ywidth*detector_pixels->ypixelwidth)/par.RePixSize;
	 detector_pixels->xpixelwidth=par.RePixSize;
	 detector_pixels->ypixelwidth=par.RePixSize;
	 
	 PixAmount=24;

	 //get source image with correct axes, since xpixelwidth has changed
	 float delta = atan(par.RePixSize/distance);
	 struct SourceImageParameters sip = {
	   .naxis1 = 2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth),
	   .naxis2 = 2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth),
	   //due to repix the axes get even numbers -> have to be shifted half a pixel
	   //since one former pixel (for MIRAX) corresponds to 6 now -> shift: 4.0
	   .crpix1 = (float)(sip.naxis1/2.+1),
	   .crpix2 = (float)(sip.naxis2/2.+1),
	   .cdelt1 = delta,
	   .cdelt2 = delta,
	   .crval1 = ra*M_PI/180.,
	   .crval2 = dec*M_PI/180.
	 };
	sky_pixels=getEmptySourceImage(&sip, &status);
        CHECK_STATUS_BREAK(status);

	}
    // END of SKY IMAGE CONFIGURATION SETUP 

    //SOURCE-POSITION DETERMINATION:
    double pixval=0.;
    int threshold=0; //for first source threshold is '0' after that: '1' if still bright enough, '2' else
    //get empty PixPositionList structure (contains pointer to current PixPosition-element
    //and count for found sources)
    position_list=getPixPositionList(sky_pixels);
    //memory-allocation for median_list
    median_list=getMedian_list(sky_pixels, &status);

    //initialization of wcs parameter structure for deteremining ra/dec of source
    struct wcsprm wcs = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis=2;
    wcs.crpix[0]=sky_pixels->crpix1;
    wcs.crpix[1]=sky_pixels->crpix2; 
    wcs.crval[0]=sky_pixels->crval1*180./M_PI;
    wcs.crval[1]=sky_pixels->crval2*180./M_PI;
    wcs.cdelt[0]=sky_pixels->cdelt1*180./M_PI; //in deg
    wcs.cdelt[1]=sky_pixels->cdelt2*180./M_PI;

    //initialization of wcs parameter structure for getting mask shadow
    struct wcsprm wcs2 = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs2)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs2.naxis=2;
    wcs2.crpix[0]=(detector_pixels->xwidth)/2;
    wcs2.crpix[1]=(detector_pixels->ywidth)/2;
    wcs2.crval[0]=ra;  //in deg
    wcs2.crval[1]=dec;
    wcs2.cdelt[0]=atan(detector_pixels->xpixelwidth/distance)*180./M_PI; //in deg
    wcs2.cdelt[1]=atan(detector_pixels->ypixelwidth/distance)*180./M_PI;


    //telescope coordinate system     //TODO: USE ATTITUDE
       //telescope pointing direction
       //Vector nz=normalize_vector(unit_vector(ra*M_PI/180.0,dec* M_PI/180.0));
       //unit-vector in z-direction:
       //Vector vz = {0.,0.,1.};

       // Vector nx= normalize_vector(vector_product(nz,vz));   
       // Vector ny= normalize_vector(vector_product(nz,nx)); 
    
    // --- END of Initialization ---


    // --- Beginning of Reconstruction Process ---

    // Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start image reconstruction process ...\n");

    //Get empty event array object (type: ReadEvent)
       //sizes are equal to those of the ReconArray -> mask-width but detector-pixelsize
       //mask has to be >= detector
    int ea_size1=2*(mask->naxis1*mask->cdelt1/xpixelsize_beforeRepix);
    int ea_size2=2*(mask->naxis2*mask->cdelt2/ypixelsize_beforeRepix);
    ea=getEventArray(ea_size1,ea_size2,&status);
     
       //detector size <= mask size
       //if mask size = det size -> shift is zero
       xdiff=((ea->naxis1)/2-xdetsize_beforeRepix)/2;
       ydiff=((ea->naxis2)/2-ydetsize_beforeRepix)/2;

       // Loop over all events in the FITS file.
       while (0==EventListEOF(&eventfile->generic)) {

	 status=readEventList_nextRow(eventfile, ea);
	 CHECK_STATUS_BREAK(status);

	 //Get the 2d-EventArray, padded to the upper right corner
	 ea->EventArray[ea->rawx+ea->naxis1/2+xdiff][ea->rawy+ea->naxis2/2+ydiff]+=ea->charge;

       } // END of scanning the event list.
       CHECK_STATUS_BREAK(status);
       
       //createTestImg(&ea->EventArray,2,detector_pixels->xwidth,detector_pixels->ywidth,
       // ea->naxis1/2+xdiff,ea->naxis2/2+ydiff,"eventArray.fits",&status);
       

       if(par.RePixSize!=0.){
	 double pixelwidth_big=xpixelsize_beforeRepix;//width of the former EventArray pixels (the real det-pix-size)

	 //new pointer for EventArray with more entries, since smaller pix-size (EventArrayRepix)
	 int ear_size1=2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth);
	 int ear_size2=2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth);

	 ear=getEventArray(ear_size1, ear_size2, &status);

	 repixNoReminder(ea,ear,2,ea->naxis1,ea->naxis2,pixelwidth_big,detector_pixels->xpixelwidth);

	 ReadEvent* ea_temp=NULL;
	 ea_temp=ea;
	 FreeEventArray(ea_temp);
	 ea=ear; //set pointer to EventArray-data to just re-pixeled array

	 xdiff=((ea->naxis1)/2-detector_pixels->xwidth)/2;
	 ydiff=((ea->naxis2)/2-detector_pixels->ywidth)/2;
       }//end re-pixel EventArray to smaller size given by RePixSize

       //Get the reconstruction array:
       //type: 1: balanced cross correlation (rnd pattern), 2: MURA
       /*int type;
       if(par.protoMirax == 1){
	 type=2;
       }else{
	 type=1;
	 }*/ //TODO
       
       recon=getReconArray(mask,2,detector_pixels,&status);
       int Size1 = recon->naxis1;
       int Size2 = recon->naxis2;

        //Get the 1d image of the reconstruction array -> needed by FFTW
       ReconImage1d=SaveReconArray1d(recon, &status);  
      
       //perform a fft with the ReconArray       
       fftReconArray=FFTOfArray_1d(ReconImage1d, Size1, Size2, -1);

       //get repixeled mask from ReconArray, which is needed later for building the mask shadow during IROS
       //basic constructor for both,the whole re-pixeled mask&/shadow element
       mask_shadow=getMaskShadowElement(Size1/2, Size2/2, Size1, Size2, &status); 
       //gets re-pixeled mask as big as EventArray with values betw. 0...1
       getMaskRepix(recon, mask_shadow, 1);
  
       do{ //search for sources as long as pixval is above certain value
	 //run as long as threshold==1

	 //Get the 1d image of the event array -> needed by FFTW
	 EventImage1d=SaveEventArray1d(ea, &status);
	 //Check whether the ReconArray and the EventArray have the same size
	 if ((recon->naxis1 != ea->naxis1) || (recon->naxis2 != ea->naxis2)){
	   printf ("Error: ReconArrray and EventArray must have the same size!\n");
	   break;
	 }
   
	 //perform a fft with the EventArray       
	 fftEventArray=FFTOfArray_1d(EventImage1d, Size1, Size2, -1);
    
       //multiply fftEventArray with komplex conjugate of fftReconArray
       //Re-part: E(Re)*R(Re)+E(Im)*R(Im); Im-part: E(Re)*R(Im)-E(Im)*R(Re)       
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
       fftInvMultiply=FFTOfArray(Multiply, Size1, Size2, +1);

       //save real part of inverse fft in sky image
       for(ii=0; ii<Size1; ii++){
	 for(jj=0; jj<Size2; jj++){
	   sky_pixels->pixel[ii][jj]=fftInvMultiply[ii+Size1*jj][0]/(Size1*Size2);
	 }
       }

       //for testing:
       char name_image[MAXFILENAME];
       sprintf(name_image,"image_%lu", position_list->entryCount);

       // Write the reconstructed source function to the output FITS file.
       if(position_list->entryCount <=2){
       saveSourceImage(sky_pixels, name_image, &status);
       CHECK_STATUS_BREAK(status);
       }
   
       //finds current brightest pixel coordinates and saves PixPosition; returns current brightest pixval
       pixval=findBrightestPix(threshold, PixAmount, sky_pixels, pixval, position_list, &wcs, &status);
       threshold=getThresholdForSources(pixval, position_list, sky_pixels, median_list, par.Sigma);

       //get mask shadow for current source
       getMaskShadow2(mask_shadow,&wcs2,position_list,sky_pixels->crpix1,sky_pixels->crpix2,detector_pixels,Size1/2,Size2/2,1,&status);
       double norm=getNormalization2(mask_shadow, ea, detector_pixels, xdiff, ydiff);
       
       //new event array: method two
       if(norm>1.){

	 for(ii=0; ii<detector_pixels->xwidth; ii++){
	   for(jj=0; jj<detector_pixels->ywidth; jj++){
	     if(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]!=0.){
	       ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]-=norm*mask_shadow->shadow[ii][jj];
	       if(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]<0.){
		 ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]=0.;
	       }
	     }
	     
	   }
	 }

       }else{
	 threshold=2; 
       }
       FreeEventArray1d(EventImage1d);
       fftw_free(fftEventArray);
       fftw_free(fftInvMultiply);
        }while(threshold==1);

    //create FITS-file with all pix-coordinates
     savePositionList(position_list, par.PositionList, &status);

  // --- END of Reconstruction Process ---

  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Free the detector and sky image pixels.
  
    //set detector_pixels to original size again to be able to call destroy-fct from 'squarepixels.c'
  detector_pixels->xwidth=xdetsize_beforeRepix;
  destroySquarePixels(&detector_pixels);

  destroyCodedMask(&mask);
  FreeReconArray(&recon);
  FreeReconArray1d(ReconImage1d);
  FreeEventArray(ea);
  FreePixPositionList(position_list);
  FreeMaskShadow(mask_shadow,Size1);
  fftw_free(fftReconArray); 
  wcsfree(&wcs);
  wcsfree(&wcs2);
  free_SourceImage(sky_pixels);
  } while(0);  // END of the error handling loop.

  // Close the FITS files.
  status=closeCoMaEventFile(eventfile);
  free(eventfile);

  if (EXIT_SUCCESS==status) headas_chat(5, "finished successfully!\n\n");
  return(status);
}
//////////////////////////////////////////////////////////////////////////////////////////

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

 //Get the filename of the position list file (FITS output file).
   else if ((status=PILGetFname("PositionList", par->PositionList))) {
  SIXT_ERROR("failed reading the filename of the output position list");
  }

  //Get the right-ascension of the telescope pointing [degrees].
 else if ((status=PILGetReal("RA", &par->RA))) {
    SIXT_ERROR("failed reading the right ascension of the telescope pointing");
  }

  //Get the declination of the telescope pointing [degrees].
 else if ((status=PILGetReal("DEC", &par->DEC))) {
    SIXT_ERROR("failed reading the declination of the telescope pointing");
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

  //Read length of DCU [m].
  else if ((status=PILGetReal("DCU_length", &par->DCU_length))) {
    SIXT_ERROR("failed reading the length of DCU");
  }

  //Read length of DCU_gap [m].
  else if ((status=PILGetReal("DCU_gap", &par->DCU_gap))) {
    SIXT_ERROR("failed reading the length of DCU_gap");
  }

  //Read length of DCA_gap [m].
  else if ((status=PILGetReal("DCA_gap", &par->DCA_gap))) {
    SIXT_ERROR("failed reading the length of DCA_gap");
  }

  //Read sigma value.
  else if ((status=PILGetReal("Sigma", &par->Sigma))) {
    SIXT_ERROR("failed reading value of Sigma");
  }
  CHECK_STATUS_RET(status, status);

  return(status);
}
