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

#include "comabackpro.h"

//////////////////////////////////////////////////////////////////////////////////////////
//RECONSTRUCTION: via backprojection:                                                   //
//                                                                                      //
//                                                                                      //
//Input:event-list(t,charge,RAWX,RAWY),mask-file,detector-setup(width of det in pixels, //
//      pixelwidth,distance,pointing,gaps)                                              //
//Output:sky-image, PositionList                                                        //
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////
/** Main procedure. */
int comabackpro_main() {
  struct Parameters par;

  //struct Telescope telescope; //Telescope coordinate system
  CoMaEventFile* eventfile=NULL;
  CodedMask* mask=NULL; 
  SquarePixels* detector_pixels=NULL;
  ProjectedMask* proj_mask=NULL;
  ProjectedMask* proj_mask_repix=NULL;
  SourceImage* sky_chart=NULL;
  //SkyImage* sky_image=NULL;
  CoMaEvent* event=NULL;
  //Attitude* ac=NULL;
  PixPositionList* position_list=NULL;
  double* median_list=NULL; //temp array of all background pix for determination of median 
  //ReadEvent* ea=NULL;
  //ReconArray* recon=NULL;
  //MaskShadow* mask_shadow=NULL;

  int status=EXIT_SUCCESS; // Error status.

  int ii,jj;//kk,ll;   //counts
  int Size1, Size2;  //Sizes of ProjectedMask in pixels 
  int Size1_RePix, Size2_RePix;  //Sizes of re-pixeled ProjectedMask in pixels
  //int lastEvent=0;
  int PixAmount=10;
  // float minRA,maxRA,minDEC,maxDEC; //dimensions of skyImg in RA/DEC [deg]
  double pixelsize1, pixelsize2; //pixelsizes in meters of projectedMask
  double RePixValue; //pixelsize in meters to which the ProjectedMask is re-pixeled to
  //double att_start, att_stop; //start and stop time for current pointing from attitude-file
  //double velocity_ra,velocity_dec; //velocity of telescope motion according to current interval between att_start, att_stop
  // double timeInterval; //interval for which current pointing can be treated as constant. sub-interval for interval 
                       //between att_start, att_stop.current pointing has to be approximated accordingly
  double pixval=0.;
  int threshold=0; //for first source threshold is '0' after that: '1' if still bright enough, '2' else

  // Register HEATOOL:
  set_toolname("comabackpro");
  set_toolversion("0.01");

  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using the PIL library.
    if ((status=comabackpro_getpar(&par))) break;

    // Open the event file.
    eventfile=openCoMaEventFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Load the coded mask from the file.
    //mask->map with naxis1 x naxis2 accessed like: mask->map[x][y]
    mask=getCodedMaskFromFile(par.Mask, &status);
    CHECK_STATUS_BREAK(status);

    // DETECTOR setup.
    double ra=par.ra;
    double dec=par.dec;
    double distance=par.MaskDistance;
     
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

    //get telescope resolution in degrees
    // double res=par.Resolution/60.;
    // END of DETECTOR CONFIGURATION SETUP 

    //Set up the TELESCOPE ATTITUDE:
    //Buffer for attitude-filename.
     char att_buffer[MAXFILENAME];
    //Copy attitude-filename from par-file to buffer.
    strcpy(att_buffer, par.Attitude);
    //Make all letters capital to compare with any spelling of 'none'.
    strtoupper(att_buffer);
    /*if (0==strcmp(att_buffer, "NONE")){
      //Set up attitude(no attitude-file to use).

      //Memory-allocation:

      //Allocates memory for struct Attitude.
      //Initializes nentries, current_entry to 0;entry to NULL.
      ac=getAttitude(&status);
      CHECK_STATUS_BREAK(status);

      //Allocates memory for struct AttitudeEntry.
      ac->entry=(AttitudeEntry*)malloc(sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status=EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeEntry failed");
	break;
      }

      //Set the values of the AttitudeEntry.

      ac->nentries=1;
      ac->entry[0]=initializeAttitudeEntry();  
      //ac->entry[0].time=0;

      //Telescope pointing direction:
      ac->entry[0].nz=normalize_vector(unit_vector(ra*M_PI/180.0,dec* M_PI/180.0));
      //Unit-vector in z-direction:
      Vector vz = {0.,0.,1.};
      //Vector perpendicular to nz-vz-plane:
      ac->entry[0].nx=vector_product(vz,ac->entry[0].nz);
      //initialize telescope coordinate sytem
      telescope.nz=ac->entry[0].nz;
      telescope.nx=ac->entry[0].nx;

    }else{
      //Load the attitude from file.
      ac=loadAttitude(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      Vector initializey={0.,0.,0.};
      telescope.ny =initializey;
      telescope.nx=ac->entry[0].nx;
      }//END of setting up the TELESCOPE ATTITUDE.*/

    // PROJECTED MASK setup
    //size of axes in pixels (ProjectedMask): 2*amount_of_mask_pixels+1
    Size1=2*mask->naxis1+1;
    Size2=2*mask->naxis2+1;
    //get the different alternating pixelsizes in meters
    //NOTE: only for square pixels (of mask/det)!
    pixelsize1=mask->cdelt1-detector_pixels->xpixelwidth; //odd pixels (former mask-pix)
    pixelsize2=detector_pixels->xpixelwidth;        //even pixels ('inbetween'-pix)
    //get repix-value (smallest out of DCU_gap,DCA_gap,pixelsize1,pixelsize2 -> divided by 'factor')
    RePixValue=getRepixValue(detector_pixels,pixelsize1,/*pixelsize2,*/4);
    //size of axes in pixels (ProjectedMaskRePix): Size1(2) in meters divided by RePixSize
    Size1_RePix=(int)((mask->naxis1*pixelsize1+(mask->naxis1+1)*pixelsize2)/RePixValue)+1; //TODO:be careful if division
    //gives reminder and (int) cuts it off...here NOT the case! +1 needed for rePix-fct -> otherwise throws error!!!
    Size2_RePix=(int)((mask->naxis2*pixelsize1+(mask->naxis2+1)*pixelsize2)/RePixValue)+1;
    // END of PROJECTED MASK setup

    // SKY CHART setup.
    float delta = atan(RePixValue/distance); //TODO:probably:RePix*d/d -> ONLY RePix!
    struct SourceImageParameters sip = {
      .naxis1 = (int)((mask->naxis1*mask->cdelt1+detector_pixels->xwidth*detector_pixels->xpixelwidth)
		     /RePixValue)+2, //+1 -> (int) rounds down
      .naxis2 = (int)((mask->naxis2*mask->cdelt2+detector_pixels->ywidth*detector_pixels->ypixelwidth)
		      /RePixValue)+2, //further +1 -> for rePix-fct TODO: correct!!!
      .crpix1 = (float)((sip.naxis1-1)/2.), //-1 ->correction of above +1 TODO: correct!!!
      .crpix2 = (float)((sip.naxis2-1)/2.),
      .cdelt1 = delta, //in rad
      .cdelt2 = delta,
      .crval1 = ra*M_PI/180., //in rad
      .crval2 = dec*M_PI/180.
    };
    sky_chart=getEmptySourceImage(&sip, &status);
    CHECK_STATUS_BREAK(status);
    // END of SKY CHART CONFIGURATION SETUP 

    //initialization of WCS PARAMETER STRUCTURE for deteremining ra/dec of skyChart-pixels
       struct wcsprm wcs = {
      .flag=-1
    }; //flag has to be set only at 1st init
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis=2;
    wcs.crpix[0]=sky_chart->crpix1;
    wcs.crpix[1]=sky_chart->crpix2; 
    wcs.crval[0]=sky_chart->crval1*180./M_PI; //in deg
    wcs.crval[1]=sky_chart->crval2*180./M_PI;
    wcs.cdelt[0]=sky_chart->cdelt1*180./M_PI;
    wcs.cdelt[1]=sky_chart->cdelt2*180./M_PI;


    // SKY IMAGE setup.

    //TODO: new input parameters: possibility to set sizeOfSkyImg per hand
    //TODO: read in parameters, hide, if zero: use maxSize according to events and FOV
    //TODO: only use those events that fit into sizeOfSkyImg

    //TODO: write functions 'setSizeOfSkyImg' and 'setSizeOfSkyImgDefault'
    /*   float firstRA,lastRA,firstDEC,lastDEC;

    if (0==strcmp(att_buffer, "NONE")){//pointed observation
      minRA=ra-5.; maxRA=ra+5.;
      minDEC=dec-5.; maxDEC=dec+5.;
    }else{ //BEGIN attitude-file
      AttitudeFile* af=NULL;
      //open the attitude file:
      af=open_AttitudeFile(par.Attitude,READONLY,&status);
      CHECK_STATUS_VOID(status);
      AttitudeFileEntry afe;

      //get 1st entry
      af->row=0;
      afe=read_AttitudeFileEntry(af,&status);
      CHECK_STATUS_VOID(status);
      firstRA=afe.ra;
      firstDEC=afe.dec;

      //TODO: get last entry according to last events' time and FOV
      //(if attitude-file exceeds exposure in time -> way too huge)

      //get last entry
      af->row=ac->nentries-1;
      afe=read_AttitudeFileEntry(af,&status);
      CHECK_STATUS_VOID(status);
      lastRA=afe.ra;
      lastDEC=afe.dec;

      if (NULL!=af) {
	if (af->fptr) fits_close_file(af->fptr,&status);
	free(af);
      }

      if(firstRA > lastRA){ //in deg  //TODO: write extra 'sort-fct.'
	minRA=lastRA;
	maxRA=firstRA;
      }else{
	minRA=firstRA; 
	maxRA=lastRA;
      }

      if(firstDEC > lastDEC){ //in deg
	minDEC=lastDEC;
	maxDEC=firstDEC;
      }else{
	minDEC=firstDEC;
	maxDEC=lastDEC;
      }

    }//END Attitude-file
    //get max. FOV
    float fov = det_phi_max(distance,mask->naxis1*mask->cdelt1,mask->naxis2*mask->cdelt2,
			    detector_pixels->xwidth*detector_pixels->xpixelwidth, 
			    detector_pixels->ywidth*detector_pixels->ypixelwidth);
    fov=fov*180./M_PI; //in deg
    minRA-=fov/2.; //RA: 0...360
    maxRA+=fov/2.;
    if(minRA < 0.){minRA=360.+ minRA;}
    if(maxRA > 360.){maxRA=maxRA-360.;}
    firstRA=minRA;lastRA=maxRA;
    if(firstRA > lastRA){ //in deg
      minRA=lastRA;
      maxRA=firstRA;
    }else{
      minRA=firstRA; 
      maxRA=lastRA;
    }
    minDEC-=fov/2.; //DEC: -180...180
    maxDEC+=fov/2.;
    if(minDEC < -180.){minDEC=360.+ minDEC;}
    if(maxDEC > 180.){maxDEC=maxDEC-360.;}
    firstDEC=minDEC;lastDEC=maxDEC;
    if(firstDEC > lastDEC){ //in deg
      minDEC=lastDEC;
      maxDEC=firstDEC;
    }else{
      minDEC=firstDEC;
      maxDEC=lastDEC;
    }

    struct SkyImageParameters skp = {
      .naxis1 = (int)(fabs(maxRA-minRA)/res),
      .naxis2 = (int)(fabs(maxDEC-minDEC)/res),
      .crpix1 = (float)(skp.naxis1/2.+1), 
      .crpix2 = (float)(skp.naxis2/2.+1),
      .cdelt1 = res, //in deg
      .cdelt2 = res,
      .crval1 = minRA+skp.crpix1*res, //in deg
      .crval2 = minDEC+skp.crpix2*res,
      .minra  = minRA,  //in deg
      .maxra  = maxRA,
      .mindec = minDEC, //in deg 
      .maxdec = maxDEC
    };
    sky_image=getEmptySkyImage(&skp, &status);
    CHECK_STATUS_BREAK(status);

    //fill in ra- and dec-array for skyImage
    fillRaDecArrays(sky_image);
    
    // END of SKY IMAGE CONFIGURATION SETUP

    */
    //SOURCE-POSITION DETERMINATION:
    //get empty PixPositionList structure (contains pointer to current PixPosition-element
    //and count for found sources)
    position_list=getPixPositionList(sky_chart);
    //memory-allocation for median_list
    median_list=getMedian_list(sky_chart, &status);


    //initialization of wcs parameter structure for getting mask shadow
    /*   struct wcsprm wcs2 = {
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
    wcs2.cdelt[1]=atan(detector_pixels->ypixelwidth/distance)*180./M_PI;*/


    // --- END of Initialization ---


    // --- Beginning of Reconstruction Process ---

    // Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start image reconstruction process ...\n");
    
    //empty ProjectedMask object
    proj_mask=getEmptyProjectedMask(Size1,Size2,pixelsize1,pixelsize2,&status);
    //filled ProjectedMask object derived from mask pattern
    getProjectedMask(mask,proj_mask);

    getOpenPixels(proj_mask); //amount of open pixels in PM
    //empty ProjectedMask object for re-pixeled ProjectedMask
    proj_mask_repix=getEmptyProjectedMask(Size1_RePix,Size2_RePix,RePixValue,RePixValue,&status);
    repixWithReminder(proj_mask,proj_mask_repix,4,Size1,Size2,pixelsize2,pixelsize1,RePixValue,0.);
    getOpenPixels(proj_mask_repix); //amount of open pixels in re-pixeled PM



    //Get the reconstruction array:
    /*    recon=getReconArray(mask, detector_pixels, &status);
    int Size1 = recon->naxis1;
    int Size2 = recon->naxis2;
    //get repixeled mask from ReconArray, which is needed later for building the mask shadow during IROS
    //basic constructor for both,the whole re-pixeled mask&/shadow element
    mask_shadow=getMaskShadowElement(Size1/2,Size2/2,detector_pixels->xwidth,detector_pixels->ywidth,&status); 
    //gets re-pixeled mask as big as EventArray with values betw. 0...1
    getMaskRepix(recon,mask_shadow,0);*/
    
    
    //get CoMaEvent; pointer for current event in eventfile
    event=getCoMaEvent(&status);
    //get 1st event
    status=CoMaEventFile_getNextRow(eventfile,event);
    CHECK_STATUS_BREAK(status);
   

    /*    ea=getEventArray(detector_pixels->xwidth,detector_pixels->ywidth,&status);
    // Loop over all events in the FITS file.
    while (0==EventListEOF(&eventfile->generic)) {

      status=readEventList_nextRow(eventfile, ea);
      CHECK_STATUS_BREAK(status);

      //Get the 2d-EventArray, padded to the upper right corner
      ea->EventArray[ea->rawx][ea->rawy]+=ea->charge;

    } // END of scanning the event list.
    CHECK_STATUS_BREAK(status);*/

    if (0==strcmp(att_buffer, "NONE")){//pointed observation

      /*   do{ //search for sources as long as pixval is above certain value
      //run as long as threshold==1*/

      float count=1.;
      while (0==EventListEOF(&eventfile->generic)) {
	// printf("%lf\n%",count/eventfile->generic.nrows*100.);
	/*	for(kk=0; kk<ea->naxis1; kk++){
		for(ll=0; ll<ea->naxis2; ll++){*/

	    //add projected mask(pm) correctly to SkyChart for current photon/event
	    //since pm-/SkyChart-pixels fit without reminder into detector-pixels:
	    //detpix=(pixel of event)*(amount of smaller pixels within one detector pixels)
	    //+1 -> (int) rounds down;
	    int detpix_x=(event->rawx)*(int)(detector_pixels->xpixelwidth/RePixValue+1);
	    int detpix_y=(event->rawy)*(int)(detector_pixels->ypixelwidth/RePixValue+1);
	  
	    //scanning all projected mask pixels to fill in
	    //skyChart (still pix-values)
	 
	    for(ii=0; ii<proj_mask_repix->naxis1; ii++){ 
	      for(jj=0; jj<proj_mask_repix->naxis2; jj++){
	    
		int shift_ii=proj_mask_repix->naxis1-1-ii;
		int shift_jj=proj_mask_repix->naxis2-1-jj;

		sky_chart->pixel[shift_ii+detpix_x][shift_jj+detpix_y]+=event->charge/proj_mask_repix->OpenPixels
		  *proj_mask_repix->map[ii][jj];
	      }
	    }

	    status=CoMaEventFile_getNextRow(eventfile, event);
	    count++;
      }

	    /*	  }
		  }*/


	//for testing:
	char name_image[MAXFILENAME];
	//sprintf(name_image,"skyImg_%lf.fits", position_list->entryCount);

	// Write the reconstructed source function to the output FITS file.
	//  if(position_list->entryCount <=4){
	  saveSourceImage(sky_chart, name_image, &status);
	  CHECK_STATUS_BREAK(status);
	  //	}
   
	do{ //search for sources as long as pixval is above certain value
	//run as long as threshold==1

	//finds current brightest pixel coordinates and saves PixPosition; returns current brightest pixval
	pixval=findBrightestPix(threshold, PixAmount,sky_chart, pixval, position_list, &wcs, &status);
	threshold=getThresholdForSources(pixval, position_list, sky_chart, median_list, par.Sigma);
	//printf("%f\n",pixval);
	//get mask shadow for current source
	//mask in det-pixel-size
	/*	float sky_crpix1=mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth+1;
	float sky_crpix2=mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth+1;
	getMaskShadow2(mask_shadow,&wcs2,position_list,sky_crpix1,sky_crpix2,detector_pixels,Size1/2,Size2/2,0,&status);
	createTestImg(mask_shadow,1,detector_pixels->xwidth,detector_pixels->ywidth,0,0,"maskshadow.fits",&status);
	double norm=getNormalization2(mask_shadow, ea, detector_pixels, (-1)*ea->naxis1/2, (-1)*ea->naxis2/2);
	//new event array: method two
	if(norm>1.){

	  for(ii=0; ii<detector_pixels->xwidth; ii++){
	    for(jj=0; jj<detector_pixels->ywidth; jj++){
	      if(ea->EventArray[ii][jj]!=0.){
		ea->EventArray[ii][jj]-=norm*mask_shadow->shadow[ii][jj];
		if(ea->EventArray[ii][jj]<0.){
		  ea->EventArray[ii][jj]=0.;
		}
	      }
	     
	    }
	  }

	  int xcount,ycount;
	  //reset SkyChart to zero for next round
	  for(xcount=0; xcount<sky_chart->naxis1; xcount++){
	    for(ycount=0; ycount<sky_chart->naxis2; ycount++){
	      sky_chart->pixel[xcount][ycount] = 0.;
	    }
	    }



	}else{
	  threshold=2; 
	  }*/
      }while(threshold==1);
     
      //create FITS-file with all pix-coordinates
      savePositionList(position_list, "posList.fits", &status);


 

    }else{//BEGIN attitude-file
      //get 1st telescope pointing
    /*     telescope.nz=getTelescopeNz(ac,event->time,&status);
      CHECK_STATUS_BREAK(status);
      //update sky_chart and wcs for current pointing
      setWCScurrentPointing(par.Attitude,ac,&telescope.nz,&wcs,&status); //in deg
      sky_chart->crval1= wcs.crval[0]*M_PI/180.; //in rad
      sky_chart->crval2= wcs.crval[1]*M_PI/180.;

      //set att_start.att_stop,velocity according to 1st pointing
      att_start=ac->entry[ac->currentry].time;
      att_stop=ac->entry[ac->currentry+1].time;
      getCurrentVelocity(par.Attitude,ac,&velocity_ra,&velocity_dec,att_start,att_stop,&status);

      //set 1st timeInterval according to velocity and resolution
      double time_ra=res/velocity_ra; 
      double time_dec=res/velocity_dec;
      if(time_ra < time_dec){
	timeInterval=time_ra;
      }else{
	timeInterval=time_dec;
      }
    
      // Loop over all events in the FITS file.
      while (0==EventListEOF(&eventfile->generic)) {
      
	while(event->time < att_stop){ //current pointing in attitude-file
	
	  while(event->time < (att_start+timeInterval)){//current 'constant'-interval

	    //add projected mask(pm) correctly to SkyChart for current photon/event
	    //since pm-/SkyChart-pixels fit without reminder into detector-pixels:
	    //detpix=(pixel of event)*(amount of smaller pixels within one detector pixels)
	    //+1 -> (int) rounds down;
	    int detpix_x=(detector_pixels->xwidth-event->rawx-1)*(int)(detector_pixels->xpixelwidth/RePixValue+1);
	    int detpix_y=(detector_pixels->ywidth-event->rawy-1)*(int)(detector_pixels->ypixelwidth/RePixValue+1);
	  
	    //scanning all projected mask pixels to fill in
	    //skyChart (still pix-values)
	 
	    for(ii=0; ii<proj_mask_repix->naxis1; ii++){ 
	      for(jj=0; jj<proj_mask_repix->naxis2; jj++){
          
		sky_chart->pixel[ii+detpix_x][jj+detpix_y]+=event->charge/proj_mask_repix->OpenPixels
		  *proj_mask_repix->map[ii][jj];
	      }
	    }

	    //get all events for current constant interval
	    status=CoMaEventFile_getNextRow(eventfile, event);
	  
	    if(status == 1){lastEvent=1;status=0;}

	    //saveSourceImage(sky_chart,"skyChart.fits", &status);
	  
	    if(lastEvent == 1){break;}
	  }//END of current 'constant'-interval

	  saveSourceImage(sky_chart,"skyChart.fits", &status);

	  //TODO: fill in skyChart at corresponding position in skyImg
	  for(ii=0; ii<sky_chart->naxis1; ii++){ 
	    for(jj=0; jj<sky_chart->naxis2; jj++){
	      double pixcrd[2]={ii+1,jj+1};
	      double imgcrd[2];
	      double world[2];
	      double phi, theta;
	      int stat=0;
	      //per pix: p2s -> RA/DEC value for current skyChart-pix
	      wcsp2s(&wcs,1,2,pixcrd,imgcrd,&phi,&theta,world,&stat);
	      if(0!=stat){
		status=EXIT_FAILURE;
		HD_ERROR_THROW("wcs projection failed!\n", status);
	      }
	      double pix_ra=world[0];
	      double pix_dec=world[1];
	      
	      //(RA-minra)/cdelt -> index for corresponding skyImg-pix
	      int sky_ii=(int)((pix_ra-sky_image->minra)/sky_image->cdelt1);
	      int sky_jj=(int)((pix_dec-sky_image->mindec)/sky_image->cdelt2);
	      //add value of current skyChart-pix to newly found pix in skyImg
	      sky_image->pixel[sky_ii][sky_jj]+=sky_chart->pixel[ii][jj];
	    }
	  }
	
	  int xcount,ycount;
	  //reset SkyChart to zero for next round
	  for(xcount=0; xcount<sky_chart->naxis1; xcount++){
	    for(ycount=0; ycount<sky_chart->naxis2; ycount++){
	      sky_chart->pixel[xcount][ycount] = 0.;
	    }
	  }

	  saveSkyImage(sky_image,"SkyImg.fits",&status);

	  //get pointing inbetween (approximation), increase 'timeInterval'
	  telescope.nz=getTelescopeNz(ac,event->time,&status);
	  CHECK_STATUS_BREAK(status);
	  timeInterval+=timeInterval;

	  //update sky_chart and wcs for current pointing
	  setWCScurrentPointing(par.Attitude,ac,&telescope.nz,&wcs,&status); //in deg
	  sky_chart->crval1= wcs.crval[0]*M_PI/180.; //in rad
	  sky_chart->crval2= wcs.crval[1]*M_PI/180.;
	
	  if(lastEvent == 1){break;}
	}//END of current pointing in attitude-file

	//get new pointing (att_start,att_stop,velocity), new constant 'timeInterval'
	ac->currentry++;
	att_start=ac->entry[ac->currentry].time;
	att_stop=ac->entry[ac->currentry+1].time;
	getCurrentVelocity(par.Attitude,ac,&velocity_ra,&velocity_dec,att_start,att_stop,&status);
	double time_ra=res/velocity_ra; 
	double time_dec=res/velocity_dec;
	if(time_ra < time_dec){
	  timeInterval=time_ra;
	}else{
	  timeInterval=time_dec;
	}

	//get pointing inbetween (approximation), increase 'timeInterval'
	telescope.nz=getTelescopeNz(ac,event->time,&status);
	CHECK_STATUS_BREAK(status);
	//update sky_chart and wcs for current pointing
	setWCScurrentPointing(par.Attitude,ac,&telescope.nz,&wcs,&status); //in deg
	sky_chart->crval1= wcs.crval[0]*M_PI/180.; //in rad
	sky_chart->crval2= wcs.crval[1]*M_PI/180.;

	} // END of scanning the event list.*/
      }//END attitude-file
    
   } while(0);  // END of the error handling loop.

  //saveSkyImage(sky_image,"SkyImg_final.fits",&status); 




    // --- END of Reconstruction Process ---

    // --- Cleaning up ---
    headas_chat(5, "cleaning up ...\n");

   // Free the detector and sky image pixels.
   destroySquarePixels(&detector_pixels);
   destroyCodedMask(&mask);



  // Close the FITS files.
  status=closeCoMaEventFile(eventfile);
  free(eventfile);

  if (EXIT_SUCCESS==status) headas_chat(5, "finished successfully!\n\n");
  return(status);
}
//////////////////////////////////////////////////////////////////////////////////////////

int comabackpro_getpar(struct Parameters* par)
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

  // Get the filename of the attitude file (FITS input file).
  else if ((status=PILGetFname("Attitude", par->Attitude))) {
    SIXT_ERROR("failed reading the filename of the attitude file");
  }

  //Get the filename of the image file (FITS output file).
   else if ((status=PILGetFname("Image", par->Image))) {
  SIXT_ERROR("failed reading the filename of the output image");
  }

 //Get the filename of the position list file (FITS output file).
   else if ((status=PILGetFname("PositionList", par->PositionList))) {
  SIXT_ERROR("failed reading the filename of the output position list");
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


    // Read the width of one detector pixel in [m].
  else if ((status=PILGetReal("RA", &par->ra))) {
    SIXT_ERROR("failed reading the width of the detector pixels");
  }

    // Read the width of one detector pixel in [m].
  else if ((status=PILGetReal("DEC", &par->dec))) {
    SIXT_ERROR("failed reading the width of the detector pixels");
  }

  // Read the resolution of the telescope [arcmin].
  else if ((status=PILGetReal("Resolution", &par->Resolution))) {
    SIXT_ERROR("failed reading the resolution of the telescope");
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
