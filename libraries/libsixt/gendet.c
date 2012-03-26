#include "gendet.h"


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


GenDet* newGenDet(const char* const filename, int* const status) 
{
  // Allocate memory.
  GenDet* det=(GenDet*)malloc(sizeof(GenDet));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }

  // Initialize all pointers with NULL.
  det->telescope=NULL;
  det->pixgrid=NULL;
  det->split  =NULL;
  det->line=NULL;
  det->psf =NULL;
  det->vignetting=NULL;
  det->coded_mask=NULL;
  det->arf =NULL;
  det->rmf =NULL;
  det->clocklist =NULL;
  det->badpixmap =NULL;
  det->filepath  =NULL;
  det->filename  =NULL;

  // Set initial values.
  det->erobackground = 0;
  det->frametime     = 0.;
  det->anyphoton     = 0;

  // Get empty GenPixGrid.
  det->pixgrid = newGenPixGrid(status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Get empty ClockList.
  det->clocklist = newClockList(status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Get empty split model.
  det->split = newGenSplit(status);
  if (EXIT_SUCCESS!=*status) return(det);

  // Split the reference to the XML detector definition file
  // into path and filename. This has to be done before
  // calling the parser routine for the XML file.
  char filename2[MAXFILENAME];
  char rootname[MAXFILENAME];
  // Make a local copy of the filename variable in order to avoid
  // compiler warnings due to discarded const qualifier at the 
  // subsequent function call.
  strcpy(filename2, filename);
  fits_parse_rootname(filename2, rootname, status);
  CHECK_STATUS_RET(*status, det);

  // Split rootname into the file path and the file name.
  char* lastslash = strrchr(rootname, '/');
  if (NULL==lastslash) {
    det->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(det->filepath, *status, 
		   "memory allocation for filepath failed", det);
    det->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(det->filename, *status, 
		   "memory allocation for filename failed", det);
    strcpy(det->filepath, "");
    strcpy(det->filename, rootname);
  } else {
    lastslash++;
    det->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(det->filename, *status, 
		   "memory allocation for filename failed", det);
    strcpy(det->filename, lastslash);
      
    *lastslash='\0';
    det->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(det->filepath, *status, 
		   "memory allocation for filepath failed", det);
    strcpy(det->filepath, rootname);
  }
  // END of storing the filename and filepath.


  // Read in the XML definition of the detector.
  parseGenDetXML(det, filename, status);
  if (EXIT_SUCCESS!=*status) return(det);
    
  // Allocate memory for the pixels.
  det->line=(GenDetLine**)malloc(det->pixgrid->ywidth*sizeof(GenDetLine*));
  if (NULL==det->line) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }
  int ii;
  for (ii=0; ii<det->pixgrid->ywidth; ii++) {
    det->line[ii] = newGenDetLine(det->pixgrid->xwidth, status);
    if (EXIT_SUCCESS!=*status) return(det);
  }

  return(det);
}


void destroyGenDet(GenDet** const det, int* const status)
{
  if (NULL!=*det) {
    if (NULL!=(*det)->telescope) {
      free((*det)->telescope);
    }
    if (NULL!=(*det)->line) {
      int i;
      for (i=0; i<(*det)->pixgrid->ywidth; i++) {
	destroyGenDetLine(&(*det)->line[i]);
      }
      free((*det)->line);
    }

    if (NULL!=(*det)->filepath) {
      free((*det)->filepath);
    }

    if (NULL!=(*det)->filename) {
      free((*det)->filename);
    }

    // Destroy the ClockList.
    destroyClockList(&(*det)->clocklist);

    // Destroy the GenPixGrid.
    destroyGenPixGrid(&(*det)->pixgrid);

    // Destroy the split model.
    destroyGenSplit(&(*det)->split);

    // Destroy the BadPixMap.
    destroyBadPixMap(&(*det)->badpixmap);

    // Destroy the PSF.
    destroyPSF(&(*det)->psf);

    // Destroy the CodedMask.
    destroyCodedMask(&(*det)->coded_mask);

    // Destroy the vignetting Function.
    destroyVignetting(&(*det)->vignetting);

    // Free the cosmic ray background model.
    if (1==(*det)->erobackground) {
      eroBkgCleanUp(status);
    }

    free(*det);
    *det=NULL;
  }
}


int addGenDetPhotonImpact(GenDet* const det, 
			  const Impact* const impact, 
			  EventListFile* const elf,
			  int* const status)
{
  // Call the detector operating clock routine.
  operateGenDetClock(det, elf, impact->time, status);
  if (EXIT_SUCCESS!=*status) return(0);

  // Determine the detected energy.
  float energy;
  
  if (NULL!=det->rmf) {
    // Determine the measured detector channel (PHA channel) according 
    // to the RMF.
    // The channel is obtained from the RMF using the corresponding
    // HEAdas routine which is based on drawing a random number.
    long channel;
    returnRMFChannel(det->rmf, impact->energy, &channel);

    // Check if the photon is really measured. If the
    // PHA channel returned by the HEAdas RMF function is '-1', 
    // the photon is not detected.
    // This can happen, if the RMF actually is an RSP, i.e. it 
    // includes ARF contributions, e.g., 
    // the detector quantum efficiency and filter transmission.
    if (channel<det->rmf->FirstChannel) {
      return(0); // Break the function (photon is not detected).
    }

    // Determine the corresponding detected energy.
    // NOTE: In this simulation the collected charge is represented 
    // by the nominal photon energy [keV], which corresponds to the 
    // PHA channel according to the EBOUNDS table.
    energy=getEBOUNDSEnergy(channel, det->rmf, 0);
    assert(energy>=0.);

  } else {
    // The detector has no particular RMF. Therefore we directly
    // use the energy of the incident photon.
    energy=impact->energy;
  }

  // Create split events.
  int npixels=makeGenSplitEvents(det, &impact->position, energy, 
				 impact->ph_id, impact->src_id, 
				 impact->time, elf, status);
  CHECK_STATUS_RET(*status, npixels);

  // Set the flag that there has been a photon interaction.
  det->anyphoton=1;

  // Return the number of affected pixels.
  return(npixels);
}


void operateGenDetClock(GenDet* const det, EventListFile* const elf,
			const double time, int* const status)
{
  // Check if the detector operation setup is time-triggered.
  if (GENDET_TIME_TRIGGERED!=det->readout_trigger) return;

  // Get the next element from the clock list.
  CLType type;
  void* element=NULL;
  do {
    CLReadoutLine* clreadoutline=NULL;
    CLClearLine*   clclearline  =NULL;
    CLWait* clwait              =NULL;

    getClockListElement(det->clocklist, time, &type, &element, status);
    if (EXIT_SUCCESS!=*status) return;

    switch (type) {
    case CL_NONE:
      // No operation has to be performed. The clock list is
      // currently in a wait status.
      break;
    case CL_NEWFRAME:
      // The clock list has internally increased the frame counter and readout 
      // time. 

      // If there has been no photon interaction during the last frame
      // and if no background model is activated, jump over the next empty frames
      // until there is a new photon impact.
      if ((0==det->erobackground)&&(0==det->anyphoton)) {
	long nframes=(long)((time-det->clocklist->readout_time)/det->frametime);
	det->clocklist->time       +=nframes*det->frametime;
	det->clocklist->frame      +=nframes;
	det->clocklist->readout_time=det->clocklist->time;
      }

      // Reset the flag.
      det->anyphoton=0;

      break;
    case CL_WAIT:
      // A waiting period is finished.
      clwait = (CLWait*)element;

      // Insert cosmic ray background events, if the appropriate model is defined.
      if (1==det->erobackground) {
	// Get background events for the required time interval (has
	// to be given in [s]).
	eroBackgroundOutput* list=eroBkgGetBackgroundList(clwait->time);
	int ii;
	for(ii = 0; ii<list->numhits; ii++) {
	  // Position of the event.
	  struct Point2d position = {
	    .x = list->hit_xpos[ii] *0.001,
	    .y = list->hit_ypos[ii] *0.001
	  };

	  makeGenSplitEvents(det, &position,
			     // Energy of the event in [keV].
			     list->hit_energy[ii],
			     -1, -1, time, elf, status);
	  CHECK_STATUS_BREAK(*status);
	} 
	eroBkgFree(list);
      }

      // Apply the bad pixel map (if available) with the bad pixel values weighted 
      // by the waiting time.
      if (NULL!=det->badpixmap) {
	applyBadPixMap(det->badpixmap, clwait->time, encounterGenDetBadPix, det->line);
      }
      break;
    case CL_LINESHIFT:
      GenDetLineShift(det);
      break;
    case CL_READOUTLINE:
      clreadoutline=(CLReadoutLine*)element;
      GenDetReadoutLine(det, clreadoutline->lineindex, 
			clreadoutline->readoutindex, 
			elf,
			status);
      break;
    case CL_CLEARLINE:
      clclearline = (CLClearLine*)element;
      GenDetClearLine(det, clclearline->lineindex);
      break;
    }
    CHECK_STATUS_VOID(*status);
  } while(type!=CL_NONE);
}


void GenDetLineShift(GenDet* const det)
{
  headas_chat(5, "lineshift\n");

  // Check if the detector contains more than 1 line.
  if (2>det->pixgrid->ywidth) return;

  // Apply the Charge Transfer Efficiency.
  int ii;
  if (det->cte!=1.) { 
    for (ii=1; ii<det->pixgrid->ywidth; ii++) {
      if (0!=det->line[ii]->anycharge) {
	int jj;
	for (jj=0; jj<det->line[ii]->xwidth; jj++) {
	  if (det->line[ii]->charge[jj] > 0.) {
	    det->line[ii]->charge[jj] *= det->cte;
	  }
	}
      }
    }
  }

  // Add the charges in line 1 to line 0.
  addGenDetLine(det->line[0], det->line[1]);

  // Clear the charges in line 1, as they are now contained in line 0.
  clearGenDetLine(det->line[1]);

  // Shift the other lines in increasing order and put the newly cleared 
  // original line number 1 at the end as the last line.
  GenDetLine* buffer = det->line[1];
  for (ii=1; ii<det->pixgrid->ywidth-1; ii++) {
    det->line[ii] = det->line[ii+1];
  }
  det->line[det->pixgrid->ywidth-1] = buffer;
}


static inline void GenDetReadoutPixel(GenDet* const det, 
				      const int lineindex, 
				      const int readoutindex,
				      const int xindex,
				      const double time,
				      EventListFile* const elf,
				      int* const status)
{
  headas_chat(5, "read out line %d as %d\n", lineindex, readoutindex);

  GenDetLine* line = det->line[lineindex];

  if (line->charge[xindex]>0.) {

    Event* event=NULL;

    // Error handling loop.
    do {

      // Determine the properties of a new Event object.
      event=getEvent(status);
      CHECK_STATUS_BREAK(*status);
    
      // Readout the signal from the pixel array ...
      event->signal = line->charge[xindex];
      // ... and delete the pixel value.
      line->charge[xindex] = 0.;
    
      // Copy the information about the original photons.
      int jj;
      for(jj=0; jj<NEVENTPHOTONS; jj++) {
	event->ph_id[jj]  = line->ph_id[xindex][jj];
	event->src_id[jj] = line->src_id[xindex][jj];
	line->ph_id[xindex][jj]  = 0;
	line->src_id[xindex][jj] = 0;
      }
      
      // Apply the charge thresholds.
      if (event->signal<=det->threshold_readout_lo_keV) {
	break;
      }
      if (det->threshold_readout_up_keV >= 0.) {
	if (event->signal>=det->threshold_readout_up_keV) {
	  break;
	}
      }
    
      // Apply the detector response if available.
      if (NULL!=det->rmf) {
	event->pha=getEBOUNDSChannel(event->signal, det->rmf);
      } else {
	event->pha=0;
      }

      // Store remaining information.
      event->rawy =readoutindex;
      event->rawx =xindex;
      event->time =time;  // Time of detection.
      event->frame=det->clocklist->frame; // Frame of detection.

      // Store the event in the output event file.
      addEvent2File(elf, event, status);
      CHECK_STATUS_BREAK(*status);

    } while(0); // END of error handling loop.

    // Release memory
    freeEvent(&event);
  }
}


void GenDetReadoutLine(GenDet* const det, 
		       const int lineindex, 
		       const int readoutindex, 
		       EventListFile* const elf,
		       int* const status)
{
  headas_chat(5, "read out line %d as %d\n", lineindex, readoutindex);

  GenDetLine* line = det->line[lineindex];

  if (0!=line->anycharge) {
    int ii;
    for (ii=0; ii<line->xwidth; ii++) {

      GenDetReadoutPixel(det, lineindex, readoutindex, ii,
			 det->clocklist->readout_time, elf, status);
      CHECK_STATUS_BREAK(*status);
    }
    CHECK_STATUS_VOID(*status);
    // END of loop over all pixels in the line.

    // Reset the anycharge flag of this line.
    line->anycharge=0;
  }
}


void GenDetClearLine(GenDet* const det, const int lineindex) {
  clearGenDetLine(det->line[lineindex]);
}


void encounterGenDetBadPix(void* const data, 
			   const int x, const int y, 
			   const float value) 
{
  // Array of pointers to pixel lines.
  GenDetLine** line = (GenDetLine**)data;

  // Check if the bad pixel type.
  if (value < 0.) { // The pixel is a cold one.
    // Delete the charge in the pixel.
    line[y]->charge[x] = 0.;
  } else if (value > 0.) { // The pixel is a hot one.
    // Add additional charge to the pixel.
    addGenDetCharge2Pixel(line[y], x, value, -1, -1);
  }
}


GenSplit* newGenSplit(int* const status) 
{
  // Allocate memory.
  GenSplit* split=(GenSplit*)malloc(sizeof(GenSplit));
  if (NULL==split) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenSplit failed!\n", 
		   *status);
    return(split);
  }

  // Initialize all pointers with NULL.

  // Set default values.
  split->type=GS_NONE;
  split->par1=0.;
  split->par2=0.;

  return(split);
}


void destroyGenSplit(GenSplit** const split)
{
  if (NULL!=*split) {
    free(*split);
    *split=NULL;
  }
}


static inline int getMinimumDistance(const double array[]) 
{
  int count, index=0;
  double minimum=array[0];

  for (count=1; count<4; count++) {
    if ( (minimum < 0.) ||
	 ((array[count]<=minimum)&&(array[count]>=0.)) ) {
      minimum = array[count];
      index = count;
    }
  }

  return(index);
}


int makeGenSplitEvents(GenDet* const det,
		       const struct Point2d* const position,
		       const float signal,
		       const long ph_id, const long src_id,
		       const double time,
		       EventListFile* const elf, 
		       int* const status)
{
  // Number of affected pixels.
  int npixels=0;
  // x- and y-indices of affected pixels.
  int x[4], y[4];
  // Signal fractions in the individual pixels.
  float fraction[4];

  // The following array entries are used to transform between 
  // different array indices for accessing neighboring pixels.
  const int xe[4] = {1, 0,-1, 0};
  const int ye[4] = {0, 1, 0,-1};

  // Which kind of split model has been selected?
  if (GS_NONE==det->split->type) {
    // No split events => all events are singles.
    npixels=1;

    // Determine the affected detector line and column.
    x[0]=getGenDetAffectedColumn(det->pixgrid, position->x);
    y[0]=getGenDetAffectedLine  (det->pixgrid, position->y);

    // Check if the returned values are valid line and column indices.
    if ((x[0]<0) || (y[0]<0)) {
      return(0);
    }
    
    // The single pixel receives the total photon energy.
    fraction[0]=1.;
    
  } else if (GS_GAUSS==det->split->type) {  
    // Gaussian split model.

    // Signal cloud sigma as a function of the photon energy.
    const float ccsigma= 
      det->split->par1 + det->split->par2*sqrt(signal);

    // Signal cloud size (3 sigma).
    const float ccsize=ccsigma*3.;

    // Calculate pixel indices (integer) of the central affected pixel:
    x[0]=getGenDetAffectedColumn(det->pixgrid, position->x);
    y[0]=getGenDetAffectedLine  (det->pixgrid, position->y);
  
    // Check if the impact position lies inside the detector pixel array.
    if ((0>x[0]) || (0>y[0])) {
      return(0);
    }

    // Calculate the distances from the impact center position to the 
    // borders of the surrounding pixel (in [m]):
    double distances[4] = {
      // Distance to right pixel edge.
      (x[0]-det->pixgrid->xrpix+1.5)*det->pixgrid->xdelt + 
      det->pixgrid->xrval - position->x,
      // Distance to upper edge.
      (y[0]-det->pixgrid->yrpix+1.5)*det->pixgrid->ydelt + 
      det->pixgrid->yrval - position->y,
      // Distance to left pixel edge.
      position->x - ((x[0]-det->pixgrid->xrpix+0.5)*det->pixgrid->xdelt + 
		     det->pixgrid->xrval),
      // distance to lower edge
      position->y - ((y[0]-det->pixgrid->yrpix+0.5)*det->pixgrid->ydelt + 
		     det->pixgrid->yrval)
    };

    int mindist = getMinimumDistance(distances);
    if (distances[mindist] < ccsize) {
      // Not a single event!
      x[1] = x[0] + xe[mindist];
      y[1] = y[0] + ye[mindist];

      double mindistgauss = gaussint(distances[mindist]/ccsigma);

      // Search for the next to minimum distance to an edge.
      double minimum = distances[mindist];
      distances[mindist] = -1.;
      int secmindist = getMinimumDistance(distances);
      distances[mindist] = minimum;

      if (distances[secmindist] < ccsize) {
	// Quadruple!
	npixels = 4;

	x[2] = x[0] + xe[secmindist];
	y[2] = y[0] + ye[secmindist];
	x[3] = x[1] + xe[secmindist];
	y[3] = y[1] + ye[secmindist];

	// Calculate the different signal fractions in the 4 affected pixels.
	double secmindistgauss = gaussint(distances[secmindist]/ccsigma);
	fraction[0] = (1.-mindistgauss)*(1.-secmindistgauss);
	fraction[1] =     mindistgauss *(1.-secmindistgauss);
	fraction[2] = (1.-mindistgauss)*    secmindistgauss ;
	fraction[3] =     mindistgauss *    secmindistgauss ;

      } else {
	// Double!
	npixels = 2;

	fraction[0] = 1. - mindistgauss;
	fraction[1] =      mindistgauss;

      } // END of double or Quadruple

    } else {
      // Single event!
      npixels = 1;
      fraction[0] = 1.;
      
    } 
    // END of check for single event

    // END of Gaussian split model.

  } else if (GS_EXPONENTIAL==det->split->type) {  
    // Exponential split model.
    // None-Gaussian, exponential signal cloud model 
    // (concept proposed by Konrad Dennerl).
    npixels=4;

    // Calculate pixel indices (integer) of central affected pixel:
    x[0] = getGenDetAffectedColumn(det->pixgrid, position->x);
    y[0] = getGenDetAffectedLine  (det->pixgrid, position->y);
  
    // Check if the impact position lies inside the detector pixel array.
    if ((0>x[0]) || (0>y[0])) {
      return(0);
    }

    // Calculate the distances from the impact center position to the 
    // borders of the surrounding pixel (in units [fraction of a pixel edge]):
    double distances[4] = {
      // Distance to right pixel edge.
      x[0]-det->pixgrid->xrpix+1.5 + 
      (det->pixgrid->xrval - position->x)/det->pixgrid->xdelt,
      // Distance to upper edge.
      y[0]-det->pixgrid->yrpix+1.5 + 
      (det->pixgrid->yrval - position->y)/det->pixgrid->ydelt,
      // Distance to left pixel edge.
      (position->x - det->pixgrid->xrval)/det->pixgrid->xdelt - 
      (x[0]-det->pixgrid->xrpix+0.5),
      // distance to lower edge
      (position->y - det->pixgrid->yrval)/det->pixgrid->ydelt - 
      (y[0]-det->pixgrid->yrpix+0.5)
    };

    // Search for the minimum distance to the edges.
    int mindist = getMinimumDistance(distances);
    x[1] = x[0] + xe[mindist];
    y[1] = y[0] + ye[mindist];

    // Search for the next to minimum distance to the edges.
    double minimum = distances[mindist];
    distances[mindist] = -1.;
    int secmindist = getMinimumDistance(distances);
    distances[mindist] = minimum;
    // Pixel coordinates of the 3rd and 4th split partner.
    x[2] = x[0] + xe[secmindist];
    y[2] = y[0] + ye[secmindist];
    x[3] = x[1] + xe[secmindist];
    y[3] = y[1] + ye[secmindist];

    // Now we know the affected pixels and can determine the 
    // signal fractions according to the model exp(-(r/0.355)^2).
    // Remember that the array distances[] contains the distances
    // to the pixel borders, whereas here we need the distances from
    // the pixel center for the parameter r.
    // The value 0.355 is given by the parameter ecc->parameter.
    fraction[0] = exp(-(pow(0.5-distances[mindist],2.)+
			pow(0.5-distances[secmindist],2.))/
		      pow(det->split->par1,2.));
    fraction[1] = exp(-(pow(0.5+distances[mindist],2.)+
			pow(0.5-distances[secmindist],2.))/
		      pow(det->split->par1,2.));
    fraction[2] = exp(-(pow(0.5-distances[mindist],2.)+
			pow(0.5+distances[secmindist],2.))/
		      pow(det->split->par1,2.));
    fraction[3] = exp(-(pow(0.5+distances[mindist],2.)+
			pow(0.5+distances[secmindist],2.))/
		      pow(det->split->par1,2.));
    // Normalization to 1.
    double sum = fraction[0]+fraction[1]+fraction[2]+fraction[3];
    fraction[0] /= sum;
    fraction[1] /= sum;
    fraction[2] /= sum;
    fraction[3] /= sum;

    // END of exponential split model.

  } else {
    SIXT_ERROR("split model not supported");
    *status=EXIT_FAILURE;
    return(0);
  }


  // Add signal to all valid pixels of the split event.
  int ii, nvalidpixels=0;
  for(ii=0; ii<npixels; ii++) {
    if ((x[ii]>=0) && (x[ii]<det->pixgrid->xwidth) &&
	(y[ii]>=0) && (y[ii]<det->pixgrid->ywidth)) {
      addGenDetCharge2Pixel(det->line[y[ii]], x[ii], signal*fraction[ii], 
			    ph_id, src_id);
      nvalidpixels++;

      // Call the event trigger routine.
      if (GENDET_EVENT_TRIGGERED==det->readout_trigger) {
	GenDetReadoutPixel(det, y[ii], y[ii], x[ii], time, elf, status);
	CHECK_STATUS_BREAK(*status);

	// In event-triggered mode each event occupies its own frame.
	det->clocklist->frame++;
      }
    }
  }
  CHECK_STATUS_RET(*status, nvalidpixels);


  // Return the number of affected pixels.
  return(nvalidpixels);
}


void setGenDetStartTime(GenDet* const det, const double t0)
{
  det->clocklist->time         = t0;
  det->clocklist->readout_time = t0;
}

