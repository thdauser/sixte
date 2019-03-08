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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "phabkg.h"

#include <math.h>


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


PHABkg* newPHABkg(const char* const filename, int* const status)
{
  // Allocate memory.
  PHABkg* phabkg=(PHABkg*)malloc(sizeof(PHABkg));
  if (NULL==phabkg) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for PHABkg failed");
    return(phabkg);
  }

  // Initialize all pointers with NULL.
  phabkg->channel     =NULL;
  phabkg->distribution=NULL;
  phabkg->vignetting  =NULL;
  phabkg->focal_length=NULL;

  // Initialize values.
  phabkg->nbins=0;
  phabkg->tnext=0.;

  // Load the specified PHA file.
  fitsfile *fptr=NULL;

  do { // Beginning of ERROR handling loop.

    headas_chat(5, "load PHA background spectrum from file '%s' ...\n", filename);

    // Open the file:
    fits_open_table(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);

    // Get the HDU type.
    int hdutype;
    fits_get_hdu_type(fptr, &hdutype, status);
    CHECK_STATUS_BREAK(*status);

    // Image HDU results in an error message.
    if (IMAGE_HDU==hdutype) {
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "no table extension available in FITS file '%s'",
	      filename);
      SIXT_ERROR(msg);
      break;
    }

    // Determine the number of rows in the file.
    fits_get_num_rows(fptr, &phabkg->nbins, status);
    CHECK_STATUS_BREAK(*status);

    // Allocate memory.
    phabkg->channel=(long*)malloc(phabkg->nbins*sizeof(long));
    CHECK_NULL_BREAK(phabkg->channel, *status,
		     "memory allocation for PHA background spectrum failed");
    phabkg->distribution=(float*)malloc(phabkg->nbins*sizeof(float));
    CHECK_NULL_BREAK(phabkg->distribution, *status,
		     "memory allocation for PHA background spectrum failed");

    // Determine the column numbers.
    int cchannel, crate;
    fits_get_colnum(fptr, CASEINSEN, "CHANNEL", &cchannel, status);
    fits_get_colnum(fptr, CASEINSEN,    "RATE",    &crate, status);
    CHECK_STATUS_BREAK(*status);

    // Read the data.
    int anynul=0;
    fits_read_col(fptr, TLONG, cchannel, 1, 1, phabkg->nbins,
		  0, phabkg->channel, &anynul, status);
    fits_read_col(fptr, TFLOAT, crate, 1, 1, phabkg->nbins,
		  0, phabkg->distribution, &anynul, status);
    CHECK_STATUS_BREAK(*status);

    // Sum up the rate values in order to obtain an accumulative distribution.
    long ii;
    for (ii=1; ii<phabkg->nbins; ii++) {
      phabkg->distribution[ii]+=phabkg->distribution[ii-1];
    }

  } while (0); // END of error handling loop

  // Clean up:
  if (fptr) fits_close_file(fptr, status);

  return(phabkg);
}


void destroyPHABkg(PHABkg** const phabkg)
{
  if (NULL!=*phabkg) {
    if (NULL!=(*phabkg)->channel) {
      free((*phabkg)->channel);
    }
    if (NULL!=(*phabkg)->distribution) {
      free((*phabkg)->distribution);
    }
    free(*phabkg);
    *phabkg=NULL;
  }
}


int getPHABkgEvent(PHABkg* const phabkg,
		   const float scaling,
		   const double tstart,
		   const double tstop,
		   double* const t,
		   long* const pha,
		   int* const status)
{
	// Check if everything is set up properly.
	if (NULL==phabkg) {
		SIXT_ERROR("no PHA background model specified");
		*status=EXIT_FAILURE;
		return(0);
	}

	if (NULL==phabkg->distribution) {
		SIXT_ERROR("no PHA background spectrum loaded");
		*status=EXIT_FAILURE;
		return(0);
	}

	// Determine the average event rate.
	double rate=phabkg->distribution[phabkg->nbins-1]*scaling;
	if (rate==0){
		SIXT_ERROR("background spectrum has rate of 0");
	    *status=EXIT_FAILURE;
	    return(0);
	}

	// Update to the start time, if time of next background event lies
	// before that.
	if (phabkg->tnext<=tstart) {
		phabkg->tnext=tstart+rndexp(1./rate, status);
		CHECK_STATUS_RET(*status, 0);
	}

	// Check if the time of the next background event lies before the
	// specified upper limit. In this case no event is produced.
	if (phabkg->tnext>tstop) {
		return(0);
	}

	// Set the time of the background event.
	*t=phabkg->tnext;

	// Determine the PHA value.
	double r=sixt_get_random_number(status)*phabkg->distribution[phabkg->nbins-1];
	CHECK_STATUS_RET(*status, 0);

	// Perform a binary search to obtain the corresponding channel.
	long min=0;
	long max=phabkg->nbins-1;
	long mid;
	while (max>min) {
		mid=(min+max)/2;
		if (phabkg->distribution[mid]<r) {
			min=mid+1;
		} else {
			max=mid;
		}
	}
	*pha=phabkg->channel[min];

	// Determine the time of the next background event.
	phabkg->tnext+=rndexp(1./rate, status);
	CHECK_STATUS_RET(*status, 0);

	return(1);
}


LinkedImpListElement* getPHABkglist(PHABkg* const phabkg, AdvDet* det, const float scaling,
		const double tstart, const double tend, int* const status){

	// Time ordered impact list
	LinkedImpListElement* newlist = NULL;

	// Points to the last (NULL) pointer in the linked list.
	LinkedImpListElement** list_next=&newlist;


	double time_current=tstart;
	do { //While there are events get and store them
		double next_time=0;
		long channel=0;
		int is_bkg_event=getPHABkgEvent(phabkg, scaling, time_current, tend,
				&next_time, &channel, status);

		if (!is_bkg_event){ //When no more events, get out and return the list
			break;
		}

		*list_next=newLinkedImpListElement(status);
		CHECK_STATUS_BREAK(*status);

		//Allocating the next element
		Impact* imp=&((*list_next)->impact);
		list_next = &((*list_next)->next);

		//Generate the random location in which this occurs
		//First we get a random pixel
		long int pixid=floor(sixt_get_random_number(status)*det->npix);

		//Get the associated value of x and y
		imp->position.x=det->sx+det->pix[pixid].sx+(sixt_get_random_number(status)-0.5)*det->pix[pixid].width;
		imp->position.y=det->sy+det->pix[pixid].sy+(sixt_get_random_number(status)-0.5)*det->pix[pixid].height;

		//Setting the value of the next events
		imp->time=next_time;
		imp->energy=getEBOUNDSEnergy(channel,det->pix[pixid].grades[0].rmf,status);
		imp->ph_id=0;
		imp->src_id=PHABKG_SRC_ID;
		time_current=next_time;
	} while(1);

	return(newlist);
}

int phabkggen(PHABkg* const phabkg[2], AdvDet* det, const float scaling, const double t0, const double tend,
		const double dt, Impact* impact,int* const status) {

	// Impact list buffer.
	static LinkedImpListElement* implist=NULL;

	// Counter for the photon IDs.
	static long long bkg_ph_id=0;
	// Current time.
	static double bkg_time=0.;
	if (bkg_time<=t0) {
	  bkg_time=t0;
	}

	// If the photon list is empty generate new photons from the
	// given source catalog.
	while((NULL==implist)&&(bkg_time<tend)) {
	  // Generate new photons for all specified catalogs.
	  double t1=MIN(bkg_time+dt, tend);
	  unsigned int ii;
	  for (ii=0; ii<MAX_PHABKG; ii++) { //Two files max
	    if (NULL==phabkg[ii]) continue;

	    // Get photons for all sources in the catalog, giving first rmf as
	    // only the channels, not the redistribution are used i.e. all rmfs
	    // are equivalent
	    LinkedImpListElement* newlist=getPHABkglist(phabkg[ii], det, scaling,
	    		bkg_time, t1, status);
	    CHECK_STATUS_BREAK(*status);

	    // Merge the photon lists.
	    implist=mergeLinkedImpLists(implist, newlist);
	  }

	  // Increase the time.
	  bkg_time+=dt;
	}

	// If there is no photon in the buffer.
	if (NULL==implist) return(0);

	// Take the first photon from the list and return it.
	copyImpact(impact, &implist->impact);

	// Delete the processed element.
	LinkedImpListElement* next=implist->next;
	free(implist);
	implist=next;

	// Set the photon ID.
	 impact->ph_id=--bkg_ph_id;

	  return(1);

}
