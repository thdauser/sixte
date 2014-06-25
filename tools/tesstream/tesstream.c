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


   Copyright 2014 Jelle de Plaa, SRON, Thorsten Brand, FAU
*/

#include "tesstream.h"

////////////////////////////////////
/** Main procedure. */
int tesstream_main() {
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  
  // Error status.
  int status=EXIT_SUCCESS;
  int Nstreams=0;
  
  PixImpFile* impfile=NULL;
  TESProfiles* profiles=NULL;
  TESFitsStream** fitsstream=NULL;
  TESDataStream* stream=NULL;
  AdvDet *det=NULL;
  
  fitsfile *ofptr=NULL;
  
  int *activearray=NULL;
  
  // Keywords
  char telescop[MAXMSG];
  char instrume[MAXMSG];
  char filter[MAXMSG];
  char ancrfile[MAXMSG];
  char respfile[MAXMSG];
  double mjdref;
  double timezero;
  double tstart;
  double tstop;
  
  int ii;
  
  // Register HEATOOL:
  set_toolname("tesstream");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at 
       // most be run once).
       
    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);
       
    // Open the pixel impact file
    impfile=openPixImpFile(par.PixImpList, READONLY, &status);
    CHECK_STATUS_BREAK(status);
    
    // Read keywords from input file
    sixt_read_fits_stdkeywords(impfile->fptr,
			       telescop,
			       instrume,
			       filter,
			       ancrfile,
			       respfile,
			       &mjdref,
			       &timezero,
			       &tstart,
			       &tstop, 
			       &status);
    CHECK_STATUS_BREAK(status);
    printf("Pixel impact file reaches from %lfs-%lfs .\n", tstart, tstop);
    if(tstart>par.tstart){
      puts("Program parameter tstart smaller as in input file.");    
    }else{
      tstart=par.tstart;
    }
    if(tstop<par.tstop){
      puts("Program parameter tstop larger as in input file.");
    }else{
      tstop=par.tstop;
    }
    printf("Simulate from %lfs-%lfs .\n", tstart, tstop);
  
    // Load the detector structure
    det=loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    
    // construct array of active pixels
    activearray=(int*)malloc(det->npix*sizeof(int));
    if(activearray==NULL){
      status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for array of actove pixels failed");
      CHECK_STATUS_BREAK(status);
    }
    if(par.Nactive==-1){
      par.Nactive=det->npix;
      par.nlo=0;
      par.nhi=det->npix-1;
    }
    int act=0;
    for(ii=0; ii<det->npix; ii++){
      if(ii<par.nlo || ii>par.nhi){
	activearray[ii]=-1;
      } else {
	activearray[ii]=act;
	act++;
      }
    }
    if(act!=par.Nactive){
      status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "Number of active pixels has been corrupted. Counted %d instead of %d pixels.", act, par.Nactive);
      SIXT_ERROR(msg);
      CHECK_STATUS_BREAK(status);
    }
    
    profiles=newTESProfiles(&status);
    CHECK_STATUS_BREAK(status);
    
    for(ii=0; ii<det->npix; ii++){
      // Test if profile is already loaded, if not, load it
      int versionindex=findTESProfileVersionIndex(profiles, det->pix[ii].version);
      if(activearray[ii]>-1 && versionindex<0){
	char profilename[MAXFILENAME];
	sprintf(profilename, "%s%s", det->filepath, det->pix[ii].tesproffilename);
	readTESProfiles(profilename,
			det->pix[ii].version, 
			profiles, 
			&status);
	CHECK_STATUS_BREAK(status);
	versionindex=findTESProfileVersionIndex(profiles, det->pix[ii].version);
	if(versionindex<0){
	  status=EXIT_FAILURE;
	  SIXT_ERROR("New profile not loaded.");
	  CHECK_STATUS_BREAK(status);
	}
	det->pix[ii].profVersionID=versionindex;
      } else {
	det->pix[ii].profVersionID=versionindex;
      }
    }
    if(status!=EXIT_SUCCESS){
      SIXT_ERROR("Failed reading pulse profile templates.");
      CHECK_STATUS_BREAK(status);
    }
    
  
    // Generate the data      
    stream=newTESDataStream(&status);
    CHECK_STATUS_BREAK(status);
    
    getTESDataStream(stream, 
		     impfile, 
		     profiles,
		     det,
		     tstart, 
		     tstop,
		     det->npix,
		     par.Nactive,
		     activearray,
		     &status);
    CHECK_STATUS_BREAK(status);
    
    // Get output data arrays
    Nstreams=(par.Nactive+TESFITSMAXPIX-1)/TESFITSMAXPIX;
    fitsstream=(TESFitsStream**)malloc(Nstreams*sizeof(TESFitsStream*));
    if(fitsstream==NULL){
      status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for FITS streams failed");
      CHECK_STATUS_BREAK(status);
    }
    for(ii=0; ii<Nstreams; ii++){
      fitsstream[ii]=newTESFitsStream(&status);
      CHECK_STATUS_BREAK(status);
    }
    createTESFitsStreamFile(&ofptr, 
			    par.streamname,
			    telescop,
			    instrume,
			    filter,
			    ancrfile,
			    respfile,
			    par.XMLFile,
			    par.PixImpList,
			    mjdref,
			    timezero,
			    tstart,
			    tstop,
			    par.clobber, 
			    &status);
    CHECK_STATUS_BREAK(status);
  
    // Write the output file
    
    int restpix=par.Nactive;
    for(ii=0; ii<Nstreams; ii++){
      // print extension name
      sprintf(fitsstream[ii]->name, "ADC%03d", ii+1);
      
      // calculate number of pixels in the extension
      int extpix=TESFITSMAXPIX;
      if(extpix>restpix){
	extpix=restpix;
      }
      
      // allocate memory for fits stream
      allocateTESFitsStream(fitsstream[ii], 
			   stream->Ntime, 
			   extpix, 
			   &status);
      CHECK_STATUS_BREAK(status);
      
      // transfer values into fits stream
      int pp, tt, ll, id;
      // time array
      for(tt=0; tt<stream->Ntime; tt++){
	fitsstream[ii]->time[tt]=stream->time[tt];
      }
      // find the ofiginal pixel id (reverse of the activearray)
      for(pp=0; pp<extpix; pp++){
	id=-1;
	for(ll=0; ll<det->npix; ll++){
	  if(activearray[ll]==par.Nactive-restpix+pp){
	    id=ll;
	  }
	}
	if(id==-1){
	  status=EXIT_FAILURE;
	  SIXT_ERROR("Pixel ID not found.");
	  CHECK_STATUS_BREAK(status);
	}
	CHECK_STATUS_BREAK(status);
	fitsstream[ii]->pixID[pp]=id;

	// adc_value array
	for(tt=0; tt<stream->Ntime; tt++){
	  fitsstream[ii]->adc_value[pp][tt]=stream->adc_value[tt][par.Nactive-restpix+pp];
	}
      }
      CHECK_STATUS_BREAK(status);
  
      double timeres=1./det->SampleFreq;
      writeTESFitsStream(ofptr, 
			  fitsstream[ii],
			  tstart,
			  tstop,
			  timeres,
			  &status);
      CHECK_STATUS_BREAK(status);
      restpix=restpix-TESFITSMAXPIX;
    }
    CHECK_STATUS_BREAK(status);
    
    fits_close_file(ofptr, &status);
    CHECK_STATUS_BREAK(status);
       
  } while(0); // END of the error handling loop.
  
  destroyTESProfiles(profiles);
  destroyAdvDet(&det);
  freePixImpFile(&impfile, &status);
  destroyTESDataStream(stream);
  for(ii=0; ii<Nstreams; ii++){
    destroyTESFitsStream(fitsstream[ii]);
  }
  free(fitsstream);
  free(activearray);
  
  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
  
}

int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("PixImpList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pixel impact list");
    return(status);
  } 
  strcpy(par->PixImpList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Streamfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the Stream file");
    return(status);
  } 
  strcpy(par->streamname, sbuffer);
  free(sbuffer); 
  
  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }
  
  status=ape_trad_query_double("tstart", &par->tstart);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstart parameter");
    return(status);
  }
  
  status=ape_trad_query_double("tstop", &par->tstop);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstop parameter");
    return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }
  
  char *pix=NULL;
  status=ape_trad_query_string("pixels", &pix);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pixel parameter");
    return(status);
  }
  if(!strcmp(pix, "all")){
    par->Nactive=-1;
  } else {
    char *ptr=strchr(pix, '-');
    if(ptr!=NULL){
      sscanf(pix, "%d - %d", &(par->nlo), &(par->nhi));
      par->nlo--;
      par->nhi--;
      par->Nactive=par->nhi-par->nlo+1;
      printf("Activate pixels %d-%d (%d pixels).\n", par->nlo+1, par->nhi+1, par->Nactive);
    } else {
      sscanf(pix, "%d", &(par->nlo));
      par->nlo--;
      par->nhi=par->nlo;
      par->Nactive=1;
    }
    if(par->Nactive<0 || par->nlo<0){
      status=EXIT_FAILURE;
      SIXT_ERROR("Input values for active pixels corrupted.");
      return(status);
    }
  }
  free(pix);

  return(status);
}
