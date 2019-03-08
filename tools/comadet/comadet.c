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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "comadet.h"

/////////////////////////////////////////////////////////////
//DETECTION:generates corresp. pixel-values in det-plane   //
//discards those photons hitting det-gaps                  //
//Input:impact-list(t,E,x,y,PH_ID,SRC_ID)                  //
//      Width(of detector in pixels), Pixelwidth           //
//      length of DCU, DCU_gap, DCa_gap                    //
//Output:event-list(pixel-values RAWX,RAWY,t,charge)       //
/////////////////////////////////////////////////////////////

////////////////////////////////////
/** Main procedure. */
int comadet_main() {
  struct Parameters par;

  ImpactFile* ilf=NULL;
  CoMaDetector* detector=NULL;

  //Error status.
  int status=EXIT_SUCCESS;


  //Register HEATOOL:
  set_toolname("comadet");
  set_toolversion("0.02");


  do {  //Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    //Read the program parameters using the PIL library.
    status=comadet_getpar(&par);
    CHECK_STATUS_RET(status, status);

    //Open the impact list FITS file.
    ilf=openImpactFile(par.ImpactList, READONLY, &status);
    CHECK_STATUS_RET(status, status);

    //Set the event list template file:
    strcpy(par.EventListTemplate, SIXT_DATA_PATH);
    strcat(par.EventListTemplate, "/templates/coma.eventlist.tpl");

    //DETECTOR setup.
    //initializes from par-file
    struct CoMaDetectorParameters cdp = {
      .pixels =
      { .xwidth = par.width,          //detector-width [pixel]
	.ywidth = par.width,
	.DCU_length = par.DCU_length, //length of DCU [m]
	.DCU_gap = par.DCU_gap,       //gap between 2 DCU's [m]
	.DCA_gap = par.DCA_gap,       //gap between 2 DCA's [m]
	.xpixelwidth = par.pixelwidth,//width of one pixel [m]
	.ypixelwidth = par.pixelwidth
      },
      .eventfile_filename=par.EventList,
      .eventfile_template=par.EventListTemplate
    };
    //create new CoMaDetector-object:
    detector=getCoMaDetector(&cdp, &status);
    CHECK_STATUS_RET(status, status);

    //END of DETECTOR CONFIGURATION SETUP

    // --- END of Initialization ---


    // --- Beginning of Imaging Process ---

    //Beginning of actual detector simulation (after loading required data):
    headas_chat(5, "start detection process ...\n");
    Impact impact;

    //Loop over all impacts in the FITS file.
    while (ilf->row<ilf->nrows) {

      getNextImpactFromFile(ilf, &impact, &status);
      CHECK_STATUS_RET(status, status);

      //detection routine:determines affected pixel and adds new event
      //to the event-file
      if(par.protoMirax==1){
	status=addImpact2CoMaDetector_protoMirax(detector, &impact);
	CHECK_STATUS_RET(status, status);
      }else{
	status=addImpact2CoMaDetector(detector, &impact);
	CHECK_STATUS_RET(status, status);
      }

    } //END of scanning the impact list.
    CHECK_STATUS_RET(status, status);

  } while(0);  // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  //Free the Coded Mask Detector.
  freeCoMaDetector(detector);

  //Close the FITS files.
  freeImpactFile(&ilf, &status);

  if (EXIT_SUCCESS==status) headas_chat(5, "finished successfully!\n\n");
  return(status);
}



int comadet_getpar(struct Parameters* par)
{
  int status=EXIT_SUCCESS; // Error status.

  // Get the filename of the impact list file (FITS input file).
  if ((status=PILGetFname("Impactlist", par->ImpactList))) {
    SIXT_ERROR("failed reading the filename of the impact list");
  }

  // Get the filename of the event list file (FITS output file).
  else if ((status=PILGetFname("EventList", par->EventList))) {
    SIXT_ERROR("failed reading the filename of the event list output");
  }

  // Check whether protoMirax or not [yes=1, no=0].
  else if ((status=PILGetInt("protoMirax", &par->protoMirax))) {
    SIXT_ERROR("failed reading protoMirax-flag");
  }

  // Read the width of the detector in [pixel].
  else if ((status=PILGetInt("Width", &par->width))) {
    SIXT_ERROR("failed reading the detector width");
  }

  // Read the width of one detector pixel in [m].
  else if ((status=PILGetReal("Pixelwidth", &par->pixelwidth))) {
    SIXT_ERROR("failed reading the width of the detector pixels");
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
  CHECK_STATUS_RET(status, status);

  return(status);
}
