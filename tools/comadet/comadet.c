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
  
  ImpactListFile* ilf=NULL;
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
    ilf=openImpactListFile(par.ImpactList, READONLY, &status);
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
      status=addImpact2CoMaDetector(detector, &impact);
      CHECK_STATUS_RET(status, status);
      
    } //END of scanning the impact list.
    CHECK_STATUS_RET(status, status);

  } while(0);  // END of the error handling loop.


  // --- cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  //Free the Coded Mask Detector.
  freeCoMaDetector(detector);

  //Close the FITS files.
  freeImpactListFile(&ilf, &status);

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



