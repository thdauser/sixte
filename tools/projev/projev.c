#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif


#include "sixt.h"
#include "attitudecatalog.h"
#include "event.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "point.h"
#include "vector.h"

#define TOOLSUB projev_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  char EventList[MAXFILENAME];
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  char XMLFile[MAXFILENAME];
  char Attitude[MAXFILENAME];

  /** [deg] */
  float RA, Dec;

  double MJDREF;
  double TIMEZERO;
  double Exposure;

  int Seed;
  
  char clobber;

  char data_path[MAXFILENAME];
};


int projev_getpar(struct Parameters *par);


////////////////////////////////////
/** Main procedure. */
int projev_main() {
  struct Parameters par; // Program parameters

  // Input event list file.
  EventListFile* elf=NULL;
  // Attitude catalog.
  AttitudeCatalog* ac=NULL;
  // Detector data structure (containing the focal length, FoV, ...).
  GenDet* det=NULL;

  // Error status.
  int status=EXIT_SUCCESS;   


  // Register HEATOOL:
  set_toolname("projev");
  set_toolversion("0.01");


  do {  // Beginning of the ERROR handling loop (will at most be run once)

    // --- Initialization ---

    // Read the program parameters using PIL library.
    status=projev_getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");

    // Start time for the simulation.
    double t0 = par.MJDREF*24.*3600. + par.TIMEZERO;

    // Determine the random number seed.
    int seed;
    if (-1!=par.Seed) {
      seed = par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed = (int)time(NULL);
    }

    // Initialize HEADAS random number generator.
    HDmtInit(seed);

    // Determine the appropriate detector XML definition file.
    char xml_filename[MAXFILENAME];
    // Convert the user input to capital letters.
    strtoupper(par.Mission);
    strtoupper(par.Instrument);
    strtoupper(par.Mode);
    // Check the available missions, instruments, and modes.
    char ucase_buffer[MAXFILENAME];
    strcpy(ucase_buffer, par.XMLFile);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Determine the base directory containing the XML
      // definition files.
      strcpy(xml_filename, par.data_path);
      strcat(xml_filename, "/instruments");

      // Determine the XML filename according to the selected
      // mission, instrument, and mode.
      if (0==strcmp(par.Mission, "SRG")) {
	strcat(xml_filename, "/srg");
	if (0==strcmp(par.Instrument, "EROSITA")) {
	  strcat(xml_filename, "/erosita.xml");
	} else {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("selected instrument is not supported");
	  break;
	}

      } else if (0==strcmp(par.Mission, "IXO")) {
	strcat(xml_filename, "/ixo");
	if (0==strcmp(par.Instrument, "WFI")) {
	  strcat(xml_filename, "/wfi");
	  if (0==strcmp(par.Instrument, "FULLFRAME")) {
	    strcat(xml_filename, "/fullframe.xml");
	  } else {
	    status=EXIT_FAILURE;
	    SIXT_ERROR("selected mode is not supported");
	    break;
	  }
	} else {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("selected instrument is not supported");
	  break;
	}

      } else if (0==strcmp(par.Mission, "GRAVITAS")) {
	strcat(xml_filename, "/gravitas");
	if (0==strcmp(par.Instrument, "HIFI")) {
	  strcat(xml_filename, "/hifi.xml");
	} else {
	  status=EXIT_FAILURE;
	  SIXT_ERROR("selected instrument is not supported");
	  break;
	}

      } else {
	status=EXIT_FAILURE;
	SIXT_ERROR("selected mission is not supported");
	break;
      }
	    
    } else {
      // The XML filename has been given explicitly.
      strcpy(xml_filename, par.XMLFile);
    }
    // END of determine the XML filename.

    // Load the detector configuration.
    det=newGenDet(xml_filename, &status);
    CHECK_STATUS_BREAK(status);

    // Set up the Attitude.
    strcpy(ucase_buffer, par.Attitude);
    strtoupper(ucase_buffer);
    if (0==strcmp(ucase_buffer, "NONE")) {
      // Set up a simple pointing attitude.

      // First allocate memory.
      ac=getAttitudeCatalog(&status);
      CHECK_STATUS_BREAK(status);

      ac->entry=(AttitudeEntry*)malloc(2*sizeof(AttitudeEntry));
      if (NULL==ac->entry) {
	status = EXIT_FAILURE;
	SIXT_ERROR("memory allocation for AttitudeCatalog failed");
	break;
      }

      // Set the values of the entries.
      ac->nentries=2;
      ac->entry[0] = defaultAttitudeEntry();
      ac->entry[1] = defaultAttitudeEntry();
      
      ac->entry[0].time = t0;
      ac->entry[1].time = t0 + par.Exposure;

      ac->entry[0].nz = unit_vector(par.RA*M_PI/180., par.Dec*M_PI/180.);
      ac->entry[1].nz = ac->entry[0].nz;

      Vector vz = {0., 0., 1.};
      ac->entry[0].nx = vector_product(vz, ac->entry[0].nz);
      ac->entry[1].nx = ac->entry[0].nx;

    } else {
      // Load the attitude from the given file.
      ac=loadAttitudeCatalog(par.Attitude, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the required time interval for the simulation
      // is a subset of the time described by the attitude file.
      if ((ac->entry[0].time > t0) || 
	  (ac->entry[ac->nentries-1].time < t0+par.Exposure)) {
	status=EXIT_FAILURE;
	char msg[MAXMSG];
	sprintf(msg, "attitude data does not cover the "
		"specified period from %lf to %lf!", t0, t0+par.Exposure);
	HD_ERROR_THROW(msg, status);
	break;
      }
    }
    // END of setting up the attitude.

    // Set the input event file.
    elf=openEventListFile(par.EventList, READWRITE, &status);
    if (EXIT_SUCCESS!=status) break;


    // For very small angles tan(x) \approx x.
    float radpermeter = 1./det->focal_length; 

    // --- END of Initialization ---


    // --- Beginning of the Sky Projection Process ---

    // Beginning of actual simulation (after loading required data):
    headas_chat(5, "start sky projection process ...\n");

    // LOOP over all events in the FITS table.
    long row;
    for (row=0; row<elf->nrows; row++) {

      // Read the next event from the file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      // Determine the Position of the source on the sky:
      // First determine telescope pointing direction at the current time.
      Vector nx, ny, nz;
      getTelescopeAxes(ac, &nx, &ny, &nz, event.time, &status);
      CHECK_STATUS_BREAK(status);

      // Determine RA and DEC of the photon origin.
      // Exact position on the detector:
      struct Point2d detpos;
      detpos.x = // in [m]
	(event.rawx*1.-det->pixgrid->xrpix+0.5+sixt_get_random_number())*
	det->pixgrid->xdelt + 
	det->pixgrid->xrval;
      detpos.y = // in [m]
	(event.rawy*1.-det->pixgrid->yrpix+0.5+sixt_get_random_number())*
	det->pixgrid->ydelt + 
	det->pixgrid->yrval;
      double d = sqrt(pow(detpos.x,2.)+pow(detpos.y,2.)); // in [m]

      // Determine the off-axis angle corresponding to the detector position.
      double offaxis_angle = d * radpermeter; // [rad]

      // Determine the source position on the sky using the telescope 
      // axis pointing vector and a vector from the point of the intersection 
      // of the optical axis with the sky plane to the source position.
      // Determine the length of this vector (in the sky projection plane).
      double r = tan(offaxis_angle); 

      Vector srcpos;
      srcpos.x = nz.x + r*(detpos.x/d*nx.x+detpos.y/d*ny.x);
      srcpos.y = nz.y + r*(detpos.x/d*nx.y+detpos.y/d*ny.y);
      srcpos.z = nz.z + r*(detpos.x/d*nx.z+detpos.y/d*ny.z);
      srcpos = normalize_vector(srcpos);

      // Determine the equatorial coordinates RA and DEC
      // (RA and DEC are in the range [-pi:pi] and [-pi/2:pi/2] respectively).
      calculate_ra_dec(srcpos, &event.ra, &event.dec);

      // Update the data in the Event List FITS file.
      updateEventInFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);
    } 
    CHECK_STATUS_BREAK(status);
    // END of LOOP over all events

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(5, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Destroy the detector data structure.
  destroyGenDet(&det, &status);
  
  // Close the files.
  freeEventListFile(&elf, &status);

  // Release memory of AttitudeCatalog
  freeAttitudeCatalog(&ac);

  if (status == EXIT_SUCCESS) headas_chat(5, "finished successfully!\n\n");
  return(status);
}



int projev_getpar(struct Parameters* par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status = EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the event list!\n", status);
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mission", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the mission!\n", status);
    return(status);
  } 
  strcpy(par->Mission, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Instrument", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the instrument!\n", status);
    return(status);
  } 
  strcpy(par->Instrument, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Mode", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the instrument mode!\n", status);
    return(status);
  } 
  strcpy(par->Mode, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the XML file!\n", status);
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("Attitude", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the attitude!\n", status);
    return(status);
  } 
  strcpy(par->Attitude, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RA", &par->RA);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the right ascension of the telescope pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_float("Dec", &par->Dec);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the declination of the telescope pointing!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading MJDREF!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("TIMEZERO", &par->TIMEZERO);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading TIMEZERO!\n", status);
    return(status);
  } 

  status=ape_trad_query_double("Exposure", &par->Exposure);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  } 

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the seed for the random number generator!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }


  // Get the name of the directory containing the detector
  // XML definition files from the environment variable.
  if (NULL!=(sbuffer=getenv("SIXT_DATA_PATH"))) {
    strcpy(par->data_path, sbuffer);
    // Note: the char* pointer returned by getenv should not
    // be modified nor free'd.
  } else {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error reading the environment variable 'SIXT_DATA_PATH'!\n", 
		   status);
    return(status);
  }

  return(status);
}



