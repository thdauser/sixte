#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "sixt.h"
#include "vector.h"
#include "telescope.h"
#include "attitudecatalog.h"
#include "check_fov.h"
#include "gti.h"

#define TOOLSUB ero_gti_main
#include "headas_main.c"


/* Program parameters */
struct Parameters {
  /** Attitude file. */
  char Attitude[MAXFILENAME];

  /** GTI file. */
  char GTIfile[MAXFILENAME];
  
  double TIMEZERO;
  double Exposure;
  double dt; 

  /** Source position [rad]. */
  double ra, dec;  

  /** [rad]. */
  double fov_diameter;

  /** Diameter of the source extension [rad]. */
  double src_diameter;

  char clobber;
};


int ero_gti_getpar(struct Parameters *parameters);


int ero_gti_main() 
{
  // Program parameters.
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  GTI* gti=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ero_gti");
  set_toolversion("0.01");
  

  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    if ((status=ero_gti_getpar(&par))) break;

    // Get the telescope attitude data.
    ac=loadAttitudeCatalog(par.Attitude, &status);
    CHECK_STATUS_BREAK(status);

    // Determine a vector pointing towards the source.
    Vector srcpos=unit_vector(par.ra, par.dec);

    // Determine the diameter of the search radius (minimum cos-value).
    // (angle(telescope,source) <= 1/2 * diameter)
    const double min_align = cos(0.5*(par.fov_diameter+par.src_diameter)); 
    double field_min_align;
    if (0.5*(par.fov_diameter+par.src_diameter) > M_PI) {
      field_min_align = -1.; 
    }

    // Get a new GTI collection.
    gti=newGTI(&status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---


    // --- Beginning of GTI calculation

    headas_chat(3, "calculate the GTIs ...\n");

    // LOOP over the given time interval in steps of dt.
    double time;
    double start=0;
    double ininterval=0;
    for (time=par.TIMEZERO; time<par.TIMEZERO+par.Exposure; time+=par.dt) {
      
      // Print the current time (program status information for the user).
      headas_printf("\rtime: %.1lf s ", time);
      fflush(NULL);

      // Determine the telescope pointing direction at the current time.
      Vector telescope_nz=getTelescopeNz(ac, time, &status);
      CHECK_STATUS_BREAK(status);

      // Check if the source lies within the FOV.
      if (check_fov(&srcpos, &telescope_nz, min_align)==0) {
	// Source lies inside the FOV.
	if (0==ininterval) {
	  ininterval=1;
	  start=time;
	}
      } else {
	// Source lies outside the FOV.
	if (1==ininterval) {
	  ininterval=0;
	  appendGTI(gti, start, time, &status);
	  CHECK_STATUS_BREAK(status);
	}
      }
    }
    CHECK_STATUS_BREAK(status);
    // END of LOOP over the specified time interval.
    

    // Store the GTI in the output file.
    saveGTI(gti, par.GTIfile, par.clobber, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release memory.
  freeAttitudeCatalog(&ac);
  freeGTI(&gti);

  if (EXIT_SUCCESS==status) headas_chat(2, "finished successfully!\n\n");
  return(status);
}


int ero_gti_getpar(struct Parameters *par)
{
  int status=EXIT_SUCCESS; // Error status
  
  // Get the filename of the input attitude file (FITS file).
  if ((status = PILGetFname("Attitude", par->Attitude))) {
    SIXT_ERROR("failed reading the name of the attitude file");
  }
  
  // Get the filename of the output GTI file (FITS file).
  else if ((status = PILGetFname("GTIfile", par->GTIfile))) {
    SIXT_ERROR("failed reading the name of the GTI file");
  }

  // Read the diameter of the FOV (in arcmin).
  else if ((status = PILGetReal("fov_diameter", &par->fov_diameter))) {
    SIXT_ERROR("failed reading the diameter of the FOV");
  }

  // Get the start time.
  else if ((status = PILGetReal("TIMEZERO", &par->TIMEZERO))) {
    SIXT_ERROR("failed reading the TIMEZERO");
  }

  // Get the exposure time.
  else if ((status = PILGetReal("Exposure", &par->Exposure))) {
    SIXT_ERROR("failed reading the exposure time");
  }

  // Get the time step.
  else if ((status = PILGetReal("dt", &par->dt))) {
    SIXT_ERROR("failed reading the 'dt' parameter");
  }

  // Get the position of the source.
  else if ((status = PILGetReal("ra", &par->ra))) {
    SIXT_ERROR("failed reading the right ascension");
  }
  else if ((status = PILGetReal("dec", &par->dec))) {
    SIXT_ERROR("Error reading the declination");
  }

  // Read the diameter of the source [deg].
  else if ((status = PILGetReal("src_diameter", &par->src_diameter))) {
    SIXT_ERROR("failed reading the diameter of the source");
  }

  else if ((status=ape_trad_query_bool("clobber", &par->clobber))) {
    SIXT_ERROR("failed reading the clobber parameter");
  }

  // Convert angles from [deg] to [rad].
  par->ra  *= M_PI/180.;
  par->ra  *= M_PI/180.;
  par->src_diameter *= M_PI/180.;

  // Convert angles from [arc min] to [rad].
  par->fov_diameter *= M_PI/(60.*180.); 
  
  return(status);
}

