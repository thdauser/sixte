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
#include "simput.h"

#define TOOLSUB ero_gti_main
#include "headas_main.c"


/** Margin around the FOV [rad]. */
const double margin=0.5*M_PI/180.; 


/* Program parameters */
struct Parameters {
  /** Attitude file. */
  char Attitude[MAXFILENAME];

  /** Source catalog. */
  char Simput[MAXFILENAME];

  /** GTI file. */
  char GTIfile[MAXFILENAME];
  
  double TIMEZERO;
  double Exposure;
  double dt; 

  /** [rad]. */
  double fov_diameter;

  char clobber;
};


int ero_gti_getpar(struct Parameters *parameters);


int ero_gti_main() 
{
  // Program parameters.
  struct Parameters par;

  AttitudeCatalog* ac=NULL;
  SimputCatalog* cat=NULL;
  GTI* gti=NULL;

  // Error status.
  int status=EXIT_SUCCESS;


  // Register HEATOOL:
  set_toolname("ero_gti");
  set_toolversion("0.02");
  

  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    // Read the program parameters using PIL library.
    if ((status=ero_gti_getpar(&par))) break;

    // Get the telescope attitude data.
    ac=loadAttitudeCatalog(par.Attitude, &status);
    CHECK_STATUS_BREAK(status);

    // Load the SIMPUT source catalog.
    cat=openSimputCatalog(par.Simput, READONLY, 0, 0, 0, 0, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the catalog contains any sources.
    if (0==cat->nentries) {
      SIXT_ERROR("SIMPUT catalog is empty");
      status=EXIT_FAILURE;
      break;
    }

    // Set up a new GTI collection.
    gti=newGTI(&status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---


    // --- Beginning of source localization calculation ---

    headas_chat(3, "determine the search cone ...\n");

    // Determine a reference vector and a cone opening angle
    // including all sources in the catalog.
    Vector refpos;
    double cone_radius=0.; // [rad]

    // Get the first source in the catalog.
    SimputSource* src=loadCacheSimputSource(cat, 1, &status);
    CHECK_STATUS_BREAK(status);

    // Determine its position and angular extension.
    refpos=unit_vector(src->ra, src->dec);
    cone_radius=getSimputSourceExtension(cat, src, &status);
    CHECK_STATUS_BREAK(status);

    // Loop over all sources in the catalog.
    long ii;
    for (ii=1; ii<cat->nentries; ii++) {
      // Get the next source in the catalog.
      SimputSource* src=loadCacheSimputSource(cat, ii+1, &status);
      CHECK_STATUS_BREAK(status);

      // Determine its position and angular extension.
      Vector srcpos=unit_vector(src->ra, src->dec);
      float extension=getSimputSourceExtension(cat, src, &status);
      CHECK_STATUS_BREAK(status);

      // Determine the angle between the reference direction and 
      // the source position.
      double angular_distance=acos(scalar_product(&refpos, &srcpos));
      
      // Check if the search radius has to be enlarged.
      if (angular_distance+extension>cone_radius) {
	double delta1=extension+angular_distance-cone_radius;

	// If the difference is small, simply enlarge the opening
	// angle of the search cone.
	if (delta1 < 1./3600.*M_PI/180.) {
	  cone_radius+=delta1;
	} else {
	  double delta2=extension-angular_distance-cone_radius;
	  if (delta2<0.) {
	    delta2=0.;
	  }
	  refpos=interpolateCircleVector(refpos, srcpos, 
					 (delta1-delta2)*0.5/angular_distance);
	  cone_radius+=0.5*(delta1+delta2);
	}
      }
    }

    // Print some informational data.
    double ra, dec;
    calculate_ra_dec(refpos, &ra, &dec);
    headas_chat(5, "ra=%lf deg, dec=%lf deg, cone radius=%lf deg\n", 
		ra*180./M_PI, dec*180./M_PI, cone_radius*180./M_PI);

    // Determine the diameter of the search radius (minimum cos-value).
    // (angle(telescope,source) <= 1/2 * diameter)
    double search_angle=0.5*par.fov_diameter+cone_radius+margin;
    double min_align; 
    if (search_angle <= M_PI) {
      min_align=cos(search_angle);
    } else {
      min_align = -1.; 
    }      

    // --- Beginning of GTI calculation ---

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

      // Check if the FOV touches the search cone.
      if (check_fov(&refpos, &telescope_nz, min_align)==0) {
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
  
  // Get the filename of the SIMPUT file.
  else if ((status = PILGetFname("Simput", par->Simput))) {
    SIXT_ERROR("failed reading the name of the SIMPUT file");
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

  else if ((status=ape_trad_query_bool("clobber", &par->clobber))) {
    SIXT_ERROR("failed reading the clobber parameter");
  }

  // Convert FOV diameter from [arc min] to [rad].
  par->fov_diameter *= M_PI/(60.*180.); 
  
  return(status);
}

