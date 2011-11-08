#include "phgen.h"


void phgen(AttitudeCatalog* const ac,
	   SourceCatalog* const srccat,
	   PhotonListFile* const plf,
	   const double t0, 
	   const double exposure,
	   const double mjdref,
	   const float fov,
	   int* const status)
{
  // Step width of the time loop.
  const double dt = 1.0;

  // Loop over the specified time interval.
  double time;
  for (time=t0; time<t0+exposure; time+=dt) {
    // Determine the telescope pointing at the current point of time.
    Vector pointing = getTelescopeNz(ac, time, status);
    CHECK_STATUS_BREAK(*status);

    // Display the program progress status.
    double ra, dec;
    calculate_ra_dec(pointing, &ra, &dec);
    headas_chat(3, "\rtime: %.1lf s, telescope: (%.3lf,%.3lf)     ", 
		time, ra*180./M_PI, dec*180./M_PI);
    fflush(NULL);
    
    // Get photons for all sources in the catalog.
    double t2 = MIN(time+dt, t0+exposure);
    LinkedPhoListElement* pholist =
      genFoVXRayPhotons(srccat, &pointing, fov,
			time, t2, mjdref, status);
    CHECK_STATUS_BREAK(*status);
      
    // Process the list of generated photons.
    while(pholist) {
      // Output to the FITS file.
      *status = addPhoton2File(plf, &(pholist->photon));
      CHECK_STATUS_BREAK(*status);
	
      // Delete the processed element.
      LinkedPhoListElement* next = pholist->next;
      free(pholist);
      pholist = next;
    }
    CHECK_STATUS_BREAK(*status);
  }
  // END of loop over the requested time interval.
}

