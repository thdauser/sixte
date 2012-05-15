#include "phgen.h"


int phgen(AttitudeCatalog* const ac,
	  SourceCatalog* const srccat,
	  const double t0, 
	  const double exposure,
	  const double mjdref,
	  const double dt,
	  const float fov,
	  Photon* const ph,
	  int* const status)
{
  // Photon list buffer.
  static LinkedPhoListElement* pholist=NULL;
  // Counter for the photon IDs.
  static long long ph_id=0;

  // Current time.
  static double time=0.;
  if (time<t0) {
    time=t0;
  }

  // If the photon list is empty generate new photons from the 
  // given source catalog.
  while((NULL==pholist)&&(time<t0+exposure)) {
    // Determine the telescope pointing at the current point of time.
    Vector pointing = getTelescopeNz(ac, time, status);
    CHECK_STATUS_BREAK(*status);

    // Display the program progress status.
    double ra, dec;
    calculate_ra_dec(pointing, &ra, &dec);
    headas_chat(5, "\rtime: %.1lf s, telescope: (%.3lf,%.3lf)     ", 
		time, ra*180./M_PI, dec*180./M_PI);
    fflush(NULL);
    
    // Get photons for all sources in the catalog.
    double t1 = MIN(time+dt, t0+exposure);
    pholist =
      genFoVXRayPhotons(srccat, &pointing, fov,
			time, t1, mjdref, status);
    CHECK_STATUS_BREAK(*status);
    
    // Increase the time.
    time+=dt;
  }

  // If there is no photon in the buffer.
  if (NULL==pholist) return(0);

  // Take the first photon from the list and return it.
  copyPhoton(ph, &pholist->photon);

  // Delete the processed element.
  LinkedPhoListElement* next=pholist->next;
  free(pholist);
  pholist = next;

  // Set the photon ID.
  ph->ph_id=++ph_id;

  return(1);
}

