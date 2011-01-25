#include "interface.h"


void photon_generation(const char* const xml_filename,
		       const char* const attitude_filename,
		       const char* const simput_filename,
		       const char* const photon_filename,
		       const double t0, const double t1,
		       int* const status)
{
  // Step width of the time loop.
  const double dt = 1.0;

  // Template for the photon list file.
  char photonlist_template[] = "";

  // Read the detector definition.
  GenDet* det = newGenDet(xml_filename, status);
  CHECK_STATUS_VOID(*status);

  // Read the attitude.
  AttitudeCatalog* ac = get_AttitudeCatalog(attitude_filename,
					    t0, t1-t0, status);
  CHECK_STATUS_VOID(*status);

  // Load the SIMPUT X-ray source catalog.
  XRaySourceCatalog* srccat = 
    loadSourceCatalog(simput_filename, det, status);
  CHECK_STATUS_VOID(*status);
  
  // Remove the old photon list file.    
  remove(photon_filename);
  
  // Open the output photon list file.
  PhotonListFile plf = { .fptr=NULL };
  *status=openNewPhotonListFile(&plf, photon_filename, photonlist_template);
  CHECK_STATUS_VOID(*status);

  // Loop over the specified time interval.
  double time;
  for (time=t0; time<t1; time+=dt) {

    // Determine the telescope pointing at the current point of time.
    Vector pointing =  getTelescopePointing(ac, time, status);
    CHECK_STATUS_BREAK(*status);
    
    // Get photons for all sources in the catalog.
    LinkedPhoListElement* pholist =
      genFoVXRayPhotons(srccat, &pointing, det->fov_diameter,
			time, time+dt, status);
    CHECK_STATUS_BREAK(*status);

    // Process the list of generated photons.
    while(pholist) {
      // Output to the FITS file.
      *status = addPhoton2File(&plf, &(pholist->photon));
      CHECK_STATUS_BREAK(*status);

      // Delete the processed element.
      LinkedPhoListElement* next = pholist->next;
      free(pholist);
      pholist = next;
    }
    CHECK_STATUS_BREAK(*status);
  }
  CHECK_STATUS_VOID(*status);
  // END of loop over the requested time interval.

  // --- Clean up ---
  closePhotonListFile(&plf);
  freeXRaySourceCatalog(&srccat);
  free_AttitudeCatalog(ac);
  destroyGenDet(&det, status);
}

