#include "phdet.h"


void phdetGenDet(GenDet* const det,
		 Impact* const impact,
		 EventListFile* const elf,
		 const double t0,
		 const double exposure,
		 int* const status)
{
  // Total number of detected photons. Only the number of
  // photons absorbed by valid pixels inside the detector is
  // counted. Split events created by one photon are counted only
  // once.
  static unsigned long n_detected_photons=0;

  // Check if an impact has been given as a parameter.
  if (NULL!=impact) {
    // Add the impact to the detector array. If it is absorbed
    // by at least one valid pixel, increase the counter for
    // the number of detected photons.
    if (addGenDetPhotonImpact(det, impact, elf, status) > 0) {
      n_detected_photons++;
    }
    CHECK_STATUS_VOID(*status);

  } else {
    // If no impact has been given as parameter, finalize the GenDet. 
    // Perform the time-triggered operations without adding any new 
    // signal charges.
    operateGenDetClock(det, elf, t0+exposure, status);
    CHECK_STATUS_VOID(*status);

    // Store the number of detected photons in the FITS header of
    // the output event file.
    fits_update_key(elf->fptr, TLONG, "NDETECTD", 
		    &n_detected_photons, "number of detected photons", 
		    status);
    CHECK_STATUS_VOID(*status);
  }
}
