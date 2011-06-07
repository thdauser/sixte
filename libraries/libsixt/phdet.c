#include "phdet.h"


void phdetGenDet(GenDet* const det,
		 ImpactListFile* const ilf,
		 EventListFile* const elf,
		 const double t0,
		 const double exposure,
		 int* const status)
{
  // Total number of detected photons. Only the number of
  // photons absorbed by valid pixels inside the detector is
  // counted. Split events created by one photon are counted only
  // once.
  long n_detected_photons=0;

  // Set the start time for the detector simulation.
  setGenDetStartTime(det, t0);

  // Loop over all impacts in the FITS file.
  while (ilf->row<ilf->nrows) {

    Impact impact;
    getNextImpactFromFile(ilf, &impact, status);
    CHECK_STATUS_VOID(*status);

    // Check whether we are still within the requested time interval.
    if (impact.time < t0) continue;
    if (impact.time > t0+exposure) break;

    // Add the impact to the detector array. If it is absorbed
    // by at least one valid pixel, increase the counter for
    // the number of detected photons.
    if (addGenDetPhotonImpact(det, &impact, elf, status) > 0) {
      n_detected_photons++;
    }
    CHECK_STATUS_BREAK(*status);

    if (0==ilf->row%1000) {
      headas_chat(2, "\r %ld of %ld impacts (%.2lf%%) ", 
		  ilf->row, ilf->nrows, ilf->row*100./ilf->nrows);
      fflush(NULL);
    }

  };
  CHECK_STATUS_VOID(*status);
  // END of loop over all impacts in the FITS file.
    
  // Finalize the GenDet. Perform the time-triggered operations 
  // without adding any new charges.
  operateGenDetClock(det, elf, t0+exposure, status);
  CHECK_STATUS_VOID(*status);

  // Store the number of simulated input photons in the FITS header
  // of the output event file.
  fits_update_key(elf->fptr, TLONG, "NPHOTONS", 
		  &ilf->nrows, "number of input photons", 
		  status);
  CHECK_STATUS_VOID(*status);

  // Store the number of detected photons in the FITS header of
  // the output event file.
  fits_update_key(elf->fptr, TLONG, "NDETECTD", 
		  &n_detected_photons, "number of detected photons", 
		  status);
  CHECK_STATUS_VOID(*status);

}

