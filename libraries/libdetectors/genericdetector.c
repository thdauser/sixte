#include "genericdetector.h"


int initGenericDetector(GenericDetector* gd, struct GenericDetectorParameters* parameters) 
{
  int status = EXIT_SUCCESS;

  // Set the charge cloud dimensions.
  // Gaussian:
  gd->gcc.ccsigma =    parameters->ccsigma;
  gd->gcc.ccsize  = 3.*parameters->ccsigma;
  headas_chat(5, "Split events: assuming Gaussian Charge Cloud model:\n");
  headas_chat(5, " charge cloud size: %lf mum (3 sigma)\n", 
	      gd->gcc.ccsize*1.e6);
  // Exponential model for eROSITA according to Konrad Dennerl:
  gd->ecc.parameter = 0.355;
  
  // Set the event thresholds:
  gd->pha_threshold = parameters->pha_threshold;
  gd->energy_threshold = parameters->energy_threshold;

  // Read the detector RMF and EBOUNDS from the specified file and 
  // assign them to the Detector data structure.
  gd->rmf = loadRMF(parameters->rmf_filename, &status);
  if(EXIT_SUCCESS!=status) return(status);

  return(status);
}



void cleanupGenericDetector(GenericDetector* gd)
{
  // Actually one has to release the energy allocated for the internal
  // objects in the RMF data structure. But unfortunately there exists
  // no routine with this functionality in the HEASP library.

  if (NULL==gd->rmf) return;
  //freeRMF(gd->rmf);
  //gd->rmf=NULL;
}



struct RMF* loadRMF(char* filename, int* status) 
{
  fitsfile* fptr=NULL;

  struct RMF* rmf = (struct RMF*)malloc(sizeof(struct RMF));
  if (NULL==rmf) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for RMF!\n", *status);
    return(rmf);
  }

  // Load the RMF from the FITS file using the HEAdas RMF access routines
  // (part of libhdsp).
  fits_open_file(&fptr, filename, READONLY, status);
  if (*status!=EXIT_SUCCESS) return(rmf);
  
  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  if ((*status=ReadRMFMatrix(fptr, 0, rmf))!=EXIT_SUCCESS) return(rmf);

  // Print some information:
  headas_chat(5, "RMF loaded with %ld energy bins and %ld channels\n",
	      rmf->NumberEnergyBins, rmf->NumberChannels);

#ifdef NORMALIZE_RMF
  // Normalize the RMF:
  headas_printf("### Warning: RMF is explicitly renormalized! ###\n");
  headas_printf(" ARF contributions must be contained in the input spectrum!\n");
  NormalizeRMF(rmf);
#else
  // Check if the RSP file contains matrix rows with a sum of more than 1.
  // In that case the RSP probably also contains the mirror ARF, what should 
  // normally not be the case for this simulation.
  long chancount, bincount;
  double maxsum = 0.;
  for (bincount=0; bincount<rmf->NumberEnergyBins; bincount++) {
    double sum = 0.;
    for (chancount=0; chancount<rmf->NumberChannels; chancount++) {
      sum += ReturnRMFElement(rmf, chancount, bincount);
    }
    if (sum > maxsum) {
      maxsum = sum;
    }
  }
  if (fabs(maxsum-1.)>1.e-3) {
    headas_printf("### Warning: RSP probably contains ARF "
		  "(row-sum = %lf)! ###\n", maxsum);
  }
#endif

  // Read the EBOUNDS extension:
  if ((*status=ReadRMFEbounds(fptr, 0, rmf))!=EXIT_SUCCESS) return(rmf);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  
  return(rmf);
}



void freeRMF(struct RMF* rmf) 
{
  
  if (NULL!=rmf) {
    free(rmf);
    rmf=NULL;
  }
}



long getChannel(float energy, struct RMF* rmf)
{
  // Check if the charge is outside the range of the energy bins defined
  // in the EBOUNDS table. In that case the return value of this function is '-1'.
  if (rmf->ChannelLowEnergy[0] > energy) {
    return(0); // TODO
  } else if (rmf->ChannelHighEnergy[rmf->NumberChannels-1] < energy) {
    return(rmf->NumberChannels - 1 + rmf->FirstChannel);
  }
  
  // Perform a binary search to obtain the detector PHA channel 
  // that corresponds to the given detector charge.
  long min, max, row;
  min = 0;
  max = rmf->NumberChannels-1;
  while (max-min > 1) {
    row = (long)(0.5*(min+max));
    if (rmf->ChannelHighEnergy[row] < energy) {
      min = row;
    } else {
      max = row;
    }
  }
  // Take the final decision wheter max or min is right:
  if (rmf->ChannelLowEnergy[max] < energy) {
    row = max;
  } else {
    row = min;
  }
  
  // Return the PHA channel.
  return(row + rmf->FirstChannel);
}



float getEnergy(long channel, struct RMF* rmf)
{
  // Subtract the channel offset (EBOUNDS may either start at 0 or at 1).
  channel -= rmf->FirstChannel;
  if ((channel < 0) || (channel >= rmf->NumberChannels)) {
    return(-1.);
  }

  // Return the mean of the energy that corresponds to the specified PHA channel
  // according to the EBOUNDS table.
  return(rmf->ChannelLowEnergy[channel] +
	 sixt_get_random_number()*(rmf->ChannelHighEnergy[channel]-
				   rmf->ChannelLowEnergy[channel]));
}



inline double gaussint(double x) 
{
  return(gsl_sf_erf_Q(x));
}
