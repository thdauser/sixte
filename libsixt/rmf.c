#include "rmf.h"


struct RMF* loadRMF(char* filename, int* const status) 
{
  struct RMF* rmf=(struct RMF*)malloc(sizeof(struct RMF));
  CHECK_NULL_RET(rmf, *status, "memory allocation for RMF failed", rmf);

  // Load the RMF from the FITS file using the HEAdas access routines
  // (part of libhdsp).
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_RET(*status, rmf);
  
  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  *status=ReadRMFMatrix(fptr, 0, rmf);
  CHECK_STATUS_RET(*status, rmf);

  // Print some information:
  headas_chat(5, "RMF loaded with %ld energy bins and %ld channels\n",
	      rmf->NumberEnergyBins, rmf->NumberChannels);

  // Check if the RMF file contains matrix rows with a sum of more than 1.
  // In that case the RSP probably also contains the mirror ARF, what should 
  // not be the case for this simulation. Row sums with a value of less than
  // 1 should actually also not be used, but can be handled by the simulation.
  long bincount;
  double min_sum=1.;
  for (bincount=0; bincount<rmf->NumberEnergyBins; bincount++) {
    double sum=0.;
    long chancount;
    for (chancount=0; chancount<rmf->NumberChannels; chancount++) {
      sum+=ReturnRMFElement(rmf, chancount, bincount);
    }
    if (sum>1.000001) {
      SIXT_ERROR("RMF contains rows with a sum > 1.0 (probably contains ARF)");
      *status=EXIT_FAILURE;
      return(rmf);
    }
    if (sum<min_sum) {
      min_sum=sum;
    }
  }

  if (min_sum<0.999999) {
    SIXT_WARNING("RMF is not normalized");
  }


  // Read the EBOUNDS extension:
  *status=ReadRMFEbounds(fptr, 0, rmf);
  CHECK_STATUS_RET(*status, rmf);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, rmf);

  return(rmf);
}


void loadArfRmfFromRsp(char* filename, 
		       struct ARF** arf,
		       struct RMF** rmf,
		       int* const status) 
{
  *rmf=(struct RMF*)malloc(sizeof(struct RMF));
  CHECK_NULL_VOID(*rmf, *status, "memory allocation for RSP failed");

  // Load the RSP from the FITS file using the HEAdas access routines
  // (part of libhdsp).
  fitsfile* fptr=NULL;
  fits_open_file(&fptr, filename, READONLY, status);
  CHECK_STATUS_VOID(*status);
  
  // Read the 'SPECRESP MATRIX' or 'MATRIX' extension:
  *status=ReadRMFMatrix(fptr, 0, *rmf);
  CHECK_STATUS_VOID(*status);

  // Print some information:
  headas_chat(5, "RSP loaded with %ld energy bins and %ld channels\n",
	      (*rmf)->NumberEnergyBins, (*rmf)->NumberChannels);

  // Allocate memory for the ARF.
  *arf=getARF(status);
  CHECK_STATUS_VOID(*status);

  // Produce an ARF from the RSP data.
  (*arf)->NumberEnergyBins=(*rmf)->NumberEnergyBins;
  (*arf)->LowEnergy=(float*)malloc((*arf)->NumberEnergyBins*sizeof(float));
  CHECK_NULL_VOID((*arf)->LowEnergy, *status, 
		  "memory allocation for energy bins failed");
  (*arf)->HighEnergy=(float*)malloc((*arf)->NumberEnergyBins*sizeof(float));
  CHECK_NULL_VOID((*arf)->HighEnergy, *status, 
		  "memory allocation for energy bins failed");
  (*arf)->EffArea=(float*)malloc((*arf)->NumberEnergyBins*sizeof(float));
  CHECK_NULL_VOID((*arf)->EffArea, *status, 
		  "memory allocation for effective area failed");
  strcpy((*arf)->Telescope, (*rmf)->Telescope);
  strcpy((*arf)->Instrument, (*rmf)->Instrument);
  strcpy((*arf)->Detector, (*rmf)->Detector);
  strcpy((*arf)->Filter, (*rmf)->Filter);

  // Calculate the row sums.
  long bincount;
  int notrmf=0;
  for (bincount=0; bincount<(*rmf)->NumberEnergyBins; bincount++) {
    double sum=0.;
    long chancount;
    for (chancount=0; chancount<(*rmf)->NumberChannels; chancount++) {
      sum+=ReturnRMFElement(*rmf, chancount, bincount);
    }
    
    // Check if this might be an RMF and not an RSP, as it should be.
    if ((0==notrmf)&&(fabs(sum-1.0)<0.001)) {
      SIXT_WARNING("response matrix declared as RSP file looks like RMF");
    } else {
      notrmf=1;
    }
    
    (*arf)->LowEnergy[bincount]=(*rmf)->LowEnergy[bincount];
    (*arf)->HighEnergy[bincount]=(*rmf)->HighEnergy[bincount];
    (*arf)->EffArea[bincount]=(float)sum;
  }

  // Normalize the RSP to obtain an ARF.
  NormalizeRMF(*rmf);

  // Read the EBOUNDS extension:
  *status=ReadRMFEbounds(fptr, 0, *rmf);
  CHECK_STATUS_VOID(*status);

  // Close the open FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


void freeRMF(struct RMF* const rmf) 
{
  if (NULL!=rmf) {
    free(rmf);
  }
}


long getEBOUNDSChannel(const float energy, const struct RMF* const rmf)
{
  // In case there is no RMF, just return a negativ number (-1).
  if (NULL==rmf) return(-1);

  // Check if the charge is outside the range of the energy bins defined
  // in the EBOUNDS table. In that case the return value of this function 
  // is '-1'.
  if (rmf->ChannelLowEnergy[0] > energy) {
    return(-1);
  } else if (rmf->ChannelHighEnergy[rmf->NumberChannels-1] < energy) {
    return(-1);
    //return(rmf->NumberChannels - 1 + rmf->FirstChannel);
  }
  
  // Perform a binary search to obtain the detector PHA channel 
  // that corresponds to the given detector charge.
  long min, max, mid;
  min = 0;
  max = rmf->NumberChannels-1;
  while (max > min) {
    mid = (min+max)/2;
    if (rmf->ChannelHighEnergy[mid] < energy) {
      min = mid+1;
    } else {
      max = mid;
    }
  }
  
  // Return the PHA channel.
  return(min + rmf->FirstChannel);
}


float getEBOUNDSEnergy(long channel, 
		       const struct RMF* const rmf, 
		       const int boundary,
		       int* const status)
{
  // Subtract the channel offset (EBOUNDS may either start at 0 or at 1).
  channel -= rmf->FirstChannel;
  if ((channel < 0) || (channel >= rmf->NumberChannels)) {
    return(-1.);
  }

  // Return the randomized mean / lower / upper energy value corresponding 
  // to the specified PHA channel according to the EBOUNDS table.
  if (0==boundary) {
    return(rmf->ChannelLowEnergy[channel] +
	   sixt_get_random_number(status)*
	   (rmf->ChannelHighEnergy[channel]-rmf->ChannelLowEnergy[channel]));
  } else if (-1==boundary) {
    return(rmf->ChannelLowEnergy[channel]);
  } else {
    return(rmf->ChannelHighEnergy[channel]);
  }
}


void returnRMFChannel(struct RMF *rmf, 
		      const float energy, 
		      long* const channel)
{
  ReturnChannel(rmf, energy, 1, channel);
}
