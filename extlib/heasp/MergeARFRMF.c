/* Routine to multiply an ARF into an RMF */

#include <string.h>
#include <stdio.h>
#include "heasp.h"


long MergeARFRMF(struct ARF *arf, struct RMF *rmf)
{
  long i, j, k;
  float effarea;

  /* First check that the ARF and RMF energy binnings are compatible.
     At a later time I could generalize this to allow interpolation */

  if (arf->NumberEnergyBins != rmf->NumberEnergyBins) return(-1);

  for (i=0; i<arf->NumberEnergyBins; i++) {
    if (arf->LowEnergy[i] != rmf->LowEnergy[i]) return(i);
    if (arf->HighEnergy[i] != rmf->HighEnergy[i]) return(i);
  }

  /* loop round energy bins multiplying appropriate elements of the RMF by
     the effective area for this energy from the ARF */

  for (i=0; i<arf->NumberEnergyBins; i++) {

    effarea = arf->EffArea[0];

    /* loop round response groups for this energy */

    for (j=rmf->FirstGroup[i]; j<rmf->FirstGroup[i]+rmf->NumberGroups[i]; j++) {

      /* loop round matrix elements for this response group */

      for (k=rmf->FirstChannelGroup[j]; 
	   k<rmf->FirstChannelGroup[j]+rmf->NumberChannelGroups[j]; k++) {

	rmf->Matrix[k] *= effarea;

      }

    }

  }
	
  return(0);

}
