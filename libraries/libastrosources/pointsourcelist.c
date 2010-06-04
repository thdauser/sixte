#include "pointsourcelist.h"


void freePointSourceList(PointSourceList* psl)
{
  if (NULL!=psl) {
    if (NULL!=psl->sources) {
      long count;
      for(count=0; count<psl->nsources; count++) {
	freeLinLightCurve(psl->sources[count].lc);
      }
      free(psl->sources);
    }
    free(psl);
  }
}



