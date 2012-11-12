#include "ladimpact.h"


LADImpact* getLADImpact(int* const status)
{
  LADImpact* imp = (LADImpact*)malloc(sizeof(LADImpact));
  CHECK_NULL_RET(imp, *status, 
		 "memory allocation for LADImpact failed", imp);
  
  // Initalize.
  imp->panel  =0;
  imp->module =0;
  imp->element=0;
  imp->position.x=0.;
  imp->position.y=0.;
  imp->energy =0.;
  imp->time   =0.;
  imp->ph_id  =0;
  imp->src_id =0;

  return(imp);
}


void freeLADImpact(LADImpact** const impact)
{
  if (NULL!=*impact) {
    free(*impact);
    *impact=NULL;
  }
}

