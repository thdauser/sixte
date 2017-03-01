//
// tessim routines for FITS-based impactlist handling
//

#include "tessim.h"

tes_impactfile_info *tes_init_impactlist(char *impactfile,int *status) {
  CHECK_STATUS_RET(*status,NULL);

  tes_impactfile_info *imp=malloc(sizeof(tes_impactfile_info));
  CHECK_STATUS_RET(*status,NULL);
  imp->impactlist=strdup(impactfile);
  imp->impfile=openPixImpFile(imp->impactlist, READONLY,status);
  CHECK_STATUS_RET(*status,NULL);

  imp->keywords=newSixtStdKeywords(status);
  sixt_read_fits_stdkeywords(imp->impfile->fptr,
			     imp->keywords,
			     status);
  CHECK_STATUS_RET(*status,NULL);

  return(imp);
}

int tes_photon_from_impactlist(PixImpact *impact, void *dataptr, int *status) {
  CHECK_STATUS_RET(*status,0);
  // read a photon from an impact list and return it

  tes_impactfile_info *data=(tes_impactfile_info *) dataptr;

  if (data->impfile==NULL) {
    *status=EXIT_FAILURE;
    return(0);
  }

  // are there still photons in the list?
  if (data->impfile->row >= data->impfile->nrows) {
    // no
    return 0;
  } 

  getNextImpactFromPixImpFile(data->impfile,impact,status);
  CHECK_STATUS_RET(*status,0);

  return(1);
}

void tes_free_impactlist(tes_impactfile_info **data, int *status) {
  free((*data)->impactlist);
  (*data)->impactlist=NULL;
  freePixImpFile(&((*data)->impfile),status);
  freeSixtStdKeywords((*data)->keywords);

  free(*data);
  *data=NULL;
}


tes_impactbuffer_info *tes_init_impactbuffer(int nImpacts, int *status){
  CHECK_STATUS_RET(*status,NULL);
  tes_impactbuffer_info *imp=malloc(sizeof(tes_impactbuffer_info));
  CHECK_STATUS_RET(*status,NULL);

  imp->numimpacts = nImpacts;
  imp->nextimpact = 0;

  imp->impacts = (PixImpact*) malloc(nImpacts*sizeof(PixImpact));
  CHECK_STATUS_RET(*status,NULL);

  // initially set all impacts to E=0 at t=0
  for (int ii=0; ii<nImpacts; ii++) {
    imp->impacts[ii].pixID = 0;
    imp->impacts[ii].time = 0;
    imp->impacts[ii].energy = 0;
    struct Point2d fakepoint = {0., 0.};
    imp->impacts[ii].detposition = fakepoint;
    imp->impacts[ii].pixposition = fakepoint;
    imp->impacts[ii].ph_id = ii;
    imp->impacts[ii].src_id = 1;
  }
  
  CHECK_STATUS_RET(*status,NULL);
  return(imp);
}

int tes_photon_from_impactbuffer(PixImpact *photon, void *dataptr, int *status){
  CHECK_STATUS_RET(*status,0);
  // read a photon from an impact buffer and return it

  tes_impactbuffer_info *data=(tes_impactbuffer_info *) dataptr;

  if (data->impacts==NULL) {
    *status=EXIT_FAILURE;
    return(0);
  }

  // are there still photons in the buffer?
  if (data->nextimpact >= data->numimpacts) {
    // no
    return 0;
  } 
  
  // copy over the contents of the impact
  PixImpact source = data->impacts[data->nextimpact]; // dereference it
  photon->pixID = source.pixID;
  photon->time = source.time;
  photon->energy = source.energy;
  photon->detposition = source.detposition;
  photon->pixposition = source.pixposition;
  photon->ph_id = source.ph_id;
  photon->src_id = source.src_id;

  data->nextimpact++;
  CHECK_STATUS_RET(*status,0);

  return(1);

}

// switch off the warning of status being unused
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void tes_free_impactbuffer(tes_impactbuffer_info **data, int *status){
  free((*data)->impacts);
  (*data)->impacts=NULL;

  free(*data);
  *data=NULL;
}
#pragma GCC diagnostic pop
