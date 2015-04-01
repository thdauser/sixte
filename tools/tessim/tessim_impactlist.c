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
