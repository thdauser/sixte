//
// TESFitsStream based writing routines for the TES simulator
//

#include "tessim.h"

// initialize the internal TESDataStream based memory management 
tes_datastream_info *tes_init_datastream(double tstart, double tstop,tesparams *tes, int *status) {
  tes_datastream_info *data=(tes_datastream_info *)malloc(sizeof(tes_datastream_info));
  CHECK_NULL_RET(data,*status,"Memory allocation failed for TES Datastream structure",NULL);

  data->Nt=(unsigned long) ((tstop-tstart)*tes->sample_rate);
  data->streamind=0;
  data->tstart=0;
  data->tstop=0;
  data->imin=tes->imin;
  data->imax=tes->imax;
  data->aducnv=tes->aducnv;
  data->stream=newTESDataStream(status);
  if (*status==EXIT_FAILURE) {
    SIXT_ERROR("memory allocation failed for stream structure");
    CHECK_STATUS_RET(*status,NULL);
  }
  allocateTESDataStream(data->stream,data->Nt,1,status);
  CHECK_STATUS_RET(*status,NULL);
  return(data);
}

void tes_append_datastream(double time,double pulse,void *dataptr, int *status) {
  CHECK_STATUS_VOID(*status); 

  tes_datastream_info *data=(tes_datastream_info *) dataptr;

  if (data->streamind==0) {
    data->tstart=time;
  }

  // only write sample if we still have space in the stream
  if (data->streamind < data->Nt) {
    data->stream->time[data->streamind]=time;
    data->tstop=time; 
    if ((pulse<data->imin) || (pulse>data->imax)) {
      data->stream->adc_value[data->streamind++][0]=0xFFFF;
    } else {
	data->stream->adc_value[data->streamind++][0]=(uint16_t) ((pulse-data->imin)*data->aducnv);
    }
  }
}

// write the datastream described by argument data to a file
void tes_save_datastream(char *streamfile, char *impactfile,
			 tes_datastream_info *data, tesparams *tes, 
			 SixtStdKeywords *keywords, int *status) {
  fitsfile *fptr;
  createTESFitsStreamFile(&fptr,
			  streamfile,
			  keywords->telescop,
			  keywords->instrume,
			  keywords->filter,
			  keywords->ancrfile,
			  keywords->respfile,
			  "none", // xmlfile
			  impactfile, // name of impactfile
			  keywords->mjdref,
			  keywords->timezero,
			  data->tstart,
			  data->tstop,
			  'y', //clobber
			  status);

  CHECK_STATUS_VOID(*status);
  
  // add properties of simulation to header
  tes_fits_write_params(fptr,tes,status);

  TESFitsStream *stream=newTESFitsStream(status);
  snprintf(stream->name,9,"ADC%03d",1);
  allocateTESFitsStream(stream,data->stream->Ntime,1,status);
  // copy over (this is stupid for a single pixel)
  unsigned long ii;
  stream->pixID[0]=0;
  for (ii=0; ii<stream->Ntime; ii++) {
    stream->time[ii]=data->stream->time[ii];
    stream->adc_value[0][ii]=data->stream->adc_value[ii][0];
  }
  writeTESFitsStream(fptr,
		     stream,
		     data->tstart,
		     data->tstop,
		     tes->timeres,
		     &(tes->Nevts), 
		     0,-1,status);
  fits_close_file(fptr,status);
  CHECK_STATUS_VOID(*status);

  destroyTESFitsStream(stream);
  free(stream);

}

void tes_free_datastream(tes_datastream_info **data,int *status) {
  CHECK_STATUS_VOID(*status);
  destroyTESDataStream( (*data)->stream);
  free((*data)->stream);
  free((*data));
  *data=NULL;
}
