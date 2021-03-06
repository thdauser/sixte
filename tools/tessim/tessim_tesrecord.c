/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

//
// output routines for the Tes simulator based on the tesrecord
// data structure
//

#include "tessim.h"

void tes_write_tesrecord(tesparams *tes,int *status);

// initialize the internal TESDataStream based memory management
tes_record_info *tes_init_tesrecord(double tstart, double tstop, tesparams *tes, int buffersize,
				    char *streamfile, char *impactfile, int clobber,int write_error,
				    SixtStdKeywords *keywords, int *status) {
  tes_record_info *data=(tes_record_info *)malloc(sizeof(tes_record_info));
  CHECK_NULL_RET(data,*status,"Memory allocation failed in tes_init_tesrecord: data structure",NULL);

  data->Nt=buffersize;
  data->streamind=-1;
  data->tstart=tstart; // note: tstart/tstop are overwritten in tes_append_tesrecord
  data->tstop=tstop;   // -> need to see whether we really need these parameters here
  data->imin=tes->imin;
  data->imax=tes->imax;
  data->aducnv=tes->aducnv;

  data->stream=newTesRecord(status);
  allocateTesRecord(data->stream,(int) data->Nt, 1./tes->sample_rate, 0, status);
  CHECK_STATUS_RET(*status,NULL);
  data->stream->pixid=tes->id;

  data->streamfile=strdup(streamfile);
  data->impactfile=strdup(impactfile);

  // FITS file is not yet initialized
  data->fptr=NULL;
  data->row=-1;
  data->timecol=-1;
  data->adccol=-1;
  data->curcol=-1;
  if (write_error){
	data->errcol=-1;
  } else {
	data->errcol=-2;
  }
  data->clobber=clobber;

  data->keywords=duplicateSixtStdKeywords(keywords,status);
  if (data->keywords->extname !=NULL) {
    free(data->keywords->extname);
    data->keywords->extname=NULL;
  }
  CHECK_NULL_RET(data,*status,"Failure duplicating keywords",NULL);
  return(data);
}

// write a pulse to the stream
// this routine executes the trigger
void tes_append_tesrecord(tesparams *tes,double time,double pulse, int *status) {
  CHECK_STATUS_VOID(*status);
  tes_record_info *data=(tes_record_info *) (tes->streaminfo);

  if (data->streamind==-1) {
    data->tstart=time;
    data->streamind=0;
  }

  // FIXME: should these to numbers not be of the same type by default then?
  if ((unsigned long) data->streamind==data->stream->trigger_size) {
    tes_write_tesrecord(tes,status);
    data->streamind=0;
    // avoid roundoff
    data->stream->time+=data->stream->delta_t*data->stream->trigger_size;
  }

  // remember how far we've come
  data->tstop=time;

  // save data in the fifo
  data->stream->adc_double[data->streamind]=pulse;
  data->stream->error_double[data->streamind]=tes->squid_error;
  if ((pulse<data->imin) || (pulse>data->imax)) {
    data->stream->adc_array[data->streamind++]=0xFFFF;
  } else {
    data->stream->adc_array[data->streamind++]=(uint16_t) ((pulse-data->imin)*data->aducnv);
  }

}

// save a photon in the current stream
void tes_append_photon_tesrecord(tesparams *tes,double time, long phid, int *status) {
  CHECK_STATUS_VOID(*status);
  tes_record_info *data=(tes_record_info *) (tes->streaminfo);
  appendPhID(data->stream->phid_list,phid,time,status);
}

// write the tesrecord described by argument data to a file
void tes_write_tesrecord(tesparams *tes,int *status) {

  // open FITS file?
  tes_record_info *dataptr=(tes_record_info *) (tes->streaminfo);

  fitsfile *fptr=dataptr->fptr;
  if (dataptr->fptr==NULL) {
    fits_create_file_clobber(&fptr,dataptr->streamfile,dataptr->clobber,status);
    CHECK_STATUS_VOID(*status);
    dataptr->fptr=fptr;

    //
    int logic=(int)'T';
    int bitpix=8;
    int naxis=0;
    fits_update_key(fptr, TLOGICAL, "SIMPLE", &(logic), NULL, status);
    fits_update_key(fptr, TINT, "BITPIX", &(bitpix), NULL, status);
    fits_update_key(fptr, TINT, "NAXIS", &(naxis), NULL, status);
    sixt_add_fits_stdkeywords(fptr,1,dataptr->keywords,status);
    fits_update_key(fptr,TSTRING,"XMLFILE","none",
		    "detector XML description",status);
    fits_update_key(fptr,TSTRING,"PIXFILE",dataptr->impactfile,
		    "Pixel impact file for this stream",status);
    CHECK_STATUS_VOID(*status);

    //pixel number
    long pixid=dataptr->stream->pixid;

    // create table extension
    int ncolumns;
    if (dataptr->errcol==-2){
      ncolumns=3;
    } else {
      ncolumns=4;
    }

    char *ttype[ncolumns];
    char *tform[ncolumns];
    char *tunit[ncolumns];

    // column 0: time
    ttype[0]="TIME";
    tform[0]="1D";
    tunit[0]="s";
    dataptr->timecol=1; // Column number, not array!

    // pixel
    ttype[1]=malloc(9*sizeof(char));
    CHECK_NULL_VOID(ttype[1],*status,"Cannot allocate ttype[1]");
    sprintf(ttype[1],"PXL%05ld",pixid);
    tform[1]="1U";
    tunit[1]="ADC";
    dataptr->adccol=2;

    // current
    ttype[2]=malloc(11*sizeof(char));
    CHECK_NULL_VOID(ttype[2],*status,"Cannot allocate ttype[2]");
    sprintf(ttype[2],"PULSE%05ld",pixid);
    tform[2]="1E";
    tunit[2]="A";
    dataptr->curcol=3;

    // error
    if (dataptr->errcol==-1){
      ttype[3]=malloc(11*sizeof(char));
      CHECK_NULL_VOID(ttype[3],*status,"Cannot allocate ttype[3]");
      sprintf(ttype[3],"ERROR%05ld",pixid);
      tform[3]="1E";
      tunit[3]="A";
      dataptr->errcol=4;
    }


    fits_create_tbl(fptr,BINARY_TBL,0,ncolumns,
		    ttype,tform,tunit,"TESDATASTREAM",status);
    CHECK_STATUS_VOID(*status);

    dataptr->row=1;

    // Write header keywords
    fits_update_key(fptr, TDOUBLE, "TSTART",
		    &(dataptr->tstart), "Start time of data stream", status);
    fits_update_key(fptr, TDOUBLE, "TSTOP",
		    &(dataptr->tstop), "Stop time of data stream", status);
    fits_update_key(fptr, TDOUBLE, "DELTAT",
		    &(dataptr->stream->delta_t), "Time resolution of data stream", status);
    int npix=1;
    fits_update_key(fptr, TINT, "NPIX",
		    &npix, "Number of pixel streams in extension", status);
    fits_update_key(fptr, TLONG, "FIRSTPIX",
		    &pixid, "ID of first pixel in extension", status);
    fits_update_key(fptr, TLONG, "LASTPIX",
		    &pixid, "ID of last pixel in extension", status);
    CHECK_STATUS_VOID(*status);

    tes_fits_write_params(fptr,tes,status);

    free(ttype[1]);
    free(ttype[2]);
  }

  LONGLONG nrows=dataptr->streamind; //this works because FITS is 1 based, arrays are 0 based

  // write the data
  double *time=malloc(nrows*sizeof(double));
  CHECK_NULL_VOID(time,*status,"Cannot allocate temporary time array\n");
  for (LONGLONG i=0; i<nrows;i++) {
    time[i]=dataptr->stream->time + i*dataptr->stream->delta_t;
  }
  fits_write_col(fptr, TDOUBLE, dataptr->timecol, dataptr->row,1, nrows, time, status);
  free(time);
  fits_write_col(fptr, TUSHORT, dataptr->adccol, dataptr->row, 1, nrows, dataptr->stream->adc_array, status);
  //NB this will be written as single precision!
  fits_write_col(fptr, TDOUBLE,  dataptr->curcol, dataptr->row, 1, nrows, dataptr->stream->adc_double, status);
  if (dataptr->errcol>0){
    fits_write_col(fptr, TDOUBLE,  dataptr->errcol, dataptr->row, 1, nrows, dataptr->stream->error_double,status);
  }
  dataptr->row+=nrows;

}

// write remaining tesrecord to file, close fits file
void tes_close_tesrecord(tesparams *tes,int *status) {
  CHECK_STATUS_VOID(*status);
  tes_record_info *streamptr=(tes_record_info *) tes->streaminfo;
  tes_write_tesrecord(tes,status);
  fits_close_file_chksum(streamptr->fptr,status);
  CHECK_STATUS_VOID(*status);
  streamptr->fptr=NULL;
}

void tes_free_tesrecord(tes_record_info **data,int *status) {
  CHECK_STATUS_VOID(*status);
  freeTesRecord( &((*data)->stream));
  free((*data)->streamfile);
  free((*data)->impactfile);
  free(*data);
  *data=NULL;
}
