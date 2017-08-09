//
// output routines for the Tes simulator that include 
// various trigger strategies
//

#include "tessim.h"

#include <assert.h>

// initialize the internal TESDataStream based memory management 
tes_trigger_info *tes_init_trigger(double tstart, double tstop, tesparams *tes, 
				   int strategy,unsigned long preBufferSize,
				   unsigned long triggerSize,
				   double threshold,unsigned int npts,unsigned int suppress,
				   char *streamfile, char *impactfile, int clobber,
				   SixtStdKeywords *keywords, int *status) {

  tes_trigger_info *data=(tes_trigger_info *)malloc(sizeof(tes_trigger_info));
  CHECK_NULL_RET(data,*status,"Memory allocation failed in tes_init_trigger: data structure",NULL);

  // create and initialize fifo. 
  // we initialize to "invalid" - this way we don't have to worry about a buffer that's not yet
  // fully filled
  data->fifo=createTesRecord(triggerSize,tes->delta_t,0,status);
  CHECK_STATUS_RET(*status,NULL);

  for (unsigned long i=0; i<data->fifo->trigger_size; i++) {
    data->fifo->adc_array[i]=0xFFFF;
    data->fifo->adc_double[i]=0.;
  }
  data->fifoind=0;

  // we're not yet buffering
  data->stream=NULL;
  data->streamind=-1;
  data->preBufferSize=preBufferSize;
  data->triggerSize=triggerSize;

  data->tstart=tstart;
  data->tstop=tstop;

  data->imin=tes->imin;
  data->imax=tes->imax;
  data->aducnv=tes->aducnv;

  data->strategy=strategy;
  data->threshold=threshold;
  data->npts=npts;
  data->helper[0]=0xFFFF;
  data->helper[1]=0xFFFF;
  data->helper[2]=0xFFFF;

  data->CanTrigger=0; // allow to trigger
  data->SuppressTrigger=suppress;

  // Initialize Trigger file
  data->fptr=opennewTesTriggerFile(streamfile,
				   keywords,
				   "none", // xmlfile
				   impactfile,
				   triggerSize,
				   preBufferSize,
				   1./tes->delta_t,
				   clobber,
				   status);
  tes_fits_write_params(data->fptr->fptr,tes,status);
  return(data);
}

// coefficients of the backwards differentiation equations
// note: there seems to be no way in C99 to initialize coeff as
// a variable length array -> we waste some memory here but at least
// we can do a compact initialization
//
typedef struct diff_coeff {
  int n;
  int divisor;
  int coeff[11];
} diff_coeff;


// coefficients are the one-sided differentiators of 
// Pavel Holoborodko
// http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
// dy_j/dh = 1/divisor * sum c_i y_{j-i}
diff_coeff diffcoeff[]={
  { .n=0, .divisor=-1, .coeff= {0}}, // dummy entry
  { .n=1, .divisor=-1, .coeff= {0}}, // dummy entry
  { .n=2, .divisor= 2, .coeff= { +1,  0, -1 } },
  { .n=3, .divisor= 4, .coeff= { +1, +1, -1, -1}},
  { .n=4, .divisor= 8, .coeff= { +1, +2,  0, -2, -1}},
  { .n=5, .divisor=16, .coeff= { +1, +3, +2, -2, -3, -1}},
  { .n=6, .divisor=32, .coeff= { +1, +4, +5,  0, -5, -4, -6}},
  { .n=7, .divisor=64, .coeff= { +1, +5, +9, +5, -5, -9, -5, -1}},
  { .n=8, .divisor=128, .coeff= { +1, +6,+14,+14,  0,-14,-14, -6, 1}},
  { .n=9, .divisor=256, .coeff= { +1, +7,+20,+28,+14,-14,-28,-20, -7, -1}},
  { .n=10, .divisor=512, .coeff= { +1, +8,+27,+48,+42,  0,-42,-48,-27, -7, -1}}
};


void tes_append_trigger(tesparams *tes,double time,double pulse, int *status) {
  //
  // trigger based on data->strategy
  //   0: moving average
  //   1: differentiation
  //   2: noise (no triggering, but output records)
  //   3: impact (trigger perfectly, i.e. at every photon impact)
  //
  // The following elements of the helper array are used by the
  // moving average calculation:
  //   helper[0]: sum to calculate the moving average
  //   helper[1]: index into the fifo (for moving average calculation)
  //

  CHECK_STATUS_VOID(*status); 

  tes_trigger_info *data=(tes_trigger_info *) (tes->streaminfo);

  // remember last sampled time
  data->tstop=time;

  // digitize the current signal
  uint16_t pulse16;
  if ((pulse<data->imin) || (pulse>data->imax)) {
    pulse16=0xFFFF;
  } else {
    pulse16=(uint16_t) ((pulse-data->imin)*data->aducnv);
  }

  // save data in the fifo
  data->fifo->adc_double[data->fifoind]=pulse;
  data->fifo->adc_array[data->fifoind]=pulse16;
  data->fifoind++;
  if (data->fifoind==data->fifo->trigger_size) {
    data->fifoind=0;
  }

  // initialize this trigger?
  if (data->helper[0]==0xFFFF) {
    if (data->strategy==TRIGGER_MOVAVG) {
      // initialize moving average sum (with invalid data)
      // (the CanTrigger statement below ensures that the moving
      // average only includes valid data once triggering is allowed
      // to start)
      data->helper[0]=data->npts*0xFFFF;
      data->helper[1]=data->triggerSize-data->npts;
      data->CanTrigger=data->npts; 
    }
  }

  int trigger=0; // initialize to false

  // only check trigger if we are allowed to trigger
  if (data->CanTrigger==0) {
    // only trigger if there's no overflow 
    if (pulse16!=0xFFFF) {
      if (data->strategy==TRIGGER_MOVAVG) {
	// moving average trigger
	// fast trigger: if pulse16/moving average > threshold
	// note: this is prone to integer overflows. 
	// helper[0] is a long, so the probability that this happens
	// is not that large. but we might want to revisit this 
	// and make it more stable
	trigger=(data->npts*pulse16 > data->helper[0]*data->threshold);
      } else {
	// trigger using differentiation
	if (data->strategy==TRIGGER_DIFF) {
	  // calculate filtered derivative
	  // ndx is the index of the newest sample
	  int ndx=data->fifoind-1; 
	  // this loop could be optimized making use of the symmetry
	  // of the coefficients (we might also want to code this in a
	  // different way to prevent integer overflows, although again
	  // I think by using a long while the trigger is an unsigned int 
	  // we have enough buffer
	  assert(diffcoeff[data->npts].n==(int) data->npts);
	  long diffsum=0;
	  for (unsigned int i=0; i<=data->npts; i++) {
	    if (ndx==-1) {
	      ndx=data->fifo->trigger_size-1;
	    }
	    diffsum+=diffcoeff[data->npts].coeff[i]*data->fifo->adc_array[ndx];
	    ndx--;
	  }
	  trigger=(diffsum > diffcoeff[data->npts].divisor*data->threshold);

	  //	  printf("   Derivval: %20.10f %20.10f\n",(double) diffsum/(double) diffcoeff[data->npts].divisor , data->threshold);

	} else {
	  if (data->strategy==TRIGGER_NOISE) {
	    trigger=1;
          } else {
            if (data->strategy==TRIGGER_IMPACT) {
              trigger=(tes->n_absorbed >= 1);
	    } else {
	      fprintf(stderr,"This should never happen!");
            }
	  }
	}
      }
    } else {
      // invalid data: ensure that next trigger does not contain the 0xFFFF
      data->CanTrigger+=data->npts;
    }
  } else {
    data->CanTrigger--;
  }

  // update moving average sum
  // note: we can always do this, including when pulse16==0xFFFF
  // since we prevent triggers if a 0xFFFF is in the queue with
  // the CanTrigger logic
  if (data->strategy==TRIGGER_MOVAVG) {
    data->helper[0]+=(long) pulse16 - (long) data->fifo->adc_array[data->helper[1]];
    // index for last element in moving average sum
    data->helper[1]++;
    if ((unsigned long ) data->helper[1]==data->fifo->trigger_size) {
      data->helper[1]=0;
    }
  }

  if (trigger) {

    //    printf("Trigger: %10.5f\n",time);
    
    // prevent trigger happening for the next SuppressTrigger samples
    data->CanTrigger=data->SuppressTrigger;

    // are we already accumulating data?
    // if yes: increase the record size appropriately
    if (data->stream!=NULL) {
      //      printf("Resizing record to length %ul\n",data->triggerSize+data->streamind);
      resizeTesRecord(data->stream,data->triggerSize+data->streamind,status);
    } else {
      // not accumulating yet: create new tes record
      data->stream=createTesRecord(data->triggerSize,tes->delta_t,0,status);
      data->streamind=0;
      data->stream->pixid=tes->id;
      // copy prebuffer into the new stream
      data->stream->time=time - data->preBufferSize*tes->delta_t;

      // find new starting index. In principle this is
      // just ii=fifoind-1 -preBufferSize but we have to
      // be careful to always remain positive since
      // preBufferSize etc. are unsigned longs!
      unsigned long ii=data->fifoind;
      // wrap around?
      if (ii<=data->preBufferSize) {
	ii+=data->fifo->trigger_size;
      }
      ii=ii-1-data->preBufferSize;

      for (unsigned long i=0; i<data->preBufferSize; i++) {
	if (ii==data->fifo->trigger_size) {
	  ii=0;
	}
	data->stream->adc_array[i]=data->fifo->adc_array[ii];
	data->stream->adc_double[i]=data->fifo->adc_double[ii];
	ii++;
      }
      
      data->streamind=data->preBufferSize;
    }
    CHECK_STATUS_VOID(*status);
  }

  // if we are accumulating data: save data in the stream
  if (data->stream!=NULL) {
    data->stream->adc_double[data->streamind]=pulse;
    data->stream->adc_array[data->streamind]=pulse16;
    data->streamind++;

    // is the buffer filled?
    if (data->streamind==data->stream->trigger_size) {
      // yes: write to file and forget this record
      writeRecord(data->fptr,data->stream,status);
      CHECK_STATUS_VOID(*status);

      freeTesRecord(&(data->stream)); // also sets data->stream to NULL
      data->streamind=-1;

      // TODO leave the option to write this to a buffer of records,
      // so that we can e.g. handle the record internally instead of 
      // writing it to the file
      // strategy: define a buffer struct:
      // If it exists, instead of writing the record, forward the pointer
      // data->stream to the buffer and set data->stream=NULL
      // The buffer struct then has to take care of writing and freeing!
    }
  }
}


/* // save a photon in the current stream */
/* void tes_append_photon_trigger(tesparams *tes,double time, long phid, int *status) { */
/*   CHECK_STATUS_VOID(*status); */
/*   tes_record_info *data=(tes_record_info *) (tes->streaminfo); */
/*   appendPhID(data->stream->phid_list,phid,time,status); */
/* } */


// write remaining stream to file, close fits file
void tes_close_trigger(tesparams *tes,int *status) {
  CHECK_STATUS_VOID(*status);
  tes_trigger_info *streamptr=(tes_trigger_info *) tes->streaminfo;

  // do we still have something in the trigger?
  if (streamptr->stream != NULL) {
    // only write valid data from the current trigger
    unsigned long tmpsize=streamptr->stream->trigger_size;
    streamptr->stream->trigger_size=streamptr->streamind;
    writeRecord(streamptr->fptr,streamptr->stream,status);
    // reset stream length to correct length
    streamptr->stream->trigger_size=tmpsize;
    // forget the current record (sets the pointer to NULL)
    freeTesRecord(&(streamptr->stream));
  }

  freeTesTriggerFile(&(streamptr->fptr),status);
  CHECK_STATUS_VOID(*status);
  streamptr->fptr=NULL;
}

void tes_free_trigger(tes_trigger_info **data,int *status) {
  CHECK_STATUS_VOID(*status);
  if ((*data)->stream != NULL ) {
    freeTesRecord( &((*data)->stream));
  }
  freeTesRecord( &((*data)->fifo));
  free(*data);
  *data=NULL;
}
