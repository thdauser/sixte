#include "erodetbkgrndgen.h"

static eroBackgroundInput bkginputdata;
static eroBackgroundRateFct bkgratefct;
static struct rateCurrentInterval rCurr = {0L, 0., NULL, 0L};

/* internal functions */

int calcEvents(const double* const hit_time, const long numrows) {
  int cc = 0;
  int numevents = 1;

  for(cc = 0; cc < (numrows - 1); cc++) {
    if(hit_time[cc] != hit_time[cc + 1]) {
      numevents++;
    }
  }

  return numevents;
}

void fillEventList(const double* const hit_time, const long numrows, size_t* const eventlist) {
  int cc = 0, dd = 0;

  eventlist[dd] = cc;
  for(cc = 0; cc < (numrows - 1); cc++) {
    if(hit_time[cc] != hit_time[cc + 1]) {
      dd++;
      eventlist[dd] = cc + 1;
    }
  }
}

inline double calcEventRate(const double* const hit_time,
                                   const long numrows,
                                   const int numevents,
                                   const double interval) {
  double eventrate = (double)numevents;

  eventrate /= (hit_time[numrows - 1] - hit_time[0]);
  eventrate *= interval;

  return eventrate;
}

/* ------------------ */

/* external functions */

void eroBkgInitialize(const char* const filename,
                      int* const status) {

  bkgratefct.numelements = 0L;
  bkgratefct.time = NULL;
  bkgratefct.rate = NULL;
  bkgratefct.currenttime = NULL;
  bkgratefct.currentrate = NULL;
  bkgratefct.intervalsum = 0.;

  bkginputdata.eventlist = NULL;
  bkginputdata.numrows = 0L;
  bkginputdata.timecolname = "TIME CCD";
  bkginputdata.energycolname = "PI CCD";
  bkginputdata.xcolname = "X";
  bkginputdata.ycolname = "Y";
  bkginputdata.interval = 0.;
  ftime(&time_struct);
  srand((unsigned int)time_struct.millitm - (unsigned int)time_struct.time);

  fits_open_table(&bkginputdata.inputfptr, filename, READONLY, status);
  fits_report_error(stderr, *status);

  fits_get_num_rows(bkginputdata.inputfptr, &bkginputdata.numrows, status);
  fits_report_error(stderr, *status);

  bkginputdata.hit_xpos = (double*) malloc(bkginputdata.numrows * sizeof(double));
  bkginputdata.hit_ypos = (double*) malloc(bkginputdata.numrows * sizeof(double));
  bkginputdata.hit_time = (double*) malloc(bkginputdata.numrows * sizeof(double));
  bkginputdata.hit_energy = (double*) malloc(bkginputdata.numrows * sizeof(double));

  fits_get_colnum(bkginputdata.inputfptr, CASEINSEN, bkginputdata.timecolname, &bkginputdata.timecolnum, status);
  fits_get_colnum(bkginputdata.inputfptr, CASEINSEN, bkginputdata.xcolname, &bkginputdata.xcolnum, status);
  fits_get_colnum(bkginputdata.inputfptr, CASEINSEN, bkginputdata.ycolname, &bkginputdata.ycolnum, status);
  fits_get_colnum(bkginputdata.inputfptr, CASEINSEN, bkginputdata.energycolname, &bkginputdata.energycolnum, status);
  fits_report_error(stderr, *status);
  
  fits_read_col(bkginputdata.inputfptr, TDOUBLE, bkginputdata.timecolnum, 1L, 1L, bkginputdata.numrows, 0, bkginputdata.hit_time, NULL, status);
  fits_read_col(bkginputdata.inputfptr, TDOUBLE, bkginputdata.xcolnum, 1L, 1L, bkginputdata.numrows, 0, bkginputdata.hit_xpos, NULL, status);
  fits_read_col(bkginputdata.inputfptr, TDOUBLE, bkginputdata.ycolnum, 1L, 1L, bkginputdata.numrows, 0, bkginputdata.hit_ypos, NULL, status);
  fits_read_col(bkginputdata.inputfptr, TDOUBLE, bkginputdata.energycolnum, 1L, 1L, bkginputdata.numrows, 0, bkginputdata.hit_energy, NULL, status);
  fits_report_error(stderr, *status);

  bkginputdata.numevents = calcEvents(bkginputdata.hit_time, bkginputdata.numrows);
  bkginputdata.eventlist = (size_t*) calloc(bkginputdata.numevents, sizeof(size_t));
  fillEventList(bkginputdata.hit_time, bkginputdata.numrows, bkginputdata.eventlist);

  bkginputdata.eventsperinterval = calcEventRate(bkginputdata.hit_time,
                                                 bkginputdata.numrows,
                                                 bkginputdata.numevents,
                                                 bkginputdata.interval);

  bkginputdata.randgen = gsl_rng_alloc(gsl_rng_ranlux);
  gsl_rng_set(bkginputdata.randgen, rand() * ((int)time_struct.time + (unsigned int)time_struct.millitm));
}

void eroBkgFree(eroBackgroundOutput* struct_to_free) {
	if(struct_to_free->numhits != 0) {
		free(struct_to_free->hit_xpos);
		free(struct_to_free->hit_ypos);
		free(struct_to_free->hit_time);
		free(struct_to_free->hit_energy);
	}
  free(struct_to_free);
}

void eroBkgCleanUp(int* const status) {
	fits_close_file(bkginputdata.inputfptr, status);
	fits_report_error(stderr, *status);
	bkginputdata.inputfptr = NULL;

	gsl_rng_free(bkginputdata.randgen);

	free(bkginputdata.hit_xpos);
	free(bkginputdata.hit_ypos);
	free(bkginputdata.hit_time);
	free(bkginputdata.hit_energy);
}

void eroBkgSetRateFct(const char* const filename, int* const status) {
  SimputLC* rate_lc=NULL;

  if(filename != NULL) {
    rate_lc = loadSimputLC(filename, status);
    bkgratefct.numelements = rate_lc->nentries;
    if(bkgratefct.time != NULL) {
      free(bkgratefct.time);
    }
    bkgratefct.time = (double*) malloc(rate_lc->nentries * sizeof(double));
    memcpy(bkgratefct.time, rate_lc->time, rate_lc->nentries * sizeof(double));

    if(bkgratefct.rate != NULL) {
      free(bkgratefct.rate);
    }
    bkgratefct.rate = (float*) malloc(rate_lc->nentries * sizeof(float));
    memcpy(bkgratefct.rate, rate_lc->flux, rate_lc->nentries * sizeof(float));

    bkgratefct.currenttime = bkgratefct.time;
    bkgratefct.currentrate = bkgratefct.rate;

    freeSimputLC(&rate_lc);
  } else {
    SIXT_ERROR("no simput filename specified for rate function!");
  }
}

void eroBkgGetRate(double interval) {
  double startfraction = 0.;
  const double startinterval = interval;

  /** if we are not at EOF: check if we start in the middle of a time slice */
  if((bkgratefct.intervalsum > 0) && (bkgratefct.currenttime <= &bkgratefct.time[bkgratefct.numelements - 1])) {
    startfraction = *(bkgratefct.currenttime + 1) - (bkgratefct.time[0] + bkgratefct.intervalsum);
  }
  bkgratefct.intervalsum += interval;

  /** free rate structure if it contains data from recent calculations */
  if(rCurr.rate != NULL) {
    free(rCurr.rate);
    rCurr.rate = NULL;
    rCurr.ratesize = 0L;
    rCurr.numelements = 0L;
  }

  /** if we start in the middle of a time slice, we have to take a certain fraction of the rate into account */
  if(startfraction > 0) {
    rCurr.rate = (float*) realloc(rCurr.rate, 10 * sizeof(float));
    rCurr.rate[0] = *bkgratefct.currentrate;
    if(startfraction < startinterval) {
      rCurr.rate[0] *= startfraction / startinterval;
    }
    rCurr.numelements++;
    interval -= startfraction;
  }

  /** if we are still not at EOF, begin looping over the rate entries until the requested interval has been covered */
  while((interval > 0) && (bkgratefct.currentrate <= &bkgratefct.rate[bkgratefct.numelements - 1])) {
    if(rCurr.numelements >= rCurr.ratesize) {
      rCurr.rate = (float*) realloc(rCurr.rate, (rCurr.numelements + 50) * sizeof(float));
      rCurr.ratesize = rCurr.numelements + 50;
    }
    interval -= *(bkgratefct.currenttime + 1) - *bkgratefct.currenttime;
    if(interval < 0) {
      rCurr.rate[rCurr.numelements] = *bkgratefct.currentrate;
      rCurr.rate[rCurr.numelements] *= (interval + (*(bkgratefct.currenttime + 1) - *bkgratefct.currenttime)) / startinterval;
    } else {
      rCurr.rate[rCurr.numelements] = *bkgratefct.currentrate;
      rCurr.rate[rCurr.numelements] *= (*(bkgratefct.currenttime + 1) - *bkgratefct.currenttime) / startinterval;
      bkgratefct.currenttime++;
      bkgratefct.currentrate++;
    }
    rCurr.numelements++;
  }

  if(rCurr.numelements == 0) {
    rCurr.numelements = 1L;
    rCurr.rate = (float*) realloc(rCurr.rate, rCurr.numelements * sizeof(float));
    if((interval != 0) && (bkgratefct.currentrate <= &bkgratefct.rate[bkgratefct.numelements - 1])) {
      rCurr.rate[0] = *bkgratefct.currentrate;
    } else {
      rCurr.rate[0] = 1;
      rCurr.ratesize = rCurr.numelements;
    }
  } else {
    rCurr.rate = (float*) realloc(rCurr.rate, rCurr.numelements * sizeof(float));
    rCurr.ratesize = rCurr.numelements;
  }
}

eroBackgroundOutput* eroBkgGetBackgroundList(double interval) {
  int cc = 0;
	eroBackgroundOutput* bkgresultlist = NULL;
  bkgresultlist = (eroBackgroundOutput*) realloc(bkgresultlist, sizeof(eroBackgroundOutput));
  bkgresultlist->numevents = 0;
  bkgresultlist->numhits = 0;

  if(interval > 0) {
    bkginputdata.interval = interval;
    bkginputdata.eventsperinterval = calcEventRate(bkginputdata.hit_time,
                                                  bkginputdata.numrows,
                                                  bkginputdata.numevents,
                                                  bkginputdata.interval);
  } else {
    SIXT_ERROR("invalid interval for background generation specified!");
    free(bkgresultlist);
    bkgresultlist = NULL;
    return bkgresultlist;
  }
  
  if(bkgratefct.numelements == 0) {
    bkgresultlist->numevents = gsl_ran_poisson(bkginputdata.randgen, bkginputdata.eventsperinterval);
  } else {
    eroBkgGetRate(interval);
    for(cc = 0; cc < rCurr.numelements; cc++) {
      bkgresultlist->numevents += gsl_ran_poisson(bkginputdata.randgen, bkginputdata.eventsperinterval * rCurr.rate[cc]);
    }
  }

  bkgresultlist->numhits = 0;
  if(bkgresultlist->numevents > 0) {
    int rand = 0, hitcnt = 0, arrsize = bkgresultlist->numevents;
    bkgresultlist->hit_energy = (double*) malloc(arrsize * sizeof(double));
    bkgresultlist->hit_time = (double*) malloc(arrsize * sizeof(double));
    bkgresultlist->hit_xpos = (double*) malloc(arrsize * sizeof(double));
    bkgresultlist->hit_ypos = (double*) malloc(arrsize * sizeof(double));

    for(cc = 0; cc < bkgresultlist->numevents; cc++) {
      rand = (int)floor((bkginputdata.numevents - 1) * gsl_ran_flat(bkginputdata.randgen, 0, 1));
      do {
        bkgresultlist->hit_energy[bkgresultlist->numhits] = bkginputdata.hit_energy[bkginputdata.eventlist[rand] + hitcnt];
        bkgresultlist->hit_time[bkgresultlist->numhits] = bkginputdata.hit_time[bkginputdata.eventlist[rand] + hitcnt];
        bkgresultlist->hit_xpos[bkgresultlist->numhits] = bkginputdata.hit_xpos[bkginputdata.eventlist[rand] + hitcnt];
        bkgresultlist->hit_ypos[bkgresultlist->numhits] = bkginputdata.hit_ypos[bkginputdata.eventlist[rand] + hitcnt];
        bkgresultlist->numhits++;
        if(&bkginputdata.hit_time[bkginputdata.eventlist[rand] + hitcnt] != &bkginputdata.hit_time[bkginputdata.numrows - 1]) {
          if(bkginputdata.hit_time[bkginputdata.eventlist[rand] + hitcnt + 1] != bkginputdata.hit_time[bkginputdata.eventlist[rand] + hitcnt]) {
            hitcnt = 0;
          } else {
            hitcnt++;
          }
        } else {
          hitcnt = 0;
        }
        if(bkgresultlist->numhits == arrsize) {
          arrsize += 50;
          bkgresultlist->hit_energy = (double*) realloc(bkgresultlist->hit_energy, arrsize * sizeof(double));
          bkgresultlist->hit_time = (double*) realloc(bkgresultlist->hit_time, arrsize * sizeof(double));
          bkgresultlist->hit_xpos = (double*) realloc(bkgresultlist->hit_xpos, arrsize * sizeof(double));
          bkgresultlist->hit_ypos = (double*) realloc(bkgresultlist->hit_ypos, arrsize * sizeof(double));
        }
      } while(hitcnt > 0);
    }

    if(bkgresultlist->numhits < arrsize) {
    	if(bkgresultlist->numhits == 0) {
    		arrsize = 0;
    		free(bkgresultlist->hit_energy);
    		free(bkgresultlist->hit_time);
    		free(bkgresultlist->hit_xpos);
    		free(bkgresultlist->hit_ypos);
    	} else {
				arrsize = bkgresultlist->numhits;
				bkgresultlist->hit_energy = (double*) realloc(bkgresultlist->hit_energy, arrsize * sizeof(double));
				bkgresultlist->hit_time = (double*) realloc(bkgresultlist->hit_time, arrsize * sizeof(double));
				bkgresultlist->hit_xpos = (double*) realloc(bkgresultlist->hit_xpos, arrsize * sizeof(double));
				bkgresultlist->hit_ypos = (double*) realloc(bkgresultlist->hit_ypos, arrsize * sizeof(double));
    	}
    }
  }
  
  return bkgresultlist;
}

/* ------------------ */
