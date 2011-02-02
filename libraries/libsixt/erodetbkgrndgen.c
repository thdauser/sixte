#include "erodetbkgrndgen.h"

static eroBackgroundInput bkginputdata;

/* internal functions */

int calcEvents(double *hit_time, long numrows) {
  int cc = 0;
  int numevents = 0;

  for(cc = 0; cc < (numrows - 1); cc++) {
    if(hit_time[cc] != hit_time[cc + 1]) {
      numevents++;
    }
  }

  return numevents;
}

inline double calcEventRate(double *hit_time,
                            long numrows,
                            int numevents,
                            double interval) {
  double eventrate = (double)numevents;

  eventrate /= hit_time[numrows - 1];
  eventrate *= interval;

  return eventrate;
}

/* ------------------ */

/* external functions */

void eroBkgInitialize(const char *filename,
                      int *status) {

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
                                       
  bkginputdata.eventsperinterval = calcEventRate(bkginputdata.hit_time,
                                                  bkginputdata.numrows,
                                                  bkginputdata.numevents,
                                                  bkginputdata.interval);

  bkginputdata.randgen = gsl_rng_alloc(gsl_rng_ranlux);
  gsl_rng_set(bkginputdata.randgen, rand() * ((int)time_struct.time + (unsigned int)time_struct.millitm));
}

void eroBkgFree(eroBackgroundOutput *struct_to_free) {
	if(struct_to_free->numhits != 0) {
		free(struct_to_free->hit_xpos);
		free(struct_to_free->hit_ypos);
		free(struct_to_free->hit_time);
		free(struct_to_free->hit_energy);
	}
  free(struct_to_free);
}

void eroBkgCleanUp(int *status) {
	fits_close_file(bkginputdata.inputfptr, status);
	fits_report_error(stderr, *status);
	bkginputdata.inputfptr = NULL;

	gsl_rng_free(bkginputdata.randgen);

	free(bkginputdata.hit_xpos);
	free(bkginputdata.hit_ypos);
	free(bkginputdata.hit_time);
	free(bkginputdata.hit_energy);
}

eroBackgroundOutput *eroBkgGetBackgroundList(double interval) {
	eroBackgroundOutput *bkgresultlist = NULL;
  bkgresultlist = (eroBackgroundOutput*) realloc(bkgresultlist, sizeof(eroBackgroundOutput));

  if(interval > 0) {
    bkginputdata.interval = interval;
    bkginputdata.eventsperinterval = calcEventRate(bkginputdata.hit_time,
                                                  bkginputdata.numrows,
                                                  bkginputdata.numevents,
                                                  bkginputdata.interval);
  } else {
    printf("background generation error: invalid interval specified: %f\nreturning NULL pointer...\n", interval);
    free(bkgresultlist);
    bkgresultlist = NULL;
    return bkgresultlist;
  }
  
  bkgresultlist->numevents = gsl_ran_poisson(bkginputdata.randgen, bkginputdata.eventsperinterval);
  bkgresultlist->numhits = 0;
  if(bkgresultlist->numevents > 0) {
    int cc = 0, rand = 0, arrsize = bkgresultlist->numevents;
    bkgresultlist->hit_energy = (double*) malloc(arrsize * sizeof(double));
    bkgresultlist->hit_time = (double*) malloc(arrsize * sizeof(double));
    bkgresultlist->hit_xpos = (double*) malloc(arrsize * sizeof(double));
    bkgresultlist->hit_ypos = (double*) malloc(arrsize * sizeof(double));

    for(cc = 0; cc < bkgresultlist->numevents; cc++) {
      rand = (int)floor((bkginputdata.numrows - 1) * gsl_ran_flat(bkginputdata.randgen, 0, 1));
      while((rand < (bkginputdata.numrows - 1)) && (bkginputdata.hit_time[rand] == bkginputdata.hit_time[rand + 1])) {
        rand++;
      }
      rand++;
      if(rand >= (bkginputdata.numrows - 1)) {
        rand = 0;
      }
      do {
        bkgresultlist->hit_energy[bkgresultlist->numhits] = bkginputdata.hit_energy[rand];
        bkgresultlist->hit_time[bkgresultlist->numhits] = bkginputdata.hit_time[rand];
        bkgresultlist->hit_xpos[bkgresultlist->numhits] = bkginputdata.hit_xpos[rand];
        bkgresultlist->hit_ypos[bkgresultlist->numhits] = bkginputdata.hit_ypos[rand];
        bkgresultlist->numhits++;
        rand++;
        if(bkgresultlist->numhits == arrsize) {
          arrsize += 50;
          bkgresultlist->hit_energy = (double*) realloc(bkgresultlist->hit_energy, arrsize * sizeof(double));
          bkgresultlist->hit_time = (double*) realloc(bkgresultlist->hit_time, arrsize * sizeof(double));
          bkgresultlist->hit_xpos = (double*) realloc(bkgresultlist->hit_xpos, arrsize * sizeof(double));
          bkgresultlist->hit_ypos = (double*) realloc(bkgresultlist->hit_ypos, arrsize * sizeof(double));
        }
      } while((rand < bkginputdata.numrows) && (bkginputdata.hit_time[rand - 1] == bkginputdata.hit_time[rand]));
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
