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


   Copyright 2016 Philippe Peille, CNES
   Copyright 2017-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "tessim.h"

#include <math.h>
#include <complex.h>

// Constructor of a timedomain_bbfb_info structure
timedomain_bbfb_info* init_timedomain_bbfb(double tclock,int delay,double carrier_frequency,
    double gbw,double phase,double bias_leakage,double bias_leakage_phase,
    double fb_leakage,int fb_leakage_lag,double fb_start,int* status){
  timedomain_bbfb_info *bbfb=(timedomain_bbfb_info*)malloc(sizeof(*bbfb));
  CHECK_NULL_RET(bbfb,*status,"Memory allocation failed in init_timedomain_bbfb",NULL);

  bbfb->delay = delay;
  bbfb->tclock = tclock;
  bbfb->time_delay = delay*tclock;
  bbfb->carrier_frequency = carrier_frequency;
  bbfb->gbw = gbw;
  bbfb->phase = phase; // phase of the carrier. Is useless until several carriers are in the same BBFB loop
  bbfb->cphase = 2*M_PI*bbfb->time_delay*carrier_frequency;
  bbfb->integral=fb_start*sqrt(2)/(bbfb->tclock*bbfb->gbw*2*M_PI*2); // necessary to compensate for artificial perfect lock at simulation start
  bbfb->bias_leakage=bias_leakage;
  bbfb->bias_leakage_phase=bias_leakage_phase;
  bbfb->fb_leakage=fb_leakage;
  bbfb->fb_leakage_lag=fb_leakage_lag;

  bbfb->fb_values = (double*)malloc(delay*sizeof(double));
  CHECK_NULL_RET(bbfb->fb_values,*status,"Memory allocation failed in init_timedomain_bbfb",NULL);

  // Initialize FB values with perfect carrier nulling. Necessary to lock FB loop at simulation start
  for (int ii=0;ii<delay;ii++){
    bbfb->fb_values[ii]=fb_start*sqrt(2)*cos(2*M_PI*bbfb->carrier_frequency*(ii+1)*tclock+bbfb->phase);
  }
  bbfb->fb_index=0;

  return (bbfb);
}

// Destructor of timedomain_bbfb_info structure
void free_timedomain_bbfb(timedomain_bbfb_info **bbfb,int *status) {
  CHECK_STATUS_VOID(*status);
  if (*bbfb!=NULL){
    free((*bbfb)->fb_values);
  }
  *bbfb=NULL;
}

// Run BBFB loop clock
// TODO: this function is effectively valid only for one pixel. It will need to be adapted for a multiple pixels loop
double run_timedomain_bbfb_loop(tesparams *tes, double time, double squid_input,double squid_noise_value,gsl_rng *rng){
  timedomain_bbfb_info* bbfb = (timedomain_bbfb_info*) tes->bbfb_info;

  // TODO: precompute lod and lor ? -> we only have ten numbers per period

  // Modulation and demodulation numbers
  double _Complex lod = cexp(I*(2*M_PI*bbfb->carrier_frequency*time+bbfb->phase));
  double _Complex lor = cexp(-I*(2*M_PI*bbfb->carrier_frequency*time+bbfb->phase+bbfb->cphase));

  // Compute random leakage phases if needed (note that this is not physical but random phase can be used as a worst case: leakage is treated as noise)
  // Signal leakage from FB line to SQUID output
  double fb_leakage_value;
  if (bbfb->fb_leakage_lag==-999){
    gsl_ran_choose(rng, &fb_leakage_value, 1, bbfb->fb_values, bbfb->delay, sizeof (double));
    fb_leakage_value*=bbfb->fb_leakage;
  } else {
    fb_leakage_value = bbfb->fb_leakage*bbfb->fb_values[(bbfb->fb_index-bbfb->fb_leakage_lag+bbfb->delay) % bbfb->delay];
  }
  // Signal leakage from bias line to SQUID output
  double bias_leakage_value;
  if (bbfb->bias_leakage_phase==-999){
    bias_leakage_value = bbfb->bias_leakage*cos(2*M_PI*bbfb->carrier_frequency*time+bbfb->phase+gsl_rng_uniform(rng)*2*M_PI)*sqrt(2);
  } else {
    bias_leakage_value = bbfb->bias_leakage*cos(2*M_PI*bbfb->carrier_frequency*time+bbfb->phase+bbfb->bias_leakage_phase)*sqrt(2);
  }

  // Modulate input
  squid_input*=cos(2*M_PI*bbfb->carrier_frequency*time+bbfb->phase)*sqrt(2); // change from rms to amplitude

  // Compute error signal at squid input and add noise
  double current_error = squid_input-bbfb->fb_values[bbfb->fb_index]+bias_leakage_value+fb_leakage_value+squid_noise_value;

  // Error signal in phi_0 units
  tes->squid_error=tes->M_in*current_error/tes->TTR;

  // SQUID transfer function simulated as a sine wave + gain correction to keep signal in units of TES current
  double error=sin(2*M_PI*tes->squid_error)/(2*M_PI)/tes->M_in*tes->TTR;

  // Integrate signal
  bbfb->integral+= error*lod; // no need to add - delay as error_index is modulo delay

  // Compute feedback and update values for next loop
  double _Complex base = bbfb->integral*bbfb->tclock*bbfb->gbw*2*M_PI*2; // x2 is to compensate modulation loss and have the correct GBW product
  //printf("t: %g e: %g phi: %g fb: %g int: %g base: %g\n",time,current_error,tes->squid_error,bbfb->fb_values[bbfb->fb_index],creal(bbfb->integral),cabs(base));
  bbfb->fb_values[bbfb->fb_index] = creal(base*lor);

  // Go to next index for feedback values
  bbfb->fb_index = (bbfb->fb_index+1) % bbfb->delay;

  // TODO: This corresponds to IQ demodulation, creal would correspond to I demodulation (might want to have this as an option using readoutMode)
  return(cabs(base));
}
