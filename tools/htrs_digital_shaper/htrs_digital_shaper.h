#ifndef HTRS_DIGITAL_SHAPER_H
#define HTRS_DIGITAL_SHAPER_H 1

#include "sixt.h"
#include "htrseventfile.h"
#include "htrsevent.h"


#define TOOLSUB htrs_digital_shaper_main
#include "headas_main.c"


//  Layout of the digital solution:
//
//    -----          --------------          -----          ---------
//   | SDD |  --->  | Preamplifier |  --->  | ADC |  --->  | Filters |
//    -----          --------------          -----          ---------


//////////////////////////////////////////////////////////////////


/** Number of pixels / channels. */
#define NCHANNELS 31

/** Number of time bins in the history of each ADC. */
#define NTIMEBINS 1000


//////////////////////////////////////////////////////////////////


/** Input voltage range of the ADC [V]. */
const double adc_input_range = 2.;

/** SDD output [V/keV]. */
const double sdd_output = 1.2e-3;

/** Number of ADC channels. */
const long adc_N_channels = 16384; // 2^14;

/** ADC reset threshold. If the ADC channel is greater or equal this
    threshold, a reset will be performed. */
const long adc_reset_threshold = 16384 - 500; // 2^14 - 500

/** Duration of a reset of the SDD [s]. */
const double reset_duration = 500.e-9;


//////////////////////////////////////////////////////////////////


struct Parameters{
  /** Filename of the input HTRS event file containing all events. */
  char input_eventlist_filename[MAXMSG];
  /** Filename of the output HTRS event file containing only events
      that are properly measured by the digital shaper. */
  char output_eventlist_filename[MAXMSG];
  /** Filename of the template for a new HTRS event file. */
  char eventlist_template[MAXMSG];

  /** Frequency of the ADC [Hz]. */
  double frequency;
  /** Gain of the preamplifier. */
  double preamp_gain;
  /** Number of sampling points. */
  int nsamplings;
};


//////////////////////////////////////////////////////////////////


int htrs_digital_shaper_getpar(struct Parameters*);


#endif /* HTRS_DIGITAL_SHAPER_H */
