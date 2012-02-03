#ifndef SIMPUTSRC_H
#define SIMPUTSRC_H 1

#include "sixt.h"
#include "simput.h"

#define TOOLSUB simputsrc_main
#include "headas_main.c"


struct Parameters {

  /** Power law. */
  float plPhoIndex;
  float plFlux;

  /** Black body temperature [keV]. */
  float bbkT;
  float bbFlux;

  /** Line dispersion [keV]. */
  float flSigma;
  float flFlux;

  float rflSpin;
  float rflFlux;

  /** Absorption column [10^22 atoms/cm^2] */
  float NH;

  /** Reference energy band [keV]. */
  float Emin;
  float Emax;

  /** PSD general parameters */
  long PSDnpt;
  float PSDfmin;
  float PSDfmax;

  /** PSD: Zero-frequency Lorentzian parameters */
  float LFQ;
  float LFrms;

  /** PSD: Horizontal branch Lorentzian parameters */
  float HBOf;
  float HBOQ;
  float HBOrms;

  /** PSD: Quasi-periodic Lorentzian parameters (1-3) */
  float Q1f;
  float Q1Q;
  float Q1rms;

  float Q2f;
  float Q2Q;
  float Q2rms;

  float Q3f;
  float Q3Q;
  float Q3rms;

  /** File name of the input ISIS parameter file containing a spectral
      model. */
  char ISISFile[MAXFILENAME];

  /** File name of the input ASCII spectrum. */
  char XSPECFile[MAXFILENAME];

  /** Source position [deg]. */
  float RA;
  float Dec;

  /** Name of the X-ray source. */
  char Src_Name[MAXMSG];

  /** File name of the output SIMPUT file. */
  char Simput[MAXFILENAME];
  
  char clobber;
};


int simputsrc_getpar(struct Parameters* const par);


#endif /* SIMPUTSRC_H */

