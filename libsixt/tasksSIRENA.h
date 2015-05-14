/************************************************************************************************

   This file is part of SIXTE/SIRENA software.

   SIXTE/SIRENA is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE/SIRENA is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

   Copyright 2014:  This file has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

************************************************************************************************/

#ifndef TASKSSIRENA_H
#define TASKSSIRENA_H 1

#include "integraSIRENA.h"
#include "pulseprocess.h"

void runDetect(TesRecord* record, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord);

int handleLibraryDetect(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, int lastRecord, fitsfile **inLibObject, const char * create);
int createLibraryDetect(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject, const char * create);
int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, const char * create, fitsfile **dtcObject);
int filderLibrary(ReconstructInitSIRENA** reconstruct_init, double samprate);
int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double);
int procRecord(ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, PulsesCollection *foundPulses);
int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl);
int writePulses(ReconstructInitSIRENA** reconstruct_init, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, fitsfile *dtcObject);
int calculateTemplate (ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight);
int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhisto, gsl_vector **yhisto);
int align(double samprate, gsl_vector **vector1, gsl_vector ** vector2);
int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m);
int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m);
int writeLibrary(ReconstructInitSIRENA* reconstruct_init, double estenergy, gsl_vector **pulsetemplate, bool appenToLibrary, fitsfile **inLibObject);

void runFilter(TesRecord* record, int nRecord, int lastRecord, ReconstructInitSIRENA** reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter);

int find_energy(double energyKeyword, gsl_vector *energygsl, long *rowFound);
int calculus_optimalFilter(int TorF, int intermediate, int mode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, gsl_vector *csdgsl, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex, double *normalizationFactor);
int interpolatePOS (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out, long *numzerosstart, long *numzerosend);
int find_matchedfilter(int runF0orB0val, double ph, gsl_vector *energiesvalues, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound);
int calculateUCEnergy (gsl_vector *vector, gsl_vector *filter, gsl_vector_complex *filterFFT, int domain, int mode, double nrmfctr, double samprate, double *calculatedEnergy);
int writeFilter(ReconstructInitSIRENA *reconstruct_init, double normalizationFactor, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject, const char *create);
int writeUCEnergy(ReconstructInitSIRENA **reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, int pulse_index, double uncE, fitsfile **dtcObject, const char *create);
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double normalizationFactor, double uncE, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject, const char * create);

void runEnergy(ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord);

int loadUCEnergies(ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, long *nz, gsl_vector **zi, double *E0z);
int calculus_bc (int calibLQ, long nx, gsl_vector *xi, double E0x, long ny, gsl_vector *yj, double E0y, double *b_cF, double *c_cF);
int calculateEnergy(ReconstructInitSIRENA *reconstruct_init, PulsesCollection **pulses);
int writeEnergy(ReconstructInitSIRENA **reconstruct_init,PulsesCollection *pulsesInRecord, const char *create);



#endif /* TASKSSIRENA_H */




