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

void runDetect(TesRecord* record, int lastRecord_runDetect, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord);

int handleLibraryDetect(ReconstructInitSIRENA* reconstruct_init_handleLibraryDetect, bool *appendToLibrary_handleLibraryDetect, int lastRecord_handleLibrary, fitsfile **inLibObject_handleLibraryDetect, const char * create_handleLibraryDetect, gsl_matrix **library_handleLibraryDetect, gsl_matrix **models_handleLibraryDetect);
int initLibraryDetect(ReconstructInitSIRENA* reconstruct_init_initLibraryDetect, gsl_matrix **library_initLibraryDetect, gsl_matrix **models_initLibraryDetect);
int loadLibraryDetect(ReconstructInitSIRENA* reconstruct_init_loadLibraryDetect, gsl_matrix **library_loadLibraryDetect, gsl_matrix **models_loadLibraryDetect);
int createLibraryDetect(ReconstructInitSIRENA* reconstruct_init_createLibraryDetect, bool *appendToLibrary_createLibraryDetect, fitsfile **inLibObject_createLibraryDetect, const char * create_createLibraryDetect, gsl_matrix **library_createLibraryDetect, gsl_matrix **models_createLibraryDetect);
int createDetectFile(ReconstructInitSIRENA* reconstruct_init_createDetectFile, double samprate, const char * create_createDetectFile, fitsfile **dtcObject_createDetectFile);
int filderLibrary(ReconstructInitSIRENA** reconstruct_init_filderLibrary, double samprate, gsl_matrix **models_filderLibrary);
int loadRecord(TesRecord* record_loadRecord, double *time_record, gsl_vector **adc_double);
int procRecord(ReconstructInitSIRENA** reconstruct_init_procRecord, double tstartRecord_procRecord, double samprate, fitsfile *dtcObject_procRecord, gsl_vector *record, gsl_matrix *library_procRecord, gsl_matrix *models_procRecord, PulsesCollection *foundPulses);
int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl);
int writePulses(ReconstructInitSIRENA** reconstruct_init_writePulses, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *pulseheights, fitsfile *dtcObject_writePulses);
int calculateTemplate (ReconstructInitSIRENA *reconstruct_init_calculateTemplate, PulsesCollection *pulsesAll_calculateTemplate, PulsesCollection *pulsesInRecord_calculateTemplate, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight);
int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhisto, gsl_vector **yhisto);
int align(double samprate, gsl_vector **vector1, gsl_vector ** vector2);
int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m);
int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m);
int writeLibrary(ReconstructInitSIRENA* reconstruct_init_writeLibrary, double estenergy, gsl_vector **pulsetemplate, bool appenToLibrary_writeLibrary, fitsfile **inLibObject_writeLibrary);

void runFilter(TesRecord* record, int nRecord_runFilter, int lastRecord_runFilter, ReconstructInitSIRENA** reconstruct_init, PulsesCollection *pulsesAll_runFilter, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter_runFilter);

int initLibraryFilter(ReconstructInitSIRENA* reconstruct_init_initLibraryFilter, gsl_vector **energylibrary_initLibraryFilter, gsl_vector **estenergylibrary_initLibraryFilter,
		gsl_matrix **templateslibrary_initLibraryFilter, gsl_matrix **templateslibraryb0_initLibraryFilter,
		gsl_matrix **matchedfilters_initLibraryFilter, gsl_matrix **matchedfiltersb0_initLibraryFilter);
int loadLibraryFilter(ReconstructInitSIRENA* reconstruct_init_loadLibraryFilter, gsl_vector **energylibrary_initLibraryFilter, gsl_vector **estenergylibrary_initLibraryFilter,
		gsl_matrix **templateslibrary_initLibraryFilter, gsl_matrix **templateslibraryb0_initLibraryFilter,
		gsl_matrix **matchedfilters_initLibraryFilter, gsl_matrix **matchedfiltersb0_initLibraryFilter);
int initNoisespectrum(ReconstructInitSIRENA* reconstruct_init_initNoisespectrum, gsl_vector **freqgsl_initNoisespectrum, gsl_vector **csdgsl_initNoisespectrum);
int loadNoisespectrum(ReconstructInitSIRENA* reconstruct_init_loadNoisespectrum, gsl_vector **freqgsl_loadNoisespectrum, gsl_vector **csdgsl_loadNoisespectrum);
int find_energy(double energyKeyword, gsl_vector *energygsl, long *rowFound);
int calculus_optimalFilter(gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val_calculusOptimalFilter, gsl_vector *freqgsl_calculusOptimalFilter, gsl_vector *csdgsl_calculusOptimalFilter, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, double *normalizationFactor_calculusOptimalFilter);
int interpolatePOS (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out, long *numzerosstart, long *numzerosend);
int getMatchedFilter(int runF0orB0val_getMatchedFilter, double pheight, gsl_vector *estenergylibrary_getMatchedFilter, gsl_matrix *matchedfilters_getMatchedFilter, gsl_matrix *matchedfiltersb0_getMatchedFilter, gsl_vector **matchedfilter_getMatchedFilter);
int find_matchedfilter(double ph, gsl_vector *energiesvalues, gsl_matrix *matchedfilters, gsl_vector **matchedfilterFound, FILE * temporalFile);
int calculateUCEnergy (gsl_vector *vector, gsl_vector *filter, int domain, double nrmfctr, double samprate, double *calculatedEnergy);
int writeFilter(ReconstructInitSIRENA *reconstruct_init_writeFilter, double normalizationFactor_writeFilter, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject_writeFilter, const char *create_writeFilter);
int writeUCEnergy(ReconstructInitSIRENA **reconstruct_init_writeUCEnergy, PulsesCollection *pulsesAll_writeUCEnergy, PulsesCollection *pulsesInRecord_writeUCEnergy, int pulse_index, double uncE_writeUCEnergy, fitsfile **dtcObject_writeUCEnergy, const char *create_writeFilterHDU);
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init_writeFilterHDU, int pulse_index, double normalizationFactor_writeFilterHDU, double uncE_writeFilterHDU, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject_writeFilterHDU, const char * create_writeFilterHDU);

void runEnergy(ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord);

int loadUCEnergies(ReconstructInitSIRENA *reconstruct_init_loadUCEnergies, PulsesCollection *pulsesAll_loadUCEnergies, long *nz_loadUCEnergies, gsl_vector **zi_loadUCEnergies, double *E0z_loadUCEnergies);
int calculus_bc (int calibLQ, long nx, gsl_vector *xi, double E0x, long ny, gsl_vector *yj, double E0y, double *b_cF, double *c_cF);
int calculateEnergy(ReconstructInitSIRENA *reconstruct_init_calculateEnergy, PulsesCollection **pulses);
int writeEnergy(ReconstructInitSIRENA **reconstruct_init_writeEnergy,PulsesCollection *pulsesInRecord_writeEnergy, const char *create_writeEnergy);



#endif /* TASKSSIRENA_H */




