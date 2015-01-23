/************************************************************************************************
   This file is part of SIXTE/SIRENA software.

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

   Copyright 2014:  Trigger has been developed by the INSTITUTO DE FISICA DE
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/*************************************************************************************************
*                 ENERGY
*
*  File:      energy.h
*  Developer: Beatriz Cobo Martí­n
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
*  19/10/06    First version
*  12/12/06	   Adding holzgauss module to energyresol module
*  29/01/06	   Adding new function cutFunction
*  24/04/07    Include MAX_CHANNEL_COUNT: a value for the maximum number of data channels in an input file
*  12/03/08	   Holzgauss reads from ers fits file 2.0.0
*  26/06/08	   Adding input/output parameters "ql" and "nBins".
* 			   New documentation in the source code.
* 			   Adding output keyword: The result of the fits "Scale Factor, Offset of Energy and fwhm"
* 			   Adding output keyword: "ENERGYRESOL"
* 			   Adding column in the hzg FITS file "Histo" and "HistoFit"
* 			   Change fuction fitGauss by getParams
* 			   Adding fuction fit_poly
*  10/09/08	   Adding input parametrer "nameLog"
* 			   Adding input parametrer "verbosity" and processing it using RIL libraries.
* 			   Deleting in/output keyword "channelcount"
*  29/09/08	   Included "CREATE" and "PROCESS" keyword in the output FITS file.
*			   Included new input parameter: deltaBase
*  14/01/09	   Included new input keywords
*  02/02/09	   Included new input keywords to MONOCHROMATIC and PSEUDO modes
* 			   Included new input keywords "ec". Energy Calibration.
* 			   Included non valid parameter
*  20/02/09	   Included, modificated and removed input keywords
*  24/02/09	   Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
*  12/03/09	   Changed documentation
*              Used the keywords libraries
*  28/04/09    The input parameter fitabstol is used
*              Added Emin and Emax and change deltaBase to deltaKbeta
*  03/06/09    Added Tmin and Tmax (input parameters)
* 			   Added analyzeDrift, fitPseudo, fcnChisq, gaussShort and trapezeRule functions
*  23/06/09    Added lorenzian, humlicek_v12, voigt, holzer
*              Deleted non-used functions
*  30/11/09    Adding "processin" and "history" variables
*  25/02/11	   Now "proccesin" variable has more memory space
*  11/03/11    "valuegsl2" added
*  25/03/11    Deleted "ql" variable
*  29/03/11    Updated .h
*  05/04/11    New variable "noPulsesTrg"
*  09/11/11    _psh input FITS file has a new column QUALITY in EUR-ENR => It is read => New variables 'qualitygsl' and 'qualitygsl2'
*  14/11/11    New input parameters: 'b_cF' and 'c_cF' (calibration factors)
*  22/07/14	   fcnChisq changed (.../F_e instead .../sigmasq=.../yhistogsl)
*  24/07/14	   Histogram without null bars used in the calculus BUT histogram with null bars in the output FITS file:
*  			   xhistoInitialgsl, yhistoInitialgsl, AeiInitialgsl
*  31/07/14    Code modified in order to work with non monochromatic input FITS files
*  	 Sep/14    Analytical calculus of b and c and sigma (no fitting)
*    Dec/14    DAL->CFITSIO migration + error routines update
*  12/01/15    Deleted 'processinA' and 'processinB'
*              Deleted 'annalsinA' and 'annalsinB'
*  22/01/15    Renaming of input parameter HorMorL --> gradesForER
******************************************************************************************/

#ifndef ENERGY_H_
#define ENERGY_H_

// Utils module

	#include <inoututils.h>

	int status = EPOK;

// CFITSIO helpers

	int  colnum=0, felem=0, keyvalint=0;
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];

// Constants

	const double PI = 4.0 * atan(1.0);

// INPUT/OUTPUT FILES

	fitsfile *psgObjectA=NULL; 		// Object which contains information of the input FITS file (inFileA)
	char psgNameA[255];				// Name of the input FITS file (inFileA)

	fitsfile *psgObjectB=NULL; 		// Object which contains information of the input FITS file (inFileB)
	char psgNameB[255];				// Name of the input FITS file (inFileB)

	FILE *fileRef;					// Pointer for file which contains errors and warnings

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Provisional => To be deleted in future
	FILE * temporalFile;
	char temporalFileName[255];
	char val[256];
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	char *unit=NULL, *comment=NULL;

// INPUT KEYWORDS

	//PSGRADE extension (inFileA)
	long eventcntA;
	long eventcntA_OK;
	long eventcntA1_OK;
	long eventcntA2_OK;
	long eventcntA3_OK;
	double energyA;

	//PSGRADE extension (inFileB)
	long eventcntB;
	long eventcntB_OK;
	long eventcntB1_OK;
	long eventcntB2_OK;
	long eventcntB3_OK;
	double energyB;

//INPUT VECTORS

	gsl_vector *tstartgslA;
	gsl_vector *pseudoenergygslA;
	gsl_vector *pseudoenergygslA_mod;
	gsl_vector *qualitygslA;
	gsl_vector *gradegslA;

	gsl_vector *tstartgslB;
	gsl_vector *pseudoenergygslB;
	gsl_vector *pseudoenergygslB_mod;
	gsl_vector *qualitygslB;
	gsl_vector *gradegslB;

//INPUT PARAMETERS

	char nameLog[255];		// Output log file name
	int verbosity;			// Verbosity level of the output log file
	char clobberStr[4];		// Clobber=yes then overwritte output files	
	int clobber=0;
	// Calibration factor to convert pseudoenergies into energies
	double b_cF, c_cF;				// E = b_cFÂ·e + c_cF^2	(e->pseudoenergy, E->energy)

	int optmode;
	int HorMorL;
	char gradesForER[4];		// Event grades to be used for Energy resolution calculation

//AUXILIARY VARIABLES

	// To process each row of pulses in inDataIterator
	int drift = 0;

	// To avoid the deprecate conversions
	char *straux = new char[255];

	// Parameters used to write output FITS files
	IOData obj;

// OUTPUT VECTORS

	gsl_vector *energygslA;
	gsl_vector *energygslA1;
	gsl_vector *energygslA2;
	gsl_vector *energygslA3;
	gsl_vector *energygslB;
	gsl_vector *energygslB1;
	gsl_vector *energygslB2;
	gsl_vector *energygslB3;

	double sigmaA1;
	double sigmaA2;
	double sigmaA3;
	double sigmaB1;
	double sigmaB2;
	double sigmaB3;

// OUTPUT KEYWORDS

	const char *create;					// Name and version of the module: name-0.0.0

//FUNCTIONS

	int initModule(int argc, char **argv);

	int inDataIteratorA(long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct);
	int inDataIteratorB(long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct);

	int calculus_bc (long nx, long nx_OK, gsl_vector *xi, double E0x, long ny, long ny_OK, gsl_vector *yj, double E0y, double *b_cF, double *c_cF);
	int calculus_sigma (long n, long n_OK, gsl_vector *energygsl, double E0, double b, double c, double *sigma);

	using namespace std;

#endif /*ENERGY_H_*/
