/***********************************************************************
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

   Copyright 2014:  INOUTUTILS has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      INOUTUTILS
*
*  File:       inoututils.h
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef INOUTUTILS_
#define INOUTUTILS_

// Utils module

	#include "genutils.h"
	
	using namespace std;

        /*#ifdef __cplusplus
        extern "C"
        #endif*/
        struct IOData
	{
		fitsfile *inObject;
		char *nameTable;
		char *nameCol;
		char *unit;
		//MC char *type;
		int type;
		int iniCol;
		int endCol;
		long iniRow;
		long endRow;
	};
	
	/*typedef struct IOData
	{
		fitsfile *inObject;
		char *nameTable;
		char *nameCol;
		char *unit;
		//MC char *type;
		int type;
		int iniCol;
		int endCol;
		long iniRow;
		long endRow;
                #ifdef __cplusplus
                    IOData();
                    IOData(const IOData& other);
                    IOData& operator=(const IOData& other);
                    ~IOData();
                #endif
	} IOdata;*/
        
	// Structure to define input parameters
	typedef struct
	{
		string name;
		string description;
		string type;
		string ValStr;
		string defValStr;
		int ValInt;
		int minValInt;
		int maxValInt;
		int defValInt;
		double ValReal;
		double minValReal;
		double maxValReal;
		double defValReal;
	} inparam;

        #ifdef __cplusplus
        extern "C"
        #endif
	int readFitsSimple(IOData obj, gsl_vector **result);
	int readFitsComplex(IOData obj, gsl_matrix **result);

	int writeFitsSimple (IOData obj,gsl_vector *vector);
	int writeFitsComplex(IOData obj, gsl_matrix *matrix);

	int toGslMatrix(void **buffer, gsl_matrix **matrix, long numCol,int numRow,int type, int eventini);
	int toGslVector(void **buffer, gsl_vector **array, long nevent, int eventini, int type);

	int fromGslVector(void **buffer, gsl_vector **array, int type);
	int fromGslMatrix(void **buffer, gsl_matrix **matrix, int type);

	int interactivePars(inparam *taskPars, int np, string task);

#endif /*INOUTUTILS_*/
