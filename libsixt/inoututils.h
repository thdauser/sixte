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

   Copyright 2014:  Trigger has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/***********************************************************************
*                      INOUTUTILS    			                      
*                                                                     
*  File:      inoututils.h
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*                                                                    
*  Revision History:                                                  
*                                                                       
*  version 1.0.0: 21/09/06     	First version   
*  version 1.0.1: 03/04/08		Changing function toGslVector: 
* 								The input parameter "type "had been include. 
* 								Type includes the kind of datas of the buffer.  
*  version 1.4.0: 16/09/08		Modify the function "writeFitsSimple2" and "int writeFitsComplex2" 
* 								Adding funcion called "writeLog" to processe of the each level message of log file. 
*  version 	1.6.0	06/09/08	Adding new library "stdlib"
*  29/03/11 Updated .h
*  08/07/14	Free of PIL and RIL
*           New EPOK
*           'writeLog' modified
***********************************************************************/

#ifndef INOUTUTILS_
#define INOUTUTILS_

// Utils module

	#include "genutils.h"

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

	int readFitsSimple(IOData obj, gsl_vector **result);
	int readFitsComplex(IOData obj, gsl_matrix **result);
	int readFitsImage(IOData obj, gsl_matrix **result);

	int writeFitsSimple (IOData obj,gsl_vector *vector);
	int writeFitsComplex(IOData obj, gsl_matrix *matrix);

	int toGslMatrix(void **buffer, gsl_matrix **matrix, long numCol,int numRow,int type, int eventini);
	int toGslVector(void **buffer, gsl_vector **array, long nevent, int eventini, int type);

	int fromGslVector(void **buffer, gsl_vector **array, int type);
	int fromGslMatrix(void **buffer, gsl_matrix **matrix, int type);

	int interactivePars(inparam *taskPars, int np, string task);

using namespace std;

#endif /*INOUTUTILS_*/
