#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fitsio.h"
#include "hdcal.h"
#include "HDgtcalf_internal.h"

#ifdef HDGTCALF_STANDALONE
#include "HDgtcalf_standalone.h"
#else
#include "headas_utils.h"
#include "headas_error.h"
#endif

/*
HISTORY
-------
  Version 1.0 written by Ziqin Pan, NASA/GSFC, MAY 2003
*/

/* Function Prototypes  */

static int gtcalf_work(char* cif, const char* tele, const char* instr,
      const char * detnam, const char* filt, const char* codenam,
      double startreftime, double stopreftime, const char *expr,
      int maxret,int fnamesize, char** file, long* extno,
      char** online, int* nret, int* nfound, int* status);
static int gtcalidx (char* mode, char* missn, char* inst, char* ciffil,
      char* instdir, int* status);
static int rdcnfgl (FILE * fptr, char * missn, char * inst, char* cifdev,
      char* cifdir, char* cif, char* datadev, char* datadir, int * status);


/*---------------------------------------------------------------------*/
int HDgtcalf(const char* tele, const char* instr,
           const char* detnam, const char* filt,
           const char* codenam, const char* strtdate,
           const char* strtime, const char* stpdate, 
           const char* stptime, const char* expr,
           int maxret, int fnamesize, char** filenam,
	   long* extno, char** online, int* nret, 
           int* nfound, int* status)
/* Get CALDB file with given criteria */
{

     char cif[FILENAME_MAX+1],instdir[FILENAME_MAX+1];
     double djm, startreftime, stopreftime, dayfrac;
     char msg[256];

     if(NULL == status)
         return HD_ERROR_THROW("Input status variable is NULL", HD_ERR_NULL_POINTER);

     if(*status) return *status;

     if( strtdate == NULL ||
         strtime == NULL ||
         stpdate == NULL ||
         stptime == NULL || 
         tele == NULL ||
         instr == NULL ||
         detnam == NULL ||
         codenam == NULL ||
         filt == NULL || 
         expr == NULL)  
     {
         *status =HD_ERR_NULL_POINTER;
         sprintf(msg,"Input parameter string is NULL");
         HD_ERROR_THROW(msg,*status); 
	 goto cleanup;
     }


     /* Get CALDB index file */

     if(gtcalidx("INST",(char*) tele,(char*) instr,cif,instdir,status)) {
        sprintf(msg,"Problem getting caldb index file");
        HD_ERROR_THROW(msg,*status); 
        goto cleanup;
     }


     /* Parse date and time string to year, month, day, hour, minute and
       second and convert it to mjd */

     if (strcmp(strtdate,"-") == 0) {
        djm =-99;
     } 
     else {
        if(dt2mjd(strtdate,&djm,status)) {
            sprintf(msg,"Unable to determine MJD from startdate string");
            HD_ERROR_THROW(msg,*status); 
            goto cleanup;
        }
     }


     if(strcmp(strtime,"-") == 0) {
        dayfrac =0;
     } 
     else {
        if(tim2df(strtime,&dayfrac,status)) {
            sprintf(msg,"Unable to determine fraction of day from startime string");
            HD_ERROR_THROW(msg,*status); 
            goto cleanup;
        }
     }

     startreftime = djm+dayfrac;

     if(strcmp(stpdate,"-") == 0) {
        djm =-99;
     } 
     else {
        if(dt2mjd(stpdate,&djm,status)) {
            sprintf(msg,"Unable to determine MJD from stopdate string");
            HD_ERROR_THROW(msg,*status); 
            goto cleanup;
        }
     }


     if(strcmp(stptime,"-") == 0) {
        dayfrac =0;
     } 
     else {
        if(tim2df(stptime,&dayfrac,status)) {
            sprintf(msg,"Unable to determine fraction of day from stoptime string");
            HD_ERROR_THROW(msg,*status); 
            goto cleanup;
        }
     }

     stopreftime = djm+dayfrac;


     /* call work function to select CALDB file with the input criteria */

     if(gtcalf_work(cif,tele,instr,detnam,filt,codenam,
         startreftime,stopreftime,expr,maxret, fnamesize,
         filenam,extno,online,nret,nfound,status)) {
         sprintf(msg,"Fail to select CALDB file with the given criteria");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
     } 

cleanup:
     return *status;
}
/*------------------------------------------------------------------------*/

static int gtcalf_work(char* cif, const char* tele, const char* instr, 
         const char * detnam, const char* filt, const char* codenam,
	 double startreftime, double stopreftime,
         const char *expr, int maxret, int fnamesize,
         char** file, long* extno,char** online,
	 int* nret, int* nfound, int* status)
{

    fitsfile* fptr=NULL;
    long nax2val=0;
    int telcol=0, inscol=0, detcol=0,codcol=0;
    int filtcol=0, cbdcol=0, refcol=0, qulcol=0;
    int devcol=0, dircol=0, filcol=0, extcol=0;

    int dw=0;
    int cbdnelem = 9;
    char *telval=NULL, *insval=NULL, *filtval=NULL, *cbdval=NULL;
    char *codval=NULL, *detval=NULL;
    char *dirval=NULL;
    double *refval=NULL,qulval=0,maxrefval=0;
    int anynull =1;
    char *tmpfile=NULL;

    CBDLIST *list1 =NULL;
    CBDLIST *list2 =NULL;

    char * caldbvar ="CALDB"; 

    int *index=NULL;
    char msg[256];

    int i;
    int k;
    int typecode;
    long repeat=1, width=1;


    if(*status) return *status;

    tmpfile=(char*) calloc(FILENAME_MAX+1,sizeof(char));
    
    
    if (parseCBD2(expr,&list1,status) ) {
        sprintf(msg,"Parse boolean expression error");
        HD_ERROR_THROW(msg,*status); 
        goto cleanup;
    }


    if (fits_open_table(&fptr, cif, READONLY, status)) {
        sprintf(msg,"Unable to open CALDB index file");
        HD_ERROR_THROW(msg,*status); 
        goto cleanup;
     }

    if (ffgkyj(fptr,"NAXIS2",&nax2val,NULL,status)) {
        sprintf(msg,"Problem reading NAXIS2 of CIF");
        HD_ERROR_THROW(msg,*status); 
        goto cleanup;
    }

    if( nax2val == 0 ) {
         *status = CIF_ERR_EMPTY_TABLE;
         sprintf(msg,"CIF table is empty");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }

    index = (int*) malloc(nax2val*sizeof(int)); 
    refval = (double*) calloc(nax2val,sizeof(double));

    for( i=0; i<nax2val; i++) {
	index[i]=0;
    }

    
    if( ffgcno(fptr,1,"TELESCOP",&telcol,status)) {
         sprintf(msg,"Unable to find TELESCOP column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if(ffgcdw(fptr,telcol,&dw,status)) {
         sprintf(msg,"Unable to find TELESCOP column display width");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }
    else {
         telval=(char*) calloc(dw+1,sizeof(char));
    }

   
    if( ffgcno(fptr,1,"INSTRUME",&inscol,status)) {
         sprintf(msg,"Unable to find INSTRUME column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if(ffgcdw(fptr,inscol,&dw,status)) {
         sprintf(msg,"Unable to find INSTRUME column display width");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }
    else {
         insval=(char*) calloc(dw+1,sizeof(char));
    }

    if( ffgcno(fptr,1,"DETNAM",&detcol,status)) {
         sprintf(msg,"Unable to find DETNAM column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if(ffgcdw(fptr,detcol,&dw,status)) {
         sprintf(msg,"Unable to find DETNAM column display width");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }
    else {
         detval=(char*) calloc(dw+1,sizeof(char));
    }

    if( ffgcno(fptr,1,"FILTER",&filtcol,status)) {
         sprintf(msg,"Unable to find FILTER column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if(ffgcdw(fptr,filtcol,&dw,status)) {
         sprintf(msg,"Unable to find FILTER column display width");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }
    else {
         filtval=(char*) calloc(dw+1,sizeof(char));
    }

    if( ffgcno(fptr,1,"CAL_CNAM",&codcol,status)) {
         sprintf(msg,"Unable to find CAL_CNAM column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }
    if(ffgcdw(fptr,codcol,&dw,status)) {
         sprintf(msg,"Unable to find CAL_CNAM column display width");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }
    else {
        codval=(char*) calloc(dw+1,sizeof(char));
    }

    if( ffgcno(fptr,1,"REF_TIME",&refcol,status)) {
         sprintf(msg,"Unable to find REF_TIME column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if( ffgcno(fptr,1,"CAL_QUAL",&qulcol,status)) {
         sprintf(msg,"Unable to find CAL_QUAL column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if( ffgcno(fptr,1,"CAL_DEV",&devcol,status)) {
         sprintf(msg,"Unable to find CAL_DEV column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }


    if( ffgcno(fptr,1,"CAL_DIR",&dircol,status)) {
         sprintf(msg,"Unable to find CAL_DIR column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if(ffgcdw(fptr,dircol,&dw,status)) {
         sprintf(msg,"Unable to find CAL_DIR column display width");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }
    else {
         dirval=(char*) calloc(dw+1,sizeof(char));
    }


    if( ffgcno(fptr,1,"CAL_FILE",&filcol,status)) {
         sprintf(msg,"Unable to find CAL_FILE column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if( ffgcno(fptr,1,"CAL_XNO",&extcol,status)) {
         sprintf(msg,"Unable to find CAL_XNO column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if( ffgcno(fptr,1,"CAL_CBD",&cbdcol,status)) {
         sprintf(msg,"Unable to find CAL_CBD column");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
    }

    if(ffgcdw(fptr,cbdcol,&dw,status)) {
         sprintf(msg,"Unable to find CAL_CBD column display width");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
    }
    else {
          cbdval=(char*) calloc(dw+1,sizeof(char));
    }

    for ( i =0; i < nax2val; i++ ) {

         if( ffgcvs(fptr,telcol,i+1,1,1," ",&telval,&anynull,status) ) {
                sprintf(msg,"Error reading TELCOL column");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
         }
         
         if( ffgcvs(fptr,inscol,i+1,1,1," ",&insval,&anynull,status) ) {
                sprintf(msg,"Error reading INSCOL column");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
         }

         if(detnam[0] == '-' ) {
		strcpy(detval,"NONE");
	 }
	 else {
		if(ffgcvs(fptr,detcol,i+1,1,1," ",&detval,
                     &anynull,status) ) {
                        sprintf(msg,"Error reading DETCOL column");
                        HD_ERROR_THROW(msg,*status); 
                        goto cleanup;
		}
	 }

         if(strcmp(filt,"-") == 0 ) {
		strcpy(filtval,"NONE");
	 }
	 else {
		if(ffgcvs(fptr,filtcol,i+1,1,1," ",&filtval,
                     &anynull,status)) {
                        sprintf(msg,"Error reading FILCOL column");
                        HD_ERROR_THROW(msg,*status); 
                        goto cleanup;
		}
	 }


         if( ffgcvs(fptr,codcol,i+1,1,1," ",&codval,&anynull,status) ) {
                sprintf(msg,"Error reading CODCOL column");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
         }

         if( ffgcvd(fptr,refcol,i+1,1,1,0.0,&(refval[i]),&anynull,status) ) {
                sprintf(msg,"Error reading REFCOL column");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
         }

         if( ffgcvd(fptr,qulcol,i+1,1,1,0,&qulval,&anynull,status) ) {
                sprintf(msg,"Error reading QULCOL column");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
         }


        if(strcasecmp(telval,tele) == 0 &&
           strcasecmp(insval,instr) == 0 &&
          (strcasecmp(detval,detnam) == 0 ||
           strcasecmp(detval,"NONE") == 0 ||
           strcasecmp(detnam,"NONE") == 0 ) &&
          (strcasecmp(filtval,filt) == 0 ||
           strcasecmp(filtval,"NONE") == 0 ||
           strcasecmp(filt,"NONE") == 0) &&
          (strcasecmp(codval,codenam) == 0 ||
           strcasecmp(codval,"NONE") == 0 ||
           strcasecmp(codenam,"NONE") ==0 ) &&
           qulval == 0 ) index[i] =1;
        else index[i] =0;

        if(fits_get_coltype(fptr,cbdcol, &typecode,&repeat,&width,status)) {
                sprintf(msg,"Error reading CBDCOL column type");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
        }

        cbdnelem = (int) (repeat/width);


        for( k=1; k<=cbdnelem; k++) {
            if( ffgcvs(fptr,cbdcol,i+1,k,1," ",&cbdval,&anynull,status) ) {
                cbdval ="none";
                sprintf(msg,"Error reading CBDCOL column");
                HD_ERROR_THROW(msg,*status); 
                goto cleanup;
            }

	    if(list2) freecbd(list2);
            list2 = NULL;
            if(parseCBD(i,cbdval,&list2,status) ) {
	       sprintf(msg,"Parse CDBnXXX error");
               HD_ERROR_THROW(msg,*status); 
               goto cleanup;
             }

             index[i] =  index[i] && cmpCBD(list1,list2);
        }


        if( startreftime != -99) {
	        if( startreftime >= refval[i] ) {
                     if(refval[i] > maxrefval && index[i]) maxrefval =refval[i];
		     index[i] =index[i] && 1;
       		}
        	else {
	             index[i] =0;
        	}
        }
        if( stopreftime != -99) {
	        if( refval[i] <= stopreftime ) {

		     index[i] =index[i] && 1;
       		}
        	else {
	             index[i] =0;
        	}
        }



         
     }
     
     for( i=0; i<nax2val ; i++) {
         
         if( refval[i] < maxrefval ) index[i] =0;
     }


     *nfound =0;
     *nret =0;
     for( i=0; i<nax2val ; i++) {
        
       if(index[i] >0 ) {
         if( *nret < maxret ) {
         if( ffgcvs(fptr,dircol,i+1,1,1," ",&dirval,&anynull,status) ) {
                sprintf(msg,"Error reading DIRCOL column");
                HD_ERROR_THROW(msg,*status); 
		goto cleanup;
         }
         if( ffgcvs(fptr,filcol,i+1,1,1," ",&tmpfile,&anynull,status) ) {
                sprintf(msg,"Error reading FILCOL column");
                HD_ERROR_THROW(msg,*status); 
		goto cleanup;
         }

         if(cpthnm(caldbvar,dirval,tmpfile,status)) {
                sprintf(msg,"Fail to get path name");
                HD_ERROR_THROW(msg,*status); 
		goto cleanup;
         }


         strncpy(file[*nret],tmpfile,fnamesize - 1);
         file[*nret][fnamesize - 1] = 0;

         if( ffgcvj(fptr,extcol,i+1,1,1,0,extno+*nret,&anynull,status) ) {
                sprintf(msg,"Error reading EXTLCOL column");
                HD_ERROR_THROW(msg,*status); 
		goto cleanup;
         }
         if( ffgcvs(fptr,devcol,i+1,1,1," ",&online[*nret],&anynull,status) ) {
                sprintf(msg,"Error reading DEVCOL column");
                HD_ERROR_THROW(msg,*status); 
		goto cleanup;
         }
         *nret =*nret+1;
         }
         *nfound =*nfound+1;
       }
      }



cleanup:
    if(telval) free(telval);
    if(insval) free(insval);
    if(filtval) free(filtval);
    if(cbdval) free(cbdval);
    if(detval) free(detval);
    if(codval) free(codval);
    if(dirval) free(dirval);
    if(refval) free(refval);
    if(tmpfile) free(tmpfile);
    if(index) free(index);
     
    return *status;

}



/*
  HDgtcalf utility routines 
*/

/*-----------------------------------------------------------*/

static int gtcalidx (char* mode, char* missn,
              char* inst, char* ciffil, 
              char* instdir, int* status) 
/* Get CALDB index */
{

     const char* caldbvar="CALDB";
     const char* cnfgvar="CALDBCONFIG";

     char *caldb;
     char *config;

     char missval[160], instval[160];
     char cifdev[160], cifdir[160];
     char insdev[160], insdir[160];

     FILE * fptr;
     char msg[256];
     
     if(*status) return *status;  
  
     caldb = getenv(caldbvar);

     if( caldb == NULL || caldb[0] == '\0' ) {
         *status = HD_ERR_NULL_POINTER;
         sprintf(msg,"CALDB environment variable not set");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
     }

     config = getenv(cnfgvar);

     if( config == NULL || config[0] =='\0') {
         *status = HD_ERR_NULL_POINTER;
         sprintf(msg,"CALDBCONFIG environment variable not set");
         HD_ERROR_THROW(msg,*status); 
         goto cleanup;
     }

      fptr = fopen(config,"r");

      if( NULL == fptr ) {
         *status = HD_ERR_NULL_POINTER;
         sprintf(msg,"CALDBCONFIG environment variable not properly set");
         HD_ERROR_THROW(msg,*status);
         goto cleanup;
      }

      do {
      		if( rdcnfgl(fptr,missval,instval,cifdev,cifdir,ciffil,
       		     insdev,insdir,status) ) {
                        sprintf(msg,"Problem reading CALDB config file");
                        HD_ERROR_THROW(msg,*status); 
			goto cleanup;
                }
	        if( strcasecmp(missn,missval) == 0 ) {
		        if(strcasecmp(instval,inst) == 0 ) {
             			if(cpthnm(cifdev,cifdir,ciffil,status)) {
                                       sprintf(msg,"Problem getting cif");
                                       HD_ERROR_THROW(msg,*status); 
			               goto cleanup;
                                 }
                                strcpy(instdir,"");
             			if(cpthnm(insdev,insdir,instdir,status)) {
                                       sprintf(msg,"Problem getting instdir");
                                       HD_ERROR_THROW(msg,*status); 
			               goto cleanup;
             			}
				break;
			}
		}

        } while(1);

                


cleanup:
 if ( fptr ) fclose( fptr );
 return *status;

}

/*----------------------------------------------------------------------*/



static int rdcnfgl (FILE * fptr, char * missn, char * inst, 
             char* cifdev, char* cifdir, char* cif,
             char* datadev, char* datadir, int * status) 
/* Read config file */
{

       char line[1000];
       char * str, *str1, *str2, *str3, *str4, *str5, *str6;


       *status=CIF_ERR_READ_CONFIG;

       while (fgets(line,1000,fptr) != NULL ) {
	   if(line[0] !='\0' && line[0] !='#') {
                str = strtok(line," " );
                str1 = strtok(NULL," " );
                str2 = strtok(NULL," " );
                str3 = strtok(NULL," " );
                str4 = strtok(NULL," " );
                str5 = strtok(NULL," " );
                str6 = strtok(NULL," " );
                if( str != NULL &&
                    str1 != NULL &&
                    str2 != NULL &&
                    str3 != NULL &&
                    str4 != NULL &&
                    str5 != NULL &&
                    str6 != NULL ) 
                {

	            strcpy(missn,str);
                    strcpy(inst,str1);
                    strcpy(cifdev,str2);
                    strcpy(cifdir,str3);
                    strcpy(cif,str4);
                    strcpy(datadev,str5);
                    strcpy(datadir,str6);
                    *status = HD_OK;
                    break;
                 }
            }
        }
        return *status;

}
