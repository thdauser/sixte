/*
  HDgtcalf utility routines 
*/


#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fitsio.h"
#include "HDgtcalf_internal.h"

#ifdef HDGTCALF_STANDALONE
#include "HDgtcalf_standalone.h"
#else
#include "headas_utils.h"
#include "headas_error.h"
#endif

/* Function prototypes */
static int ccldj (int iy, int im, int id, double * djm, int *status);

static void copycbd(CBD * des, CBD * src); 
static CBD *  newcbd(); 
static CBDLIST * newcbdlist(); 
static char * strindex(char* str, int * j, char** ct, int n ) ;
static CBDLIST * addcbdlist(CBDLIST * list, CBD * cbd);
static int isNumber (char * s );
static void parseVal (char *str, double * val1, double * val2, int * type);


/*--------------------------------------------------------------------*/

/* parse string date to mjd */
int dt2mjd (const char* date, double * mjd, int * status) {

    int id,im,iy;
    char msg[256];

    if(*status) return *status; 


    if( strcasecmp(date,"NOW") == 0 ) {
         if(ffgsdt(&id,&im,&iy,status)) {
             sprintf(msg,"Error get system date");
             HD_ERROR_THROW(msg,*status); 
             return *status;
         }
    } 
    else {
         if(ffs2dt((char*) date,&iy,&im,&id,status)) {
             sprintf(msg,"Error parse date");
             HD_ERROR_THROW(msg,*status); 
             return *status;
         }
    }
            
    if(ccldj(iy,im,id,mjd,status)) {
             sprintf(msg,"Error convert date to mjd");
             HD_ERROR_THROW(msg,*status); 
             return *status;
    }

    return *status;

}

/* parse string time to doube */
int tim2df(const char * time, double * dayfrac, int * status) {

	int hour, min;
	double second;
        char msg[256];
	char time1[80];
        int timeref=0;

        if(strcasecmp(time,"NOW") == 0 ) {
           if(ffgstm(time1,&timeref,status)) {
               sprintf(msg,"Error get system date & time string");
               HD_ERROR_THROW(msg,*status); 
               return *status;
           }
           if (ffs2tm(time1,NULL,NULL,NULL,&hour,&min,&second,status)){
               sprintf(msg,"Error parse time");
               HD_ERROR_THROW(msg,*status); 
               return *status;
           }
        }
	else {
           if(ffs2tm((char*) time,NULL,NULL,NULL,&hour,&min,&second,status)) {
               sprintf(msg,"Error parse time");
               HD_ERROR_THROW(msg,*status); 
	       return *status;
           }
	}
	

	*dayfrac = (hour*3600+min*60 +second) /86400;

	return *status;
}

/* Convert Gregorian Calendar to Modified Julian Date */
static int ccldj (int iy, int im, int id, double * djm, int *status) {

    int mtab[] ={31,28,31,30,31,30,31,31,30,31,30,31};
    char msg[256];

    if(*status) return *status;

    if( iy < -4699 ) {
	 *status =CLDJ_ERR_BAD_YEAR;
         sprintf(msg,"Error convert date & time to mjd: bad year");
         HD_ERROR_THROW(msg,*status); 
	 return *status;
    }
    else if( im >= 1 && im <= 12 ) {
	if (iy%4 == 0 ) {
	   mtab[1] =29;
        }
        if ( iy%100== 0 && iy%400 !=0 ) {
	   mtab[1]=28;
        }

        if(id < 1 || id > mtab[im-1] ) {
	    *status =CLDJ_ERR_BAD_DAY;
            sprintf(msg,"Error convert date & time to mjd: bad day");
            HD_ERROR_THROW(msg,*status); 
	    return *status;
        }

        *djm = (1461*(iy-(12-im)/10+4712))/4;
        *djm +=(306*((im+9)%12)+5)/10;
        *djm -=(3*((iy-(12-im)/10+4900)/100))/4;
        *djm +=id-2399904;
    } else {
        *status =CLDJ_ERR_BAD_MONTH;
         sprintf(msg,"Error convert date & time to mjd: bad month");
         HD_ERROR_THROW(msg,*status); 
         return *status;
    }

    return *status;

}

/* Get path name */
int cpthnm(const char* disk, const char* dir, char* file, 
             int* status)
{
	char* root;
        char path[FILENAME_MAX];
        char msg[256];

        if(*status) return *status;

        root = getenv(disk);

        if( root == NULL) {
                if( strcmp(disk,"") !=0 ) {
                    if(disk[0] == '/') {
                       strcpy(path,disk);
                    } else {
        	       strcpy(path,"/");
        	       strcat(path,disk);
                    }
                }
                else {
                       strcpy(path,""); 
                }
        }
        else {
		strcpy(path,root);
        }

        if( dir !=NULL ) {
             if( strcmp(dir,"") !=0 ) {
	             if(dir[0] =='/') {
	                 strcat(path,dir);
       		      } 
	             else {
			 strcat(path,"/");
	   	         strcat(path,dir);
             	      }
             }
         }
         else {
             *status =HD_ERR_NULL_POINTER;
             sprintf(msg,"Error get path name: dir is NULL string");
             HD_ERROR_THROW(msg,*status); 
         } 
             

  
 

	if( file != NULL ) {
                if( strcmp(file,"") != 0) {
                   if(file[0] =='/') {
		      strcat(path,file);
                   }
                   else {
		      strcat(path,"/");
		      strcat(path,file);
                   }
                       
                }
	}
        else {
            *status =HD_ERR_NULL_POINTER;
            sprintf(msg,"Error get file name: file is NULL string");
            HD_ERROR_THROW(msg,*status); 
        } 
        strcpy(file,path);

	return *status;
}


/* CBDxxx utilites routines */ 

void freecbd(CBDLIST * list) {
       CBDLIST * list1=NULL;
	  
	while (list != NULL ) {
           if(list->cbd) free (list->cbd);
	    list1 =list;
	    list = list->next;
        }

        while (list1 !=NULL) {
             list = list1->prev;
             free(list1);
	     list1 =list;
        }
	
	

}
	

int parseCBD2(const char* str, CBDLIST ** cbdlist, int* status) {

    CBD *cbd;

    int i, j;
    char * s1;
    char * s2;
    char * s3;
    int last = 0;


    int ops1[2] ={CBD_AND,CBD_OR};
    int ops2[5] ={CBD_LT,CBD_LE,CBD_EQ,CBD_GE,CBD_GT};

    int ops3;

    char * op1[2] ={".AND.",".OR."};
    char * op2[5] ={".LT.",".LE.",".EQ.",".GE.",".GT."};
     

    if(*status) return *status;

    if(str == NULL || str[0]=='\0' ) {
        *cbdlist =NULL;
	return *status;
    }


    s1 = (char*) str;
    ops3 = CBD_AND;

    while (!last) {
          s2 =strindex(s1,&i,op1,2);
	  
          if (s2 == NULL ) {
		last =1; 
                s3 = strindex(s1,&j,op2,5);
                if(s3 == NULL ) break;

                cbd = newcbd();
		strncpy(cbd->name,s1,s3-s1+1);
                cbd->name[s3-s1] ='\0';
		strcpy(cbd->sval,s3+4);
                parseVal(cbd->sval,
                      &(cbd->minval),
                      &(cbd->maxval),
                      &(cbd->type));
                cbd->op2 = ops2[j];
                cbd->op1 = ops3;
                *cbdlist = addcbdlist(*cbdlist,cbd);
           }
           else {
                s3 = strindex(s1,&j,op2,5);
		if(s3 == NULL ) break;

                cbd = newcbd();
		strncpy(cbd->name,s1,s3-s1+1);
                cbd->name[s3-s1] ='\0';
		strncpy(cbd->sval,s3+4,s2-s3-3);
                cbd->sval[s2-s3-4]='\0';
                parseVal(cbd->sval,
                      &(cbd->minval),
                      &(cbd->maxval),
                      &(cbd->type));
                cbd->op2 = ops2[j];
                cbd->op1 = ops3;
                ops3 = ops1[i];  
                if( i== 0) {
			s1 =s2+5;
                }
                else {
			s1=s2+4;
		}
                *cbdlist = addcbdlist(*cbdlist,cbd);
	}

     }

     return *status;

}

int parseCBD(int m, const char* str, CBDLIST ** cbdlist, int* status) {

    int i ;
    int j ;
    int k ;
    int len;

    CBD cbd;
    CBDLIST * list=NULL;

    if(*status) return *status;

    if (str == NULL || strstr(str,"NONE")==str) goto cleanup;

    if( *cbdlist !=NULL ) 
	    list = *cbdlist;


    i=0;
    len =strlen(str);

    while (i < len )  {
	if( str[i] == '(' ) {
             strncpy(cbd.name,str,i+1);
	     cbd.name[i] ='\0';
             *status =CBD_OK;
             break;
        }
	*status =CBD_ERR_PARSE;
        i++;
    }


    if(*status ) goto cleanup;

    j=i+1;

    while ( j < len)  {
	if( str[j] ==')') {
		strncpy(cbd.units,str+j+1,len-j );
                cbd.units[len-j-1]='\0';
                *status =CBD_OK;
                break;
        }
	*status =CBD_ERR_PARSE;
        j++;
       
    } 

   if(*status ) goto cleanup;

   

   for(k=i+1; k<=j; k++) {
       if( str[k] ==',' || str[k]==')') {
	  strncpy(cbd.sval,str+i+1,k-i);
          cbd.sval[k-i-1]='\0';
          parseVal(cbd.sval,
                      &(cbd.minval),
                      &(cbd.maxval),
                      &(cbd.type));
          i=k;
          if(cbd.sval[0] !='\0' ) {
              if(list == NULL ) {
                 list = newcbdlist();
                 list->cbd = newcbd();
                 copycbd(list->cbd, &cbd);
                 *cbdlist = list;
              }
              else {
		 list->next = newcbdlist();
                 list->next->prev = list;
                 list->next->cbd = newcbd();
                 copycbd(list->next->cbd,&cbd);
                 list = list->next;
              }
          }
       }
   }


cleanup:
	return *status; 
          
}

static void parseVal (char *str, double * val1, double * val2, int * type) { 
      
      char * s;
      char * s1;

      *val1 =0.;
      *val2 =0.;

      if( str == NULL || str[0] =='\0' ) return;

      if(str[0]=='"') {
	*type =CBD_STR;
	return;
      }

      s=strtok(str,"-");
      s1=strtok(NULL,"-");
	
      if(s1 == NULL) {
         if(isNumber(s)) {
		*val1 =atof(s);
		*val2 =atof(s);
                *type =CBD_VAL;
         }
         else {
		*type =CBD_STR;
	}
     }
     else if ( isNumber(s) && isNumber(s1) ) {
		*val1 =atof(s);
                *val2 =atof(s1);
                *type =CBD_RANGE;
     } else {
	*type =CBD_STR;
     }
}

static void copycbd(CBD * des, CBD * src) {
     strcpy(des->name,src->name);
     strcpy(des->units,src->units);
     strcpy(des->sval,src->sval);
     des->type = src->type;
     des->minval = src->minval;
     des->maxval = src->maxval;
     des->op1 = src->op1;
     des->op2 = src->op2;
}

static CBDLIST * newcbdlist() {
	CBDLIST * list;
	list = (CBDLIST*) calloc(1,sizeof(CBDLIST));
        list->cbd =NULL;
        list->prev =NULL;
        list->next =NULL;
	return list;
}
    
static CBD *  newcbd() {
       CBD * cbd;
	cbd= (CBD*) calloc(1,sizeof(CBD));

        cbd->minval =0.;
        cbd->maxval =0.;
	return cbd;
       
        
}

static CBDLIST * addcbdlist(CBDLIST * list, CBD *cbd) {
        CBDLIST * list1;

	if(cbd == NULL ) return list;

         if(list == NULL ) {
		list = newcbdlist();
		list->cbd = newcbd();
		copycbd(list->cbd,cbd); 
         }
         else {
                list1=list;
		while(list->next !=NULL) list=list->next;
		list->next = newcbdlist();
		list->next->cbd = newcbd();
		copycbd(list->next->cbd,cbd); 
		list->next->prev =list;
		list =list1;
         }

         return  list;
}


static char * hd_strcasestr(char* s1, const char* s2 ) {
	int i=0;
	int l1;
	int l2;
        int r;

	l1 = strlen(s1);
	l2 = strlen(s2);

      
        i=0;
	r =-1;
        while (i<l1-l2) {
	    if(strncasecmp(s1+i,s2,l2) == 0) {
		 r = i;
            	 break;
            }
            i++;
        }

        if( r == -1 ) {
           return NULL;
        }
	else {
           return &s1[i];
        }
}
		 
		
		
static char * strindex(char* str, int * j, char** ct, int  n ) {
    
     char* s1;
     char* s2;
     int i;
     
     *j =0;

     if( str == NULL ) return NULL;


     s1 = hd_strcasestr(str,ct[*j]);

     for (i=1; i<n; i++) {
         s2 = hd_strcasestr(str,ct[i]);
         if(s1 == NULL ) { 
 	     s1 = s2;
             *j =i;
         }
         else if( s2 != NULL && s1 > s2) {
	     s1 = s2;
             *j = i;
	 }
     }

     
     return s1;

}

int cmpCBD(CBDLIST * list1, CBDLIST * list2) {
   CBDLIST * list3;
   int result =1;
   int r;

   if( list1 == NULL || list2 == NULL ) return 1;

   result  =1;
   list3 =list2;
   r =1;
   while ( list1 != NULL ) {
      list2 =list3;
      while (list2 !=NULL ) {
	if(strcasecmp(list1->cbd->name,list2->cbd->name) != 0 ) {
	 	r =1;
	}
	else if( list1->cbd->type == CBD_STR && list2->cbd->type == CBD_STR) {
	       if(strcasecmp(list1->cbd->sval,list2->cbd->sval)== 0) r=1;
	       else r=0;
        }
        else if( list1->cbd->type == CBD_VAL && list2->cbd->type
                  == CBD_VAL ) {
             switch (list1->cbd->op2) {
		case CBD_LT:
                    if( list2->cbd->minval < list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_LE:
                    if( list2->cbd->minval <= list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_EQ:
                    if( list2->cbd->minval == list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_GE:
                    if( list2->cbd->minval >= list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_GT:
                    if( list2->cbd->minval > list1->cbd->minval) r=1;
                    else r=0;
                    break;
              }
        }
	else if(list1->cbd->type == CBD_VAL && list2->cbd->type 
                == CBD_RANGE) {
	    switch(list1->cbd->op2) {
		case CBD_LT:
                    if( list2->cbd->minval < list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_LE:
                    if( list2->cbd->minval <= list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_EQ:
                    if( list2->cbd->minval <= list1->cbd->minval &&
                        list2->cbd->maxval >= list1->cbd->minval ) r=1;
                    else r=0;
                    break;
		case CBD_GE:
                    if( list2->cbd->maxval >= list1->cbd->minval) r=1;
                    else r=0;
                    break;
		case CBD_GT:
                    if( list2->cbd->maxval > list1->cbd->minval) r=1;
                    else r=0;
		    break;
                default:
                   break;
             }
         }
	/*list2 =list2->next;*/
        list2 = (r==0) ? list2->next : NULL;
     }

   if( list1->cbd->op1 == CBD_OR ) {
	result = result || r;
   }
   else if(list1->cbd->op1 == CBD_AND)  {
	result = result && r;
   } 
	list1=list1->next;
        r =0;

   }        
       
	
   return result;
}



static int isNumber (char * s ) {

    int i;
    int result =1;
    

    if( s == NULL || s[0] =='\0' ) {
	result =0;
        return result;
    }

    for(i=0; i<strlen(s); i++) {
	if( isdigit(s[i]) || s[i] =='.') {
		result = 1;
        }
        else {
	       result = 0;
	       return result;
        }
    }

    return result;

}


