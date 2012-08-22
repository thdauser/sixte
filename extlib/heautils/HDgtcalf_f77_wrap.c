/******************************************************************************
 *   File name: HDgtcalf_f77_wrap.c                                           *
 *                                                                            *
 * Description: Cfortran wrapper for HDgtcalf.c.                              *
 *                                                                            *
 *      Author: Ziqin Pan, L3 Communications, for HEASARC/GSFC/NASA           *
 *              James Peachey, L3 Communications, for HEASARC/GSFC/NASA       *
 *                                                                            *
 *****************************************************************************/
#include "hdcal.h"
#include "fitsio.h"
#include "f77_wrap.h"

/*Fortran wrapper */
#define HDgtcalf_STRV_A13 NUM_ELEMS(velem)
#define HDgtcalf_STRV_A15 NUM_ELEMS(velem)
CFextern VOID_cfF(HDGTCALF,hdgtcalf)
CFARGT27(NCF,DCF,ABSOFT_cf2(VOID),STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,PINT,PINT,PSTRINGV,PLONG,PSTRINGV,PINT,PINT,PINT,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0));
CFextern VOID_cfF(HDGTCALF,hdgtcalf)
CFARGT27(NCF,DCF,ABSOFT_cf2(VOID),STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,PINT,PINT,PSTRINGV,PLONG,PSTRINGV,PINT,PINT,PINT,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0,CF_0))
{
    QCF(STRING,1)
    QCF(STRING,2)
    QCF(STRING,3)
    QCF(STRING,4)
    QCF(STRING,5)
    QCF(STRING,6)
    QCF(STRING,7)
    QCF(STRING,8)
    QCF(STRING,9)
    QCF(STRING,10)
    QCF(PINT,11)
    QCF(PINT,12)
    QCF(PSTRINGV,13)
    QCF(PLONG,14)
    QCF(PSTRINGV,15)
    QCF(PINT,16)
    QCF(PINT,17)
    QCF(PINT,18)

    int velem;
    char* tele, *instr;
    char* detnam, *filt;
    char* codenam, *strtdate;
    char* strtime, *stpdate, *stptime;
    char* expr;
    int maxret;
    int fnamesize;
    char** filenam;
    long* extno;
    char** online;
    int* nret; 
    int* nfound;
    int* status;

    tele = TCF(HDgtcalf,STRING,1,0);
    instr = TCF(HDgtcalf,STRING,2,0);
    detnam = TCF(HDgtcalf,STRING,3,0);
    filt = TCF(HDgtcalf,STRING,4,0);
    codenam = TCF(HDgtcalf,STRING,5,0);
    strtdate = TCF(HDgtcalf,STRING,6,0);
    strtime = TCF(HDgtcalf,STRING,7,0);
    stpdate = TCF(HDgtcalf,STRING,8,0);
    stptime = TCF(HDgtcalf,STRING,9,0);
    expr = TCF(HDgtcalf,STRING,10,0);
    maxret = TCF(HDgtcalf,INT,11,0);
    velem = maxret;

    fnamesize = TCF(HDgtcalf,INT,12,0);
    filenam = TCF(HDgtcalf,PSTRINGV,13,0);
    extno = TCF(HDgtcalf,PLONG,14,0);
    online = TCF(HDgtcalf,PSTRINGV,15,0);
    nret = TCF(HDgtcalf,PINT,16,0);
    nfound = TCF(HDgtcalf,PINT,17,0);
    status = TCF(HDgtcalf,PINT,18,0);

    HDgtcalf(tele, instr,
           detnam,  filt,
           codenam, strtdate,
           strtime, stpdate, stptime,
           expr, maxret, fnamesize,filenam,
	   extno, online, nret, 
           nfound, status);

    RCF(STRING,1);
    RCF(STRING,2);
    RCF(STRING,3);
    RCF(STRING,4);
    RCF(STRING,5);
    RCF(STRING,6);
    RCF(STRING,7);
    RCF(STRING,8);
    RCF(STRING,9);
    RCF(STRING,10);
    RCF(PINT,11);
    RCF(PINT,12);
    RCF(PSTRINGV,13);
    RCF(PLONG,14);
    RCF(PSTRINGV,15);
    RCF(PINT,16);
    RCF(PINT,17);
    RCF(PINT,18);
}
