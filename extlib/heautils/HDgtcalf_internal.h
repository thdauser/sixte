/* Enumerated error codes for HDgtcalf errors */

typedef enum hdgtcalf_error_code {

   CLDJ_ERR_BASE = 2000,
   CLDJ_ERR_BAD_YEAR,
   CLDJ_ERR_BAD_MONTH,
   CLDJ_ERR_BAD_DAY,

   CIF_ERR_BASE = 2050,
   CIF_ERR_EMPTY_TABLE,
   CIF_ERR_READ_CONFIG,

   CBD_OK  =0,
   CBD_ERR_BASE = 2100,
   CBD_ERR_PARSE

} hdgtcalf_error_code;


/* define CBDxxx  structure */

#define CBD_AND 1
#define CBD_OR 2
#define CBD_LT 3
#define CBD_LE 4
#define CBD_EQ 5
#define CBD_GE 6
#define CBD_GT 7
#define CBD_VAL 8
#define CBD_STR 9
#define CBD_RANGE 10

typedef struct cbd{
     char name[20];
     int type;
     double minval;
     double maxval;
     char sval[80];
     char units[20];
     int op1;
     int op2;

} CBD;

typedef struct cbdlist {
     CBD * cbd;
     struct cbdlist * prev;
     struct cbdlist * next;
} CBDLIST;



/* prototype functions */


int cpthnm(const char* disk, const char* dir, char* file, int* status);

int dt2mjd (const char* date, double * mjd, int * status);

int tim2df(const char * time, double * dayfrac, int * status);



int parseCBD(int i, const char* str, CBDLIST ** cbdlist, int* status);

int parseCBD2(const char* str, CBDLIST ** cbdlist, int* status);

int cmpCBD( CBDLIST * list1, CBDLIST * list2 );

void freecbd(CBDLIST * list);



