#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pil.h"
#include "headas_utils.h"
#include "headas_error.h"
#include "f77_wrap.h"

int headas_clobpar = 0;

/********************************** hdbasename **********************/
/* Returns the right-most part of the argument path which does not
   contain '/'. (May be empty.) */
const char *hdbasename(const char *path) {
  const char *p = path;
  const char *basename = path;
  while(p && *p) if('/' == *p++) basename = p;
  return basename;
}

/********************************** headas_clobber **********************/

int headas_clobberfile(char *filename) {

  /* 
     delete (clobber) the specified file, IF the file exists and IF
     the 'clobber' parameter (in the task .par file) is equal to 'yes'.
  */
  if (headas_clobpar)
    remove(filename);
  
  return(0);
}
FCALLSCFUN1(INT, headas_clobberfile, HDCLOBBER, hdclobber, STRING)
     
/********************************** hd_ran2 **********************/

/* Random number generator based on ran2() from Numerical Recipes in C, 
 * 2nd ed., p282.
 * 
 * Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
 * the endpoint values). Call with idum a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a sequence.
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float hd_ran2(long *idum){
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0){
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
FCALLSCFUN1(FLOAT, hd_ran2, HD_RAN2, hd_ran2, PLONG)

/* ----------------------------------------------------------------- */
/* 
 * expand_item_list - Expand a delimited list, or an "@" file
 *
 * Parses a delimited list or "@" file into individual strings.  An
 * array of string pointers is returned, one pointer for each item in
 * the list.
 *
 * A delimited list is broken into one string per delimited element.
 * The user can choose the delimiter character, often a comma.  An "@"
 * file is broken into one string per line of the file.
 * 
 * The user can select whether leading and trailing blanks are trimmed
 * or whether empty items are removed.  Also, for a delimited list,
 * the user can guard against delimiters embedded within CFITSIO
 * vector syntax, as in the comma within "item1[x,y]".
 *
 * If the return value of this function is non-zero, then the user is
 * responsible for free()'ing it, but the individual list elements
 * should not be free()'d.
 * 
 * If an expression such as "@filename.txt" is passed, and
 * filename.txt is not readable, then the expression is parsed as a
 * delimited list instead.
 *
 * char *liststr - string containing list expression, either
 *                 "@filename" - indicating a file whose contents are parsed,
 *                               one item per line in the file;
 *                 "expression..." - a delimited list expression to be parsed
 * int *nitems - upon return, *nitems contains number of items found.
 *               Note that (*nitems == 0) is a possible valid result,
 *               if the list is empty.
 * char fieldsep - delimiter character for delimited lists, e.g. ','
 *                 ignored for "@" files
 * int trim - Trim leading and trailing blanks? (1=yes, 0=no)
 * int skipempty - Remove empty strings from list? (1=yes, 0=no)
 * int guardparen - Guard against commas embedded in parens? (1=yes, 0=no)
 *                  also guards against quoted strings;  either type
 *                  of quotation mark is allowed, and a matching pair
 *                  is required; nested quotations are not allowed.
 *
 * int *status - error status; if non-zero, indicates failure to
 *               allocate memory or read an "@" file.
 *
 * RETURNS: 0 upon failure, or empty list.
 *     Otherwise, a pointer to an array of string pointers.
 *     User responsible to free() this pointer when done.
 *
 * EXAMPLES:
 *   result = expand_item_list("@filename.txt",0,*nitems,1,1,1, &status);
 *      - parsed according to contents of filename.txt
 *   result = expand_item_list("item1,item2",',',*nitems,1,1,1, &status);
 *      - parses to "item1","item2"
 *   result = expand_item_list("item1 ,item2",',',*nitems,1,1,1, &status);
 *      - parses to "item1","item2" (note trailing space)
 *   result = expand_item_list("item1,,item2",',',*nitems,1,1,1, &status);
 *      - parses to "item1","item2" (note skipped empty field)
 *   result = expand_item_list("item1[x,y],item2",',',*nitems,1,1,1, &status);
 *      - parses to "item1[x,y]","item2" (note CFITSIO syntax obeyed)
 * */
char **expand_item_list(char *liststr, int *nitems, char fieldsep,
			int trim, int skipempty, int guardparen,
			int *status)
{
  int allocLen, totalLen, llen;
  char *lines,line[PIL_PATH_MAX];
  char *ptr, *pdest, **cptr;
  int i;
  FILE *aFile;
  
  *nitems = 0;
  if (*status) return 0;
  
  /* Initial error checking */
  if ((liststr == 0) || (liststr[0] == 0)) return 0;
  if (nitems == 0) return 0;
  if (fieldsep == 0) fieldsep = ',';
  
  /* *************
     CASE 1: an "@" file containing a list of items, one per line */
  if (liststr[0] == '@') {
    
    /* Cribbed from ffimport_list() in cfitsio/cfileio.c */
    totalLen =    0;
    allocLen = 1024;
    lines    = (char *)malloc( allocLen * sizeof(char) );
    if( !lines ) {
      fprintf(stderr,"Couldn't allocate memory to hold ASCII file contents.");
      *status = MEMORY_ALLOCATION;
      return 0;
    }
    lines[0] = '\0';
    
    if( (aFile = fopen( liststr+1, "r" )) == NULL ) {
      /* Couldn't open @file, so just parse it explicitly */
      goto DELIMITED_LIST;
    }
    
    while( fgets(line,PIL_PATH_MAX,aFile)!=NULL ) {
      /* Make sure null terminated string */
      line[PIL_PATH_MAX-1] = 0;
      llen = strlen(line);
      
      /* Skip blank lines */
      if (skipempty && (llen == 0)) continue;
      
      /* Trim trailing spaces, tabs, newlines */
      ptr = line+llen-1;
      while ( (ptr >= line) && 
	      ((trim && ((*ptr == ' ')  || (*ptr == '\t')))
	       || (*ptr == '\n') || (*ptr == '\r')) ) {
	*ptr-- = 0;
	llen--;
      }
      if (skipempty && (llen == 0)) continue;
      
      ptr = line;
      /* Trim leading spaces and tabs */
      if (trim) {
	while ((*ptr == ' ') || (*ptr == '\t')) {
	  ptr++; 
	  llen--;
	}
	if (skipempty && (llen == 0)) continue;
      }
      
      if (totalLen + llen + 3 >= allocLen) {
	allocLen += 256;
	lines = (char *)realloc(lines, allocLen * sizeof(char) );
	if( ! lines ) {
	  fprintf(stderr, 
		  "Couldn't re-allocate memory to hold ASCII file contents.");
	  *status = MEMORY_ALLOCATION;
	  break;
	}
      }
      
      /* Copy string to output, be sure to add 1 to length to include null */
      strcpy( lines+totalLen, ptr );
      totalLen += llen+1;
      
      (*nitems) ++;
    }
    fclose(aFile);
    
    if (*nitems == 0) {
      free(lines);
      return 0;
    }
    
    lines = (char *)realloc(lines, totalLen*sizeof(char) 
			    + (*nitems)*sizeof(char *));
    if (lines == 0) {
      fprintf(stderr, 
	      "Couldn't re-allocate memory to hold pointer array.");
      *status = MEMORY_ALLOCATION;
      return 0;
    }

    /* Move char data to end of array, starting from the end */
    pdest = (char *) ( ((char **) lines) + (*nitems) ) + totalLen - 1;
    ptr   = lines + totalLen - 1;
    while(ptr >= lines) {
      *pdest-- = *ptr--;
    }

    cptr = (char **) lines;
    ptr = (char *) ( ((char **) lines) + (*nitems) );
    for (i=0; i< (*nitems); i++) {
      cptr[i] = ptr;       /* Put pointer to string i */
      while (*ptr) ptr++;  /* Skip chars of string */
      ptr++;               /* Skip null of string */
    }

  } else {
  /* *************
     CASE 2: an explicit delimited string */
  DELIMITED_LIST:

    /* Check for empty string */
    llen = strlen(liststr);
    if (llen == 0) return 0;

    /* Scan for number of commas */
    *nitems = 1;
    for(ptr=liststr; *ptr; ptr++) {
      if (*ptr == fieldsep) {
	(*nitems) ++;
      }
    }

    /* Allocate enough memory to store pointer array and strings */
    allocLen = (llen+1)*sizeof(char) + (*nitems)*sizeof(char *);
    lines = (char *)malloc(allocLen);
    if( ! lines ) {
      fprintf(stderr, 
	      "Couldn't allocate memory to hold file list.");
      *status = MEMORY_ALLOCATION;
      return 0;
    }
    memset(lines, 0, allocLen);   /* ... zero it */

    cptr = (char **) lines;
    pdest = (char *) ( ((char **) lines) + (*nitems) );

    ptr = liststr;
    do {
      /* Skip leading spaces */
      if (trim) {
	while ((*ptr) && (*ptr == ' ')) ptr++;
      }
      /* Copy up to next comma */
      if ((*ptr != fieldsep) || (skipempty == 0)) {
	int depth = 0, inquote = 0;

	*cptr = pdest;
	/* Handle case of vector columns with comma inside expression */
	while ((*ptr) && ((*ptr != fieldsep) || (depth > 0) || inquote)) {
	  if (guardparen) {

	    /* Check for balanced parentheses (outside of quoted strings) */
	    if (! inquote ) {
	      if ((*ptr == '[') || (*ptr == '(') || (*ptr == '{')) depth ++;
	      if ((*ptr == '}') || (*ptr == ')') || (*ptr == ']')) depth --;
	    }
	    
	    /* Check for balanced quotation marks */
	    if (! inquote && ((*ptr == '\'') || (*ptr == '"'))) {
	      inquote = *ptr;              /* Save entry quote mark ... */
	    } else if (*ptr == inquote) {  /* ... to compare with exit  */
	      inquote = 0;
	    }
	  }
	  *pdest++ = *ptr++;
	}
	/* Trim any trailing blanks */
	if (trim) {
	  if (pdest > *cptr) pdest--;
	  while ((pdest >= *cptr) && (*pdest == ' ')) pdest--;
	  if ((*pdest) && (*pdest != ' ')) pdest++;
	}

	/* Check for an empty string */
	if (skipempty && (pdest == *cptr)) continue;

	/* Null terminate dest string */
	*pdest++ = 0;
	cptr++;
      }
      if (*ptr == fieldsep) ptr++;
    } while (*ptr);

    /* Check for no entries */
    *nitems = cptr - (char **) lines;
    if (cptr == (char **) lines) {
      *nitems = 0;
      free(lines);
      return 0;
    }
  }

  /* Success, return pointer array (and char data that follows) */
  return (char **) lines;
}
