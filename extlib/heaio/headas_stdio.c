#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cfortran.h"
#include "headas_stdio.h"
#include "fitsio2.h"

#ifdef printf
#undef printf
#endif

#ifdef WIN32
int strcasecmp(const char * p1, const char * p2) {
  /* TODO Implement this. */
  return 0;
}
#endif

#ifdef __cplusplus
  typedef std::iostream hd_cxx_std_iostream;
#else
  typedef void hd_cxx_std_iostream;
#endif

struct hd_FILE_s {
  FILE **cStream;
  hd_cxx_std_iostream **cxxStream;
};

/* Define global variables related to I/O. */
FILE *heaerr = 0;
FILE *heain = 0;
FILE *heaout = 0;
FILE *heaprom = 0;
FILE *heanullstream = 0;

static hd_FILE_t sHDerrStream = { &heaerr, 0 };
static hd_FILE_t sHDinStream = { &heain, 0 };
static hd_FILE_t sHDoutStream = { &heaout, 0 };
static hd_FILE_t sHDpromptStream = { &heaprom, 0 };
static hd_FILE_t sHDnullStream = { &heanullstream, 0 };

hd_FILE *hd_err = &sHDerrStream;
hd_FILE *hd_in = &sHDinStream;
hd_FILE *hd_out = &sHDoutStream;
hd_FILE *hd_prompt = &sHDpromptStream;
hd_FILE *hd_null = &sHDnullStream;

int headas_chatpar = 0;

int HDIO_init(void){
    int status=0;
    const char *hderror = 0;
    const char *hdoutput = 0;
    const char *hdprompt = 0;

    heanullstream = tmpfile();

    /* Make sure streams have some kind of sensible default values. */
    heaerr = stderr;
    heain = stdin;
    heaout = stdout;
    heaprom = heaerr;

    /* Attempt to redirect stderr first if needed. This is
       the most important redirection. */
    if(0 == status) {
      /* Redirect stderr if HEADASERROR variable is set to a file name. */
      hderror = HD_getenv("HEADASERROR");
      if(hderror) {
        /* Special cases null, stderr and stdout handled below so do
           nothing for the moment. */
        if(strcasecmp("null",hderror) &&
            strcasecmp("stderr",hderror) &&
            strcasecmp("stdout",hderror)) {
          heaerr = freopen(hderror,"a+",stderr);
          if(NULL == heaerr) {
            /* Write error message. */
            printf("Unable to redirect stderr to the location given by "
                "the environment variable:\n\tHEADASERROR == %s (at %s:%d)\n",
                hderror, __FILE__, __LINE__);
            
            status = errno;

            /* Reopen failed for stderr. The error must be reported
               *somehow*, so redirect heaerr to point to stdout. */
            heaerr = stdout;

            /* Unbuffer stdout, like stderr. */
            setbuf(stdout, 0);
          }
        }
      }
    }

    if(0 == status) {
      /* Redirect stdout if HEADASOUTPUT variable is set to a file name. */
      hdoutput = HD_getenv("HEADASOUTPUT");
      if(hdoutput) {
        /* Special cases null, stderr and stdout handled below so do
           nothing for the moment. */
        if(strcasecmp("null",hdoutput) &&
            strcasecmp("stderr",hdoutput) &&
            strcasecmp("stdout",hdoutput)) {
          heaout = fopen(hdoutput,"a+");
          if(NULL == heaout) {
            /* Write error message. */
            printf("Unable to redirect stdout to the location given by "
                "the environment variable:\n\tHEADASOUTPUT == %s (at %s:%d)\n",
                hdoutput, __FILE__, __LINE__);

            status = errno;
            heaout = heanullstream;
          }
        }
      }
    }

    if(0 == status) {
      if (HD_getenv("HEADASNOQUERY")) heaprom = heanullstream;
      else {
        /* Redirect prompts if HEADASPROMPT variable is set to a file name. */
        hdprompt = HD_getenv("HEADASPROMPT");
        if(hdprompt) {
          /* Special cases null, stderr and stdout handled below so do
             nothing for the moment. */
          if(strcasecmp("null",hdprompt) &&
              strcasecmp("stderr",hdprompt) &&
              strcasecmp("stdout",hdprompt)) {
            heaprom = fopen(hdprompt,"a+");
            if(NULL == heaprom) {
              /* Write error message. */
              printf("Unable to redirect prompts to the location given by "
                "the environment variable:\n\tHEADASPROMPT == %s (at %s:%d)\n",
                hdprompt, __FILE__, __LINE__);

              status = errno;
              heaprom = heanullstream;
            }
          }
        } else {
          heaprom = fopen("/dev/tty","a+");
          if(NULL == heaprom) {
            status = errno;
            heaprom = heanullstream;

            /* Write error message. */
            printf("Unable to redirect prompts to the /dev/tty (at %s:%d)\n",
                __FILE__, __LINE__);
          }
        }
      }
    }

    if (0 == status) {
      /* Handle symbolic redirections, which were skipped above so
         that they would have their intended effects. */
      if(hderror) {
        if(0 == strcasecmp("null", hderror)) {
          heaerr = heanullstream;
          hd_err->cStream = &heanullstream;
        } else if(0 == strcasecmp("stdout",hderror)) {
          heaerr = heaout;
          hd_err->cStream = &heaout;
        }
      }

      if(hdoutput) {
        if(0 == strcasecmp("null", hdoutput)) {
          heaout = heanullstream;
          hd_out->cStream = &heanullstream;
        } else if(0 == strcasecmp("stderr",hdoutput)) {
          heaout = heaerr;
          hd_out->cStream = &heaerr;
        }
      }

      if(hdprompt) {
        if(0 == strcasecmp("null", hdprompt)) {
          heaprom = heanullstream;
          hd_prompt->cStream = &heanullstream;
        } else if(0 == strcasecmp("stderr",hdprompt)) {
          heaprom = heaerr;
          hd_prompt->cStream = &heaerr;
        } else if(0 == strcasecmp("stdout",hdprompt)) {
          heaprom = heaout;
          hd_prompt->cStream = &heaout;
        }
      }
    }

    return status;
}

/* This is a C-callable replacement for printf which
   also takes a chatter threshold argument */
int headas_chat(int threshold, const char* fmt, ...){
    va_list ap;
    int status=0;

    if (headas_chatpar < 0){
	fprintf(stderr,"Error: headas_chat(): no chatter parameter present\n");
	status = 1;
	return status;
    }

    if (headas_chatpar >= threshold){
	va_start(ap, fmt);
	status = vfprintf(heaout, fmt, ap);
	va_end(ap);
    }

    return status;
}


/* This can be #define'd to replace printf or
   called directly by a C application */
int headas_printf(const char* fmt, ...) {
  va_list ap;
  int status=1;

  va_start(ap, fmt);
  status = vfprintf(heaout, fmt, ap);
  va_end(ap);

  return status;
}

/* This is intended to be a printf replacement
   substituted when compiling PIL (sets up a
   dedicated stream via HEADASPROMPT environment var) */
int pil_printf(const char* fmt, ...) {
  va_list ap;
  int status=1;

  va_start(ap, fmt);
  status = vfprintf(heaprom, fmt, ap);
  va_end(ap);

  return status;
}

static const size_t sEnvVarArrayStep = 128u;

typedef struct hd_env_var {
  char *Name;
  char *Value;
} hd_env_var;

typedef struct hd_env_var_array {
  hd_env_var *Begin;
  hd_env_var *End;
} hd_env_var_array;

static hd_env_var_array sEnvVars = { 0, 0};

static char *CpyStr(const char *s) {
  size_t l;
  char *r;

  if (s) {
    l = strlen(s);
    r = (char *) malloc((l + 1) * sizeof(char));
    if (r) strcpy(r, s);
  } else {
    r = 0;
  }

  return r;
}

static hd_env_var *FindEnvVar(const char *name) {
  hd_env_var *v_p;

  if (!sEnvVars.Begin) {
    /* First time called; initialize array of environment variables. */
    sEnvVars.Begin = (hd_env_var *) malloc(
        sEnvVarArrayStep * sizeof(hd_env_var)
    );

    /* Initialize pointers. */
    if (sEnvVars.Begin) {
      sEnvVars.End = sEnvVars.Begin + sEnvVarArrayStep;
      v_p = sEnvVars.Begin;
    }
    else return 0;

    /* Initialize space. */
    memset(sEnvVars.Begin, 0, sEnvVarArrayStep * sizeof(hd_env_var));

  } else {
    /* See if the requested environment variable is already defined,
       and if so, simply return it. */
    for (v_p = sEnvVars.Begin; v_p != sEnvVars.End; ++v_p) {
      if (!v_p->Name) break; /* Stop looking when first empty slot is found. */
      if (0 == strcmp(v_p->Name, name)) return v_p;
    }
  }

  /* Environment variable is not yet defined. See if there is
     a free spot in the environment variable array. */
  while (v_p->Name && v_p != sEnvVars.End) ++v_p;

  /* If the end of the array was reached without finding a free spot,
     increase the size of the array table. */
  if (v_p == sEnvVars.End) {
    size_t offset = v_p - sEnvVars.Begin; /* Current number of elements. */

    /* Allocate more space. */
    sEnvVars.Begin = (hd_env_var *) realloc(
        sEnvVars.Begin, (offset + sEnvVarArrayStep) * sizeof(hd_env_var)
    );

    if (sEnvVars.Begin) {
      /* Space was successfully realloc'ed. Update pointers. */
      v_p = sEnvVars.Begin + offset;
      sEnvVars.End = v_p + sEnvVarArrayStep;

      /* Initialize the newly-allocated portion of the space. */
      memset(v_p, 0, sEnvVarArrayStep * sizeof(hd_env_var));
    } else return 0;
  }

  /* Finally, copy the name to the new environment variable. */
  v_p->Name = CpyStr(name);

  return v_p;
}

const char *HD_getenv(const char *name) {
  hd_env_var *v_p;
  if (!name) return 0;
  v_p = FindEnvVar(name);
  if (!v_p->Value) v_p->Value = CpyStr(getenv(name));
  return v_p->Value;
}

const char *HD_setenv(const char *name, const char *value) {
  hd_env_var *v_p;
  if (!name) return 0;
  v_p = FindEnvVar(name);
  if (v_p->Value) free(v_p->Value);
  v_p->Value = CpyStr(value);
  return v_p->Value;
}

/* ------- f77 wrappers -------- */

/* This is a f77-callable routine to write
   to the dedicated stdout stream (heaout)
   controlled by HEADASOUTPUT env var
   Usage: call hdecho('my text') */
void headas_f77echo(const char *txt){

    fprintf(heaout,"%s\n",txt);

    return;
}
FCALLSCSUB1(headas_f77echo, HDECHO, hdecho, STRING)

/* This is a f77-callable routine to write to
   stderr (subject to the HEADASERROR env var)
   Usage: call hderr('my error message') */
void headas_f77err(const char *txt){

    fprintf(heaerr,"%s\n",txt);

    return;
}
FCALLSCSUB1(headas_f77err, HDERR, hderr, STRING)


/* This is a f77-callable wrapper for headas_chat()
   Usage: call hdchat(3,'my text') */
#define HDSTRLIM 80
void headas_f77chat(int threshold, const char *txt){
    char msg[HDSTRLIM+2];

    sprintf(msg, "%s\n", txt);
    headas_chat(threshold, msg);

    return;
}
FCALLSCSUB2(headas_f77chat, HDCHAT, hdchat, INT, STRING)

#define GetFilePtr(A) (A && A->cStream ? *A->cStream : 0)

  static hd_FILE *InitStream(FILE * restrict fptr) {
    hd_FILE *retval;
    do {
      if(!fptr) {
        retval = 0;
        continue;
      }

      retval = (hd_FILE *) malloc(sizeof(hd_FILE));
      if(!retval) continue;

      retval->cStream = (FILE **) malloc(sizeof(FILE *));

      if(retval->cStream) {
        *retval->cStream = fptr;
        retval->cxxStream = 0;
      } else {
        free(retval);
        retval = 0;
        errno = ENOMEM;
      }
    } while(0);
    return retval;
  }

  static FILE *CloseStream(hd_FILE * restrict theStream) {
    FILE *retval = 0;
    do {
      if(!theStream) continue;

      if(theStream && theStream->cxxStream);
      else {
        FILE **tmp = theStream->cStream;
        if(tmp) retval = *tmp;
        theStream->cStream = 0;
        if(theStream != hd_err && theStream != hd_in && theStream != hd_out) {
          free(tmp);
          free(theStream);
        }
      }
    } while(0);
    return retval;
  }

  /* The order of, and prototypes for stdio replacements are taken
     from the ISO C standard, Chapter 7. */
  /* 7.19.4 Operations on files. */
  int HD_remove(const char *fileName) {
    return remove(fileName);
  }

  int HD_rename(const char *oldName, const char *newName) {
    return rename(oldName, newName);
  }

  hd_FILE *HD_tmpfile(void) {
    return InitStream(tmpfile());
  }

  char *HD_tmpnam(char *fileName) {
   return tmpnam(fileName);
  }

  /* 7.19.5 File access functions. */
  int HD_fclose(hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fclose(CloseStream(theStream));
    } while(0);
    return retval;
  }

  int HD_fflush(hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fflush(GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  hd_FILE *HD_fopen(const char * restrict fileName,
      const char * restrict mode) {
    return InitStream(fopen(fileName, mode));
  }

  hd_FILE *HD_freopen(const char * restrict fileName,
      const char * restrict mode, hd_FILE * restrict theStream) {
    return InitStream(freopen(fileName, mode, CloseStream(theStream)));
  }

  void HD_setbuf(hd_FILE * restrict theStream, char * restrict buf) {
    do {
      if(theStream && theStream->cxxStream);
      else setbuf(GetFilePtr(theStream), buf);
    } while(0);
  }

  void HD_setvbuf(hd_FILE * restrict theStream, char * restrict buf,
      int mode, size_t size) {
    do {
      if(theStream && theStream->cxxStream);
      else setvbuf(GetFilePtr(theStream), buf, mode, size);
    } while(0);
  }

  /* 7.19.6 Formatted input/output functions. */
  int HD_fprintf(hd_FILE * restrict theStream, const char * restrict fmt, ...) {
    int retval = -1;
    va_list ap;
    va_start(ap, fmt);
    do {
      if(theStream && theStream->cxxStream);
      else retval = vfprintf(GetFilePtr(theStream), fmt, ap);
    } while(0);
    va_end(ap);
    return retval;
  }

  int HD_fscanf(hd_FILE * restrict theStream, const char * restrict fmt, ...) {
    int retval = EOF;
    va_list ap;
    va_start(ap, fmt);
    do {
      if(theStream && theStream->cxxStream);
      else {
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
        retval = vfscanf(GetFilePtr(theStream), fmt, ap);
#else
        fprintf(GetFilePtr(hd_err), "HD_fscanf not yet implemented.\n");
        assert(0);
#endif
      }
    } while(0);
    va_end(ap);
    return retval;
  }

  int HD_printf(const char * restrict fmt, ...) {
    int retval = -1;
    va_list ap;
    va_start(ap, fmt);
    do {
      if(hd_out && hd_out->cxxStream);
      else retval = vfprintf(GetFilePtr(hd_out), fmt, ap);
    } while(0);
    va_end(ap);
    return retval;
  }

  int HD_scanf(const char * restrict fmt, ...) {
    int retval = EOF;
    va_list ap;
    va_start(ap, fmt);
    do {
      if(hd_in && hd_in->cxxStream);
      else {
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
        retval = vfscanf(GetFilePtr(hd_in), fmt, ap);
#else
        fprintf(GetFilePtr(hd_err), "HD_scanf not yet implemented.\n");
        assert(0);
#endif
      }
    } while(0);
    va_end(ap);
    return retval;
  }

  int HD_snprintf(char * restrict s, size_t n, const char * restrict fmt, ...) {
    int retval = -1;
    va_list ap;
    va_start(ap, fmt);
    fprintf(GetFilePtr(hd_err), "HD_snprintf not implemented.\n");
    assert(0);
    va_end(ap);
    return retval;
  }

  int HD_sprintf(char * restrict s, const char * restrict fmt, ...) {
    int retval = -1;
    va_list ap;
    va_start(ap, fmt);
    retval = vsprintf(s, fmt, ap);
    va_end(ap);
    return retval;
  }

  int HD_sscanf(const char * restrict s, const char * restrict fmt, ...) {
    int retval = EOF;
    va_list ap;
    va_start(ap, fmt);
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
    retval = vsscanf(s, fmt, ap);
#else
    fprintf(GetFilePtr(hd_err), "HD_sscanf not yet implemented.\n");
    assert(0);
#endif
    va_end(ap);
    return retval;
  }

  int HD_vfprintf(hd_FILE * restrict theStream, const char * restrict fmt,
      va_list ap) {
    int retval = -1;
    do {
      if(theStream && theStream->cxxStream);
      else retval = vfprintf(GetFilePtr(theStream), fmt, ap);
    } while(0);
    return retval;
  }

  int HD_vfscanf(hd_FILE * restrict theStream, const char * restrict fmt,
      va_list ap) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else {
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
        retval = vfscanf(GetFilePtr(theStream), fmt, ap);
#else
        fprintf(GetFilePtr(hd_err), "HD_vfscanf not yet implemented.\n");
        assert(0);
#endif
      }
    } while(0);
    return retval;
  }

  int HD_vprintf(const char * restrict fmt, va_list ap) {
    int retval = -1;
    do {
      if(hd_out && hd_out->cxxStream);
      else retval = vfprintf(GetFilePtr(hd_out), fmt, ap);
    } while(0);
    return retval;
  }

  int HD_vscanf(const char * restrict fmt, va_list ap) {
    int retval = EOF;
    do {
      if(hd_in && hd_in->cxxStream);
      else {
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
        retval = vfscanf(GetFilePtr(hd_in), fmt, ap);
#else
        fprintf(GetFilePtr(hd_err), "HD_vscanf not yet implemented.\n");
        assert(0);
#endif
      }
    } while(0);
    return retval;
  }

  int HD_vsnprintf(char * restrict s, size_t n, const char * restrict fmt,
      va_list ap) {
    int retval = -1;
    fprintf(GetFilePtr(hd_err), "HD_vsnprintf not implemented.\n");
    assert(0);
    return retval;
  }

  int HD_vsprintf(char * restrict s, const char * restrict fmt, va_list ap) {
    return vsprintf(s, fmt, ap);
  }

  int HD_vsscanf(const char * restrict s, const char * restrict fmt,
      va_list ap) {
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
    return vsscanf(s, fmt, ap);
#else
    fprintf(GetFilePtr(hd_err), "HD_vsscanf not yet implemented.\n");
    assert(0);
    return EOF;
#endif
  }

  /* 7.19.7 Character input/output functions. */
  int HD_fgetc(hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fgetc(GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  char *HD_fgets(char * restrict s, int n, hd_FILE * restrict theStream) {
    char *retval = 0;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fgets(s, n, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  int HD_fputc(int c, hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fputc(c, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  int HD_fputs(const char * restrict s, hd_FILE * restrict theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fputs(s, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  int HD_getc(hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = getc(GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  int HD_getchar(void) {
    int retval = EOF;
    do {
      if(hd_in && hd_in->cxxStream);
      else retval = getc(GetFilePtr(hd_in));
    } while(0);
    return retval;
  }

  char *HD_gets(char *s) {
    char *retval = s;
    fprintf(heaerr, "Do not use HD_gets. Use HD_fgets instead.\n");
    assert(0);
    return retval;
  }

  int HD_putc(int c, hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = putc(c, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  int HD_putchar(int c) {
    int retval = EOF;
    do {
      if(hd_out && hd_out->cxxStream);
      else retval = putc(c, GetFilePtr(hd_out));
    } while(0);
    return retval;
  }

  int HD_puts(const char *s) {
    int retval = EOF;
    do {
      if(hd_out && hd_out->cxxStream);
      else retval = fputs(s, GetFilePtr(hd_out));
    } while(0);
    return retval;
  }

  int HD_ungetc(int c, hd_FILE *theStream) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = ungetc(c, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  /* 7.19.8 Direct input/output functions. */
  size_t HD_fread(void * restrict ptr, size_t size, size_t nmemb,
      hd_FILE * restrict theStream) {
    size_t retval = 0u;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fread(ptr, size, nmemb, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  size_t HD_fwrite(const void * restrict ptr, size_t size, size_t nmemb,
      hd_FILE * restrict theStream) {
    size_t retval = 0u;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fwrite(ptr, size, nmemb, GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  /* 7.19.9 File positioning functions. */
  int HD_fgetpos(hd_FILE * restrict theStream, fpos_t * restrict pos) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fgetpos(GetFilePtr(theStream), pos);
    } while(0);
    return retval;
  }

  int HD_fseek(hd_FILE *theStream, long int offset, int whence) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fseek(GetFilePtr(theStream), offset, whence);
    } while(0);
    return retval;
  }

  int HD_fsetpos(hd_FILE *theStream, const fpos_t *pos) {
    int retval = EOF;
    do {
      if(theStream && theStream->cxxStream);
      else retval = fsetpos(GetFilePtr(theStream), pos);
    } while(0);
    return retval;
  }

  long int HD_ftell(hd_FILE *theStream) {
    long int retval = -1L;
    do {
      if(theStream && theStream->cxxStream);
      else retval = ftell(GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  void HD_rewind(hd_FILE *theStream) {
    do {
      if(theStream && theStream->cxxStream);
      else rewind(GetFilePtr(theStream));
    } while(0);
  }

  /* 7.19.10 Error-handling functions. */
  void HD_clearerr(hd_FILE *theStream) {
    do {
      if(theStream && theStream->cxxStream);
      else clearerr(GetFilePtr(theStream));
    } while(0);
  }

  int HD_feof(hd_FILE *theStream) {
    int retval = -1;
    do {
      if(theStream && theStream->cxxStream);
      else retval = feof(GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  int HD_ferror(hd_FILE *theStream) {
    int retval = -1;
    do {
      if(theStream && theStream->cxxStream);
      else retval = ferror(GetFilePtr(theStream));
    } while(0);
    return retval;
  }

  void HD_perror(const char *s) {
    perror(s);
  }
