/******************************************************************************
 *   File name: headas_stdio.h                                                *
 *                                                                            *
 * Description: Public C/C++ callable HEADAS replacemenets for C standard     *
 *     I/O facilities.                                                        *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: James Peachey, for HEASARC/GSFC/NASA                          *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/
#ifndef HEADAS_STDIO_H
#define HEADAS_STDIO_H

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#ifdef __cplusplus
#include <iostream>
#endif

#include <stdarg.h>
#include <stdio.h>
/******************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Constants.                                                               *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Type declarations/definitions.                                           *
   ****************************************************************************/
  struct hd_FILE_s;
  typedef struct hd_FILE_s hd_FILE_t;
  typedef hd_FILE_t hd_FILE;

#ifndef ANSI_ISO_IEC_9899_1999_COMPLIANT
#define restrict
#endif

  /* On Windows, global variables which are exported must also be imported*/

#ifndef IMPSYM
#ifdef WIN32
#define IMPSYM __declspec(dllimport)
#else
#define IMPSYM
#endif
#endif

  /****************************************************************************/

  /****************************************************************************
   * Global variable forward declarations.                                    *
   ****************************************************************************/
  extern IMPSYM hd_FILE *hd_err;
  extern IMPSYM hd_FILE *hd_in;
  extern hd_FILE *hd_out;
  extern IMPSYM hd_FILE *hd_prompt;
  extern IMPSYM hd_FILE *hd_null;
  extern IMPSYM FILE *heaerr;
  extern IMPSYM FILE *heain;
  extern IMPSYM FILE *heaout;
  extern IMPSYM FILE *heaprom;
  extern IMPSYM FILE *heanullstream;
  extern IMPSYM int headas_chatpar;
  /****************************************************************************/

  /****************************************************************************
   * Function declarations.                                                   *
   ****************************************************************************/
  int HDIO_init(void);
  int headas_chat(int threshold, const char *fmt, ...);
  int headas_printf(const char *fmt, ...);
  int pil_printf(const char *fmt, ...);

  const char *HD_getenv(const char *name);
  const char *HD_setenv(const char *name, const char *value);

  /* The order of, and prototypes for stdio replacements are taken
     from the ISO C standard, Chapter 7. */
  /* 7.19.4 Operations on files. */
  int HD_remove(const char *fileName);
  int HD_rename(const char *oldName, const char *newName);
  hd_FILE *HD_tmpfile(void);
  char *HD_tmpnam(char *fileName);

  /* 7.19.5 File access functions. */
  int HD_fclose(hd_FILE *theStream);
  int HD_fflush(hd_FILE *theStream);
  hd_FILE *HD_fopen(const char * restrict fileName,
      const char * restrict mode);
  hd_FILE *HD_freopen(const char * restrict fileName,
      const char * restrict mode, hd_FILE * restrict theStream);
  void HD_setbuf(hd_FILE * restrict theStream, char * restrict buf);
  void HD_setvbuf(hd_FILE * restrict theStream, char * restrict buf,
      int mode, size_t size);

  /* 7.19.6 Formatted input/output functions. */
  int HD_fprintf(hd_FILE * restrict theStream, const char * restrict fmt, ...);
  int HD_fscanf(hd_FILE * restrict theStream, const char * restrict fmt, ...);
  int HD_printf(const char * restrict fmt, ...);
  int HD_scanf(const char * restrict fmt, ...);
  int HD_snprintf(char * restrict s, size_t n, const char * restrict fmt, ...);
  int HD_sprintf(char * restrict s, const char * restrict fmt, ...);
  int HD_sscanf(const char * restrict s, const char * restrict fmt, ...);
  int HD_vfprintf(hd_FILE * restrict theStream, const char * restrict fmt,
      va_list ap);
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
  int HD_vfscanf(hd_FILE * restrict theStream, const char * restrict fmt,
      va_list ap);
#endif
  
  int HD_vprintf(const char * restrict fmt, va_list ap);
  int HD_vscanf(const char * restrict fmt, va_list ap);
  int HD_vsnprintf(char * restrict s, size_t n, const char * restrict fmt,
      va_list ap);
  int HD_vsprintf(char * restrict s, const char * restrict fmt, va_list ap);
#ifdef ANSI_ISO_IEC_9899_1999_COMPLIANT
  int HD_vsscanf(const char * restrict s, const char * restrict fmt,
      va_list ap);
#endif

  /* 7.19.7 Character input/output functions. */
  int HD_fgetc(hd_FILE *theStream);
  char *HD_fgets(char * restrict s, int n, hd_FILE * restrict theStream);
  int HD_fputc(int c, hd_FILE *theStream);
  int HD_fputs(const char * restrict s, hd_FILE * restrict theStream);
  int HD_getc(hd_FILE *theStream);
  int HD_getchar(void);
  char *HD_gets(char *s);
  int HD_putc(int c, hd_FILE *theStream);
  int HD_putchar(int c);
  int HD_puts(const char *s);
  int HD_ungetc(int c, hd_FILE *theStream);

  /* 7.19.8 Direct input/output functions. */
  size_t HD_fread(void * restrict ptr, size_t size, size_t nmemb,
      hd_FILE * restrict theStream);
  size_t HD_fwrite(const void * restrict ptr, size_t size, size_t nmemb,
      hd_FILE * restrict theStream);

  /* 7.19.9 File positioning functions. */
  int HD_fgetpos(hd_FILE * restrict theStream, fpos_t * restrict pos);
  int HD_fseek(hd_FILE *theStream, long int offset, int whence);
  int HD_fsetpos(hd_FILE *theStream, const fpos_t *pos);
  long int HD_ftell(hd_FILE *theStream);
  void HD_rewind(hd_FILE *theStream);

  /* 7.19.10 Error-handling functions. */
  void HD_clearerr(hd_FILE *theStream);
  int HD_feof(hd_FILE *theStream);
  int HD_ferror(hd_FILE *theStream);
  void HD_perror(const char *s);
  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

#endif

/******************************************************************************
 * $Log: headas_stdio.h,v $
 * Revision 1.6  2005/08/09 18:55:14  elwin
 * Fixes to import global variables properly in the Win32 build. -- LEB
 *
 * Revision 1.5  2003/02/13 18:44:53  peachey
 * Added full support for prompt and null streams. Refined logic
 * concerning the prompt stream; if HEADASNOQUERY is set prompts
 * will be redirected to the null stream. Refined logic for other
 * streams to use the null stream when errors occur during
 * redirection.
 *
 * Revision 1.4  2003/02/12 23:31:06  peachey
 * Add HD_getenv and HD_setenv to handle getting and setting env. variables.
 *
 * Revision 1.3  2003/02/12 23:26:18  peachey
 * hdIOInit code moved here from headas_init.c, and renamed to HDIO_init.
 *
 * Revision 1.2  2002/11/01 19:19:10  peachey
 * Add full replacements for all stdio functions.
 *
 * Revision 1.1  2002/10/04 18:56:18  peachey
 * New self-contained header file for heaio library.
 *
 *****************************************************************************/
