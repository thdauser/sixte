/******************************************************************************
 *   File name: headas_utils.h                                                *
 *                                                                            *
 * Description: Public interface for HEADAS utilities library.                *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: James Peachey, for HEASARC/GSFC/NASA                          *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/
#ifndef HEADAS_UTILS_H
#define HEADAS_UTILS_H

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include "fitsio.h"
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
  /****************************************************************************/

#ifndef IMPSYM
#ifdef WIN32
#define IMPSYM __declspec(dllimport)
#else
#define IMPSYM
#endif
#endif

  /****************************************************************************
   * Global variable forward declarations.                                    *
   ****************************************************************************/
  extern IMPSYM int headas_clobpar;
  /****************************************************************************/

  /****************************************************************************
   * Function declarations.                                                   *
   ****************************************************************************/
  void set_history(int hpar);
  void get_history(int *hvar);

  void set_toolname(const char *name);
  void get_toolname(char *name);
  void set_toolversion(const char *vers);
  void get_toolversion(char *vers);
  void get_toolnamev(char *str);

  const char *hdbasename(const char *path);
  int headas_parstamp(fitsfile *fptr, int hdunum);
  int HDpar_stamp(fitsfile *fptr, int hdunum, int *status);
  void HDpar_note(const char * name, const char * format, ...);
  int headas_clobberfile(char *filename);
  float hd_ran2(long *idum);
  char **expand_item_list(char *liststr, int *nitems, char fieldsep,
			  int trim, int skipempty, int guardparen,
			  int *status);
  int HDfile_system_check(const char *file_name, const char *open_mode);
  int HDfile_check(const char *file_name, const char *open_mode);

  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

#endif

/******************************************************************************
 * $Log: headas_utils.h,v $
 * Revision 1.8  2005/08/09 19:05:02  elwin
 * Added import specification for Win32 dll build. -- LEB
 *
 * Revision 1.7  2005/05/06 15:05:05  miket
 * adding new function HDpar_note() contributed by R.Wiegand
 *
 * Revision 1.6  2004/08/17 19:24:21  peachey
 * Add functions for testing for existence of files: HDfile_system_check,
 * which queries the file system, and HDfile_check, which calls
 * HDfile_system_check, but then also uses cfitsio's fits_file_exists function
 * to check for gzipped files, extended syntax, etc.
 *
 * Revision 1.5  2003/08/29 14:14:27  irby
 * Needed a semi-colon on line 56...
 *
 * Revision 1.4  2003/08/29 13:47:19  miket
 * added HDpar_stamp
 *
 * Revision 1.3  2003/04/04 21:21:14  miket
 * merging changes from swift_bld4_0 branch
 *
 * Revision 1.2.2.1  2003/03/26 21:01:55  miket
 * added Craig's expand_item_list routine and use it in headas_parstamp
 *
 * Revision 1.2  2002/10/04 21:49:43  peachey
 * Simplify hdbasename function.
 *
 * Revision 1.1  2002/10/04 18:56:47  peachey
 * New self-contained public header file for heautils library.
 *
 ******************************************************************************/
