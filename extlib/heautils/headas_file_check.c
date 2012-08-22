/******************************************************************************
 *   File name: headas_file_check.c                                           *
 *                                                                            *
 * Description: Function to check for existence, and/or access mode, of a     *
 *              file. Uses cfitsio's fits_file_exists function to determine   *
 *              this, allowing this function to be used for FITS files, which *
 *              may not literally exist with the given name, yet still may be *
 *              accessible.                                                   *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: James Peachey, for HEASARC/GSFC/NASA                          *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include <string.h>

#include "fitsio.h"
#include "longnam.h"
#include "headas_utils.h"
#include "ape/ape_error.h"
#include "ape/ape_par.h"
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
   * Type definitions.                                                        *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Global variable definitions.                                             *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Static variable definitions.                                             *
   ****************************************************************************/
  static int (*s_ape_default_file_check)(const char *, const char *) = 0;
  /****************************************************************************/

  /****************************************************************************
   * Static function definitions.                                             *
   ****************************************************************************/
  /* Warning: despite the similarity in names, this function has the reverse
     convention from HDfile_system_check, in that it returns 1 to indicate
     the file is accessible in the given mode, and 0 to indicate problems. */
  /* TODO: Remove this function if it is not used explicitly in HEADAS. Instead
     call s_ape_default_file_check directly. */
  int HDfile_system_check(const char *file_name, const char *open_mode) {
    int status = eOK;
    int access_ok = 0;

    /* Perform initialization once. */
    if (0 == s_ape_default_file_check) status = HDfile_check_init();

    if (eOK == status) {
      /* Use Ape's default file checker to do the file-system part of the
         check. */
      status = (*s_ape_default_file_check)(file_name, open_mode);
    }

    access_ok = (eOK == status) ? 1 : 0;

    return access_ok;
  }
  /****************************************************************************/

  /****************************************************************************
   * Function definitions.                                                    *
   ****************************************************************************/
  /* TODO: Modify heainit to call this instead of using PIL call. */
  int HDfile_check_init(void) {
    int status = eOK;
    if (0 == s_ape_default_file_check) {
      /* Restore ape to its default function. If it was already using the
         default this is a no-op. */
      status = ape_util_set_file_check_func(0);

      if (eOK == status) {
        /* Get pointer to the default ape function. */
        status = ape_util_get_file_check_func(&s_ape_default_file_check);
      }

      if (eOK == status && 0 == s_ape_default_file_check) {
        /* This should not happen. Ape should always have a default file
           checking function. */
        status = eNullPointer;
      }
    }

    if (eOK == status) {
      /* Change Ape to call HDfile_check for access checks. */
      status = ape_util_set_file_check_func(&HDfile_check);
    }

    return status;
  }

  /* Warning: despite the similarity in names, this function has the reverse
     convention from HDfile_system_check, in that it returns 0 to indicate
     the file is accessible in the given mode, and 1 to indicate problems. */
  int HDfile_check(const char *file_name, const char *open_mode) {
    /* First, try asking the file system directly. */
    int exists = HDfile_system_check(file_name, open_mode);

    if (0 == exists) {
      int status = 0;
      int cfitsio_exists = 0;
      /* File system failed to find the file. Try cfitsio. */
      fits_file_exists(file_name, &cfitsio_exists, &status);
      if (0 == status) {
        switch (cfitsio_exists) {
          case -1:
            /* Cfitsio believes it's a web-only file: accessible if mode is not write. */
            if (0 == open_mode || 0 == strstr(open_mode, "w")) exists = 1;
            break;
          case 0:
            /* Cfitsio believes file does not exist at all: accessible if mode is write. */
            if (0 != open_mode && 0 != strstr(open_mode, "w")) exists = 1;
            break;
          case 1:
            /* Cfitsio believes file exists on disk: must assume accessible, but who knows for sure? */
            exists = 1;
            break;
          case 2:
            /* Cfitsio believes file exists on disk as compressed file: accessible if mode is not write */
            if (0 == open_mode || 0 == strstr(open_mode, "w")) exists = 1;
            break;
          default:
            break;
        }
      }
    }

    /* Reverse the sense of this test, return Ape status OK if file check
       succeeded, not-accessible status if error. */
    return exists == 0 ? eFileNotAccessible : eOK;
  }
  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log: headas_file_check.c,v $
 * Revision 1.5  2007/12/04 15:54:55  peachey
 * Return the correct Ape status code signifying file not accessible
 * in the event of error, instead of just 0/1.
 *
 * Revision 1.4  2007/10/10 21:26:35  peachey
 * Use Ape's default file access checker for first stage of
 * file checking before performing additional checks using cfitsio.
 *
 * Revision 1.3  2006/06/23 20:55:13  peachey
 * Rework file checking mechanism to be more careful about interpreting
 * output from system access checks and cfitsio's fits_file_exists function.
 *
 * Revision 1.2  2004/09/02 20:13:30  elwin
 * Changed variable name to file_name in #ifdef WIN32 section of code. -- LEB
 *
 * Revision 1.1  2004/08/17 19:24:21  peachey
 * Add functions for testing for existence of files: HDfile_system_check,
 * which queries the file system, and HDfile_check, which calls
 * HDfile_system_check, but then also uses cfitsio's fits_file_exists function
 * to check for gzipped files, extended syntax, etc.
 *
 ******************************************************************************/
