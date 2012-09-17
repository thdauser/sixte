/******************************************************************************
 *   File name:                                                               *
 *                                                                            *
 * Description:                                                               *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: Bryan Irby, for HEASARC/GSFC/NASA                             *
 *              modified from code supplied by Bob Wiegand                    *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/
#ifndef HEADAS_LOCK_H
#define HEADAS_LOCK_H

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
/******************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Constants.                                                               *
   ****************************************************************************/

#define HD_TMPLEN 1024

  /****************************************************************************/

  /****************************************************************************
   * Type declarations/definitions.                                           *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Global variable forward declarations.                                    *
   ****************************************************************************/
  /****************************************************************************/

  /****************************************************************************
   * Function declarations.                                                   *
   ****************************************************************************/

  /* Caller should use or copy returned string immediately since it will be
   * altered on future calls: */
   const char * HD_tmpfile (const char * base, const char * ext);
  /* Caller should use or copy returned string immediately since it will be
   * altered on future calls: */
   const char * HDtmpdir ();
  /* length of buffer argument should be at least HD_TMPLEN */
   int HDtmpdir_r (char * buffer);
  /* length of buffer argument should be at least HD_TMPLEN */
   int HDtmpfile_r (char * buffer, const char * base, const char * ext);

  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

#endif

/******************************************************************************
 * $Log: headas_lock.h,v $
 * Revision 1.1  2005/03/23 17:24:31  irby
 * Add routines for temporary file and pseudo-random number generation.
 *
 ******************************************************************************/
