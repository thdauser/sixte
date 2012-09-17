/******************************************************************************
 *   File name: HDgtcalf_standalone.h                                         *
 *                                                                            *
 * Description: Declarations necessary for HDgtcalf to be used outside        *
 *     the context of Headas.                                                 *
 *                                                                            *
 *      Author: James Peachey, L3 Communications, for HEASARC/GSFC/NASA       *
 *                                                                            *
 *****************************************************************************/
#ifndef HDGTCALF_STANDALONE_H
#define HDGTCALF_STANDALONE_H

/* The following definition of HD_ERROR_THROW should only be
   used by HDgtcalf in standalone mode. */
#ifndef HD_ERROR_THROW
#define HD_ERROR_THROW(A, B) HDgtcalf_error_throw(A, B)
#endif

/* The following enum was copied from headas_error. Only the name
   was changed. */
/* Enumerated error codes for HEAdas standard errors. */
typedef enum hdgtcalf_hd_error_code {
  HD_OK = 0,
  HD_ERR_BASE = 1000,
  HD_ERR_DYN_ALLOC_FAIL,
  HD_ERR_PROGRAMMER_ERRORS = HD_ERR_BASE + 500,
  HD_ERR_MNGR_ERROR,
  HD_ERR_NULL_POINTER,
  HD_ERR_OUT_OF_BOUNDS
} hdgtcalf_hd_error_code;

/* No-op replacement for HDerror_throw. */
int HDgtcalf_error_throw(const char * msg, int status);

#endif
