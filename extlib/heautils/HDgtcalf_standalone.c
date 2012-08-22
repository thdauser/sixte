/******************************************************************************
 *   File name: HDgtcalf_standalone.c                                         *
 *                                                                            *
 * Description: Definitions necessary for HDgtcalf to be used outside         *
 *     the context of Headas.                                                 *
 *                                                                            *
 *      Author: James Peachey, L3 Communications, for HEASARC/GSFC/NASA       *
 *                                                                            *
 *****************************************************************************/
#include "HDgtcalf_standalone.h"

/* No-op replacement for HDerror_throw. */
int HDgtcalf_error_throw(const char * msg, int status) { return status; }
