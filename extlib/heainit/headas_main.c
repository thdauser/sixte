/******************************************************************************
 *   File name: headas_main.c                                                 *
 *                                                                            *
 * Description: Universal implementation of main, which may be included       *
 *     by any HEADAS binary. This code performs the necessary start-up        *
 *     operations, then calls the task's top-level function and finally       *
 *     performs the necessary shut-down operations.                           *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: Mike Tripicco, for HEASARC/GSFC/NASA                          *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 *                                                                            *
 * Usage:                                                                     *
 *                                                                            *
 *     #define TOOLSUB name_of_top_level_tool_function                        *
 *     #include "headas_main.c"                                               *
 ******************************************************************************/

#ifndef TOOLSUB
#error: TOOLSUB is not defined
#endif

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include "headas.h"
/******************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Function declarations.                                                   *
   ****************************************************************************/
  int headas_init(int argc, char **argv);
  int headas_close(int errNum);
  /****************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

int TOOLSUB(void);

/* Universal main function. */
int main(int argc, char *argv[]) {
  int status;

  status = headas_init(argc, argv);

  if(0 == status) status = TOOLSUB();

  status = headas_close(status);

  return status;
}

/******************************************************************************
 * $Log: headas_main.c,v $
 * Revision 1.3  2002/10/04 21:51:03  peachey
 * Changes to use the new error handling facility automatically.
 *
 ******************************************************************************/
