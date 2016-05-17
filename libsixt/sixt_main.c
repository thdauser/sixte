//
// sixt_main.c: modification of headas_main.c that preserves
// argc/argv as sixt_argc and sixt_argv for those SIXT tools that
// need argc/argv access
//

// Otherwise the tool works identically to headas_main

#ifndef TOOLSUB
#error: TOOLSUB is not defined
#endif

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include "headas.h"
/******************************************************************************/

int sixt_argc=-1;
char **sixt_argv;

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
int main(int argc, char **argv) {
  int status;

  sixt_argc=argc;
  sixt_argv=argv;
  
  status = headas_init(argc, argv);

  if(0 == status) status = TOOLSUB();

  status = headas_close(status);

  return status;
}
