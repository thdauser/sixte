#include <stdio.h>
#include "hdcal.h"
#include "HDgtcalf_internal.h"

#ifdef HDGTCALF_STANDALONE
#include "HDgtcalf_standalone.h"
#else
#include "headas_error.h"
#endif

/*
int HDgtcalf (const char* tele, const char* instr, const char* detnam, 
      const char* filt, const char* codenam, const char* strtdate,
      const char* strtime, const char* stpdate, const char* stptime,
      const char* expr, int maxret,int fnamesize, char** filenam,
      long* extno, char** online, int* nret, int* nfound, int* status);
*/
int main() {
  int status = HD_OK;
  int func_status = 137;

  func_status = HDgtcalf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &func_status);
  if (137 != func_status) {
    fprintf(stderr, "HDgtcalf called with a non-0 status returned %d, not 137\n", func_status);
    status = 1;
  }

  func_status = 0;

  HDgtcalf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &func_status);
  if (HD_ERR_NULL_POINTER != func_status) {
    fprintf(stderr, "HDgtcalf called with some 0s returned %d, not HD_ERR_NULL_POINTER (%d)\n", func_status,
      HD_ERR_NULL_POINTER);
    status = 1;
  }

  func_status = HDgtcalf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  if (HD_ERR_NULL_POINTER != func_status) {
    fprintf(stderr, "HDgtcalf called with all 0s returned %d, not HD_ERR_NULL_POINTER (%d)\n", func_status,
      HD_ERR_NULL_POINTER);
    status = 1;
  }

  return status;
}
