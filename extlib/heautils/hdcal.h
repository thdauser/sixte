#ifndef HDCAL_H
#define HDCAL_H

#ifdef __cplusplus
extern "C" {
#endif

int HDgtcalf (const char* tele, const char* instr, const char* detnam, 
  const char* filt, const char* codenam, const char* strtdate,
  const char* strtime, const char* stpdate, const char* stptime,
  const char* expr, int maxret,int fnamesize, char** filenam,
  long* extno, char** online, int* nret, int* nfound, int* status);

#ifdef __cplusplus
}
#endif

#endif
