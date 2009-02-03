#include "strftcpy.h"

// The following function copies a substring of specified length out of source, 
// starting with character start
void strftcpy(char *dest, char *source, int start, int length) {
  char temp[length];
  int i;

  // copy the individual characters
  for(i=0;i<length;i++) {
    temp[i] = source[start+i];
  }
  temp[length] = '\0';  // terminate string

  // copy temporary string to destination
  strcpy(dest, temp);

}



