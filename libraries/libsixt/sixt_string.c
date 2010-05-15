#include "sixt_string.h"


void strtoupper(char* string) {
  int count=0;
  while (string[count] != '\0') {
    string[count] = toupper(string[count]);
    count++;
  };
}

