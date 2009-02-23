#ifndef STRFTCPY_H
#define STRFTCPY_H 1


#include <stdlib.h>
#include <string.h>


// copies a substring of 'start' to 'destination'
void strftcpy(char *dest, char *source, int start, int length);

#endif
