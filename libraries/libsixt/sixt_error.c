#include "sixt.h"

void sixt_error(const char* const func, const char* const msg)
{
  // Print the formatted output message.
  printf("Error in %s: %s!\n", func, msg);
}
