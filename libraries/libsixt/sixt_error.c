#include "sixt.h"

void sixt_error(const char* const func, const char* const msg)
{
  // Print the formatted output message.
  //printf("Error in %s: %s!\n", func, msg);

  // Use the HEADAS error output routine.
  char output[MAXMSG];
  sprintf(output, "Error in %s: %s!\n", func, msg);
  HD_ERROR_THROW(output, EXIT_FAILURE);
}
