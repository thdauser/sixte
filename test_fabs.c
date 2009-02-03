#include "test_fabs.h"

int main(void) {
  printf("|-11.2| = %g\n", betrag(-11.2));

  return EXIT_SUCCESS;
}

double betrag(double wert) {
  printf("cos(PI) = %g\n", cos(M_PI));

  return(fabs(wert));
}
