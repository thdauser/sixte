#include "erodetbkgrndgen.h"

int main(int argc, char** argv) {
  int status=EXIT_SUCCESS;

  if(argc < 2) {
    printf("too few arguments!\n");
    exit(1);
  }

  eroBkgInitialize("merged_hitlist.fits", &status);
  eroBkgSetRateFct("bkglc.simput", &status);
  int cc;
  for(cc=1; cc<argc; cc++) {
    eroBackgroundOutput* output = (eroBackgroundOutput*) malloc(sizeof(eroBackgroundOutput));
    output = eroBkgGetBackgroundList(atof(argv[cc]));
    // int ii;
    //    for(ii = 0; ii < output->numhits; ii++) {
    //        printf("[%f|%f] @ %f : %f\n", output->hit_xpos[ii], output->hit_ypos[ii], output->hit_time[ii], output->hit_energy[ii]);
    //    }
    printf("number of events: %d\n", output->numevents);
    printf("[%d] --------------------------------------------\n", cc);
  }
  eroBkgCleanUp(&status);

  return(status);
} 
