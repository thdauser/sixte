#ifndef SPLIT_H
#define SPLIT_H 1

#include <math.h>


const float split_index = 0.35;       // The charge distribution of a split event is ~ exp(-(r/split_index)^2) .
const float split_threshold = 0.01;   // If the split charge in a pixel is below the threshold, it is not regarded as a split event.

#endif
