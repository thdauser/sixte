#include "split.h"


void split_events2(
		  double det_xa, double det_ya,     // coordinates of event [real detector coordinates]
		  double fraction[3][3]             // energy fractions falling on neighboring pixels
		  ) 
{
  // Set the split charge fractions to 0 at the beginning.
  int countx, county;
  for (countx=0; countx<3; countx++) {
    for (county=0; county<3; county++) {
      fraction[countx][county] = 0.;
    }
  }

  // TODO !!!!!!!!!!!!!!!
  // Definition of real detector coordinates ??
  // Is det_xa = 0 at the center of the pixel or at the border ??
  // !!!!!!!!!!!!!!!!

  // Calculate the index of the neighboring pixels that might be involved in the charge splitting.
  double xd, yd;
  xd = det_xa-(int)det_xa;
  yd = det_ya-(int)det_ya;
  int xneighbor, yneighbor;
  xneighbor = ceil(xd) + 1;
  yneighbor = ceil(yd) + 1;


  // Calculate the (so far not normalized !) charge fractions of the 4 neighboring pixels.
  double r, sum;
  // centering pixel
  r = sqrt(pow(xd,2.)+pow(yd,2.));
  fraction[1][1] = exp(-pow(r/split_index, 2.));
  sum = fraction[1][1];

  // neighbor in x-direction
  r = sqrt(pow(1.-fabs(xd),2.)+pow(yd,2.));
  fraction[xneighbor][1] = exp(-pow(r/split_index, 2.));
  sum += fraction[xneighbor][1];

  // neighbor in y-direction
  r = sqrt(pow(xd,2.)+pow(1.-fabs(yd),2.));
  fraction[1][yneighbor] = exp(-pow(r/split_index, 2.));
  sum += fraction[1][yneighbor];

  // neighbor in diagonal direction
  r = sqrt(pow(1.-fabs(xd),2.)+pow(1.-fabs(yd),2.));
  fraction[xneighbor][yneighbor] = exp(-pow(r/split_index, 2.));
  sum += fraction[xneighbor][yneighbor];


  // Rescale the charge fraction in the individual pixels, so that the total charge is conserved.
  fraction[1][1] = fraction[1][1]/sum;
  if (fraction[1][1] < split_threshold) { fraction[1][1] = 0; }

  fraction[xneighbor][1] = fraction[xneighbor][1]/sum;
  if (fraction[xneighbor][1] < split_threshold) { fraction[xneighbor][1] = 0; }

  fraction[1][yneighbor] = fraction[1][yneighbor]/sum;
  if (fraction[1][yneighbor] < split_threshold) { fraction[1][yneighbor] = 0; }

  fraction[xneighbor][yneighbor] = fraction[xneighbor][yneighbor]/sum;
  if (fraction[xneighbor][yneighbor] < split_threshold) { fraction[xneighbor][yneighbor] = 0; }

}
