#ifndef GAUSSIANCHARGECLOUD_H
#define GAUSSIANCHARGECLOUD_H 1

typedef struct {

  /** Sigma value for Gaussian shape charge clouds (given in [m]).
      This quantity is used to calculate size/extension of the
      Gaussian shape charge cloud. */
  double ccsigma;
  /** Size of the charge cloud (given in [m]). Quantity to estimate
      the extension of a Gaussian shape charge cloud. The value is
      defined to be three times ccsigma. It is assumed that
      approximately all charge is with in a radius of this size (3
      sigma) around the center of the charge cloud. */
  double ccsize; 

} GaussianChargeCloud;

#endif /* GAUSSIANCHARGECLOUD_H */

