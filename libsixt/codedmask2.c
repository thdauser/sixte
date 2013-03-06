int getImpactPos (struct Point2d* const position,
		  const vector* const phodir,
		  const CodedMask* const mask, 
		  const struct Telescope* const telescope,
		  const float distance,
		  int* const status)
{
  // Check if a CodedMask is specified. If not, break.
  if (NULL==mask) return(0);

  // Get a random number.
  double rand=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0); 

  //Check if photon is absorbed by an opaque pixel.
  //Use random number an compare with transparency.
   if (rand>mask->transparency) {
    // The photon is absorbed.
    return(0);
   }

   //Photon passes through transparent pixel.
   //Determine its impact position on the detection plane.

   //First:Determine the pixel(random)the photon passes through (mask plane).
   int pixel = (int)(rand/mask->transparency * mask->n_transparent_pixels);
   //Pixel points to an arbitrary pixel out of all transparent ones.

   //Second:Determine the components of the photon direction with respect to the
   //detector coordinate axes nx, ny (in mask plane).
   double x_comp = scalar_product(&phodir, &telescope->nx);
   double y_comp = scalar_product(&phodir, &telescope->ny);

   //Third:Determine the component of the photon direction within the mask plane.
   double radius = sqrt(pow(x_comp,2.)+pow(y_comp,2.));
   //And the azimuthal angle (with respect to the nx-axis)
   double alpha=atan2(y_comp, x_comp);

   //Now get the impact positon in the mask plane in meters.
   //(pixel_pos - reference_pixel)*width of one pixel+value at refernce_pixel
   position->x =
     ((double)(mask->transparent_pixels[pixel][0])-mask->crpix1)
     *mask->cdelt1+mask->crval1;
   CHECK_STATUS_RET(*status, 0);
   position->y =
     ((double)(mask->transparent_pixels[pixel][1])-mask->crpix2)
     *mask->cdelt2+mask->crval2;
   CHECK_STATUS_RET(*status, 0);

   //Finally get the impact position in the detection plane.
   //Shift the above according to the off-axis position and
   //the distance mask-detector.
   position->x -= cos(alpha) * radius * distance;
   position->y -= sin(alpha) * radius * distance;

  return(1);
}

