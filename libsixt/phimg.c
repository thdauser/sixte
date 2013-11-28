#include "phimg.h"


int phimg(const GenTel* const tel,
	  Attitude* const ac,
	  Photon* const ph,
	  Impact* const imp,
	  int* const status)
{
  // Calculate the minimum cos-value for sources inside the FOV: 
  // (angle(x0,source) <= 1/2 * diameter)
  const double fov_min_align=cos(tel->fov_diameter/2.); 

  // Determine the telescope pointing direction at the current time.
  struct Telescope telescope;
  telescope.nz=getTelescopeNz(ac, ph->time, status);
  CHECK_STATUS_RET(*status, 0);

  // Check whether the photon is inside the FOV.
  // Compare the photon direction to the direction of the telescope axis.
  Vector photon_direction=unit_vector(ph->ra, ph->dec);
  if (check_fov(&photon_direction, &telescope.nz, fov_min_align)==0) {
    // Photon is inside the FOV!
	
    // Determine telescope data like pointing direction (attitude) etc.
    // The telescope coordinate system consists of an x-, y-, and z-axis.
    getTelescopeAxes(ac, &telescope.nx, &telescope.ny, &telescope.nz, 
		     ph->time, status);
    CHECK_STATUS_RET(*status, 0);

    // Determine the photon impact position on the detector (in [m]).

    // Convolution with PSF:
    // Function returns 0, if the photon does not fall on the detector. 
    // If it hits the detector, the return value is 1.
    struct Point2d position;
    int get_psf_pos_retval=
      get_psf_pos(&position, *ph, telescope, tel->focal_length, 
		  tel->vignetting, tel->psf, status);
    CHECK_STATUS_RET(*status, 0);
    if (get_psf_pos_retval!=0) {
      // Check whether the photon hits the detector within the FOV. 
      // (Due to the effects of the mirrors it might have been scattered over 
      // the edge of the FOV, although the source is inside the FOV.)
      if (sqrt(pow(position.x,2.)+pow(position.y,2.)) <
	  tan(tel->fov_diameter)*tel->focal_length) {
	
	// New impact.
	imp->time    =ph->time;
	imp->energy  =ph->energy;
	imp->position=position;
	imp->ph_id   =ph->ph_id;
	imp->src_id  =ph->src_id;

	return(1);
      } else {
	return(0);
      }
    } else {
      return(0);
    } 
    // END get_psf_pos(...)
  } else {
    return(0);
  }
  // End of FOV check.
}

