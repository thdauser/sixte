/////////////////////////////////////////////////////////////////////////////////////////
//
// general information:
// file:             orbit_calc.c
// related to:       orbit_calc.h
// project:          eROSITA - NRTA - simulation
// author:           Christian Schmid
// date:             01/2008
//
/////////////////////////////////////////////////////////////////////////////////////////
//
// description:
// This file belongs to a library out of the NRTA-simulation for the eROSITA mission.
// It contains several functions to calculate the orbit of a satellite including
// perturbations by the particular shape of the earth.
// The routines are implemented for maximum speed, so there are some redundant variables.
//
/////////////////////////////////////////////////////////////////////////////////////////


#include "orbit_calc.h"


// This compiler flag determines, whether the transformed set of orbital elements is used for
// the calculation of the orbit, whether the standard Keplerian elements are used.
#define ORBIT_TRANS_ELEM 1



//////////////////////////////////////////////////////////////////////////////////////////////////
// This function calculates the initial orbit data and saves it to the corresponding structure. //
//////////////////////////////////////////////////////////////////////////////////////////////////
void orbit_init(double a, double e, double i, double Omega, double omega, double M, struct orbit_data *odata) {
  // save basic orbit data in structure:
  odata->a = a;
  odata->p = odata->a*(1-powl(e,2.0));
  odata->e = e;
  odata->i = i;
  odata->Omega = Omega;
  odata->omega = omega;
  odata->M = M;

  // calculate mean motion:
  odata->n = sqrtl(mu/powl(odata->a,3.0));
}





//////////////////////////////////////////////////////////////////////////////////////////////
// This function performs a timestep, i.e. it calculates the next position of the satellite //
// on its orbit, considering perturbation terms up to first order (only J2).                //
//////////////////////////////////////////////////////////////////////////////////////////////
void orbit_step_J2(struct orbit_data *odata, struct orbit_position *oposition) {
  // variables
  long double E;                 // eccentric anomaly E(t)
  long double f;                 // true anomaly f(t)
  long double rl;                // length of vector r
  long double vr;                // radial velocity component
  long double vf;                // velocity component along orbit

  // calculation acceleration constants:
  const long double cosO = cosl(odata->Omega);
  const long double sinO = sinl(odata->Omega);
  const long double cosi = cosl(odata->i);
  const long double cosi2 = powl(cosi,2.0);
  const long double sini = sinl(odata->i);


  // calculate mean anomaly at time t (M(t)):
  odata->M += odata->n * odata->dt;

  // get the eccentric anomaly E(t) by solving the Kepler equation (Flury p. 18):
  E = kepler_equation(odata->M,odata->e);

  // calculate the true anomaly from the eccentric anomaly (Flury p. 16):
  f = 2.*atanl(sqrtl((1+odata->e)/(1-odata->e))*tanl(E/2));


  // perturbation effects (according to Nitschke p. 25)
  // calculate the factor s, that contributes to J2 orbit corrections, ds = dt * s:
  const double ds = odata->dt*(-1.5*J2*odata->n*powl(R_e/(odata->p),2.)); 

  // change of the ascending node:
  odata->Omega += ds * cosi;
  if (odata->Omega < 0.) {
    odata->Omega += 2.*M_PI;
  } else if (odata->Omega > 2.*M_PI) {
    odata->Omega -= 2.*M_PI;
  }
  // change of the argument of perigee:
  odata->omega += ds/2.*(1.-5.*cosi2);
  if (odata->omega < 0.) {
    odata->omega += 2.*M_PI;
  } else if (odata->omega > 2.*M_PI) {
    odata->omega -= 2.*M_PI;
  }
  // change of the mean eccentric anomaly:
  odata->M -= ds/2. * sqrtl(1.-powl(odata->e,2.)) * (3.*cosi2 - 1.);
  if (odata->M < 0.) {
    odata->M += 2.*M_PI;
  } else if (odata->M > 2.*M_PI) {
    odata->M -= 2.*M_PI;
  }
  // end of perturbation calculation


  // now we know the argument of latitude
  // and can calculate its cosine and sine
  const long double cosu = cosl(odata->omega+f);
  const long double sinu = sinl(odata->omega+f);
  
  // now we can determine the position and velocity of the satellite (Flury p. 38)
  rl = odata->p/(1.+odata->e*cosl(f));
  vr = sqrtl(mu/odata->p)*odata->e*sinl(f);
  vf = sqrtl(mu*odata->p)/rl;
  
  // position
  oposition->r.x = rl*(cosO*cosu-sinO*cosi*sinu);
  oposition->r.y = rl*(sinO*cosu+cosO*cosi*sinu);
  oposition->r.z = rl*sini*sinu;

  // velocity
  oposition->v.x = vr*(cosO*cosu-sinO*cosi*sinu) - vf*(cosO*sinu+sinO*cosu*cosi);
  oposition->v.y = vr*(sinO*cosu+cosO*cosi*sinu) - vf*(sinO*sinu-cosO*cosu*cosi);
  oposition->v.z = vr*sini*sinu + vf*cosu*sini;
  // end of unperturbed calculation

}





//////////////////////////////////////////////////////////////////////////////////////////////
// This function performs a timestep, i.e. it calculates the next position of the satellite //
// on its orbit, considering second order perturbations (J2, J2^2, J3, J4).                 //
//////////////////////////////////////////////////////////////////////////////////////////////
void orbit_step_J234(struct orbit_data *odata, struct orbit_position *oposition) 
{
  const long double J22 = powl(J2,2.);

// variables
  // orbital elements:
  long double a;        // semimajor axis
  long double e;        // eccentricity
  long double i;        // inclination
  long double Omega;    // right ascension of ascending node
  long double omega;    // argument of perigee
  long double M;        // mean anomaly

  // elements for low-eccentricity orbits
  long double h,k,u;
  
  long double E;                   // eccentric anomaly
  long double f;                   // true anomaly
  long double r;                   // lenght of position vector |r|



//////////////////////////////////////////////////////////
// change of mean orbital elements (secular perturbations)
  // required data:
  long double p = odata->a*(1.-powl(odata->e,2.));

  // frequently used values:
  long double Rep2 = powl(R_e/p,2.);
  long double Rep3 = powl(R_e/p,3.);
  long double Rep4 = powl(R_e/p,4.);
  long double sini = sinl(odata->i);
  long double sini2 = powl(sinl(odata->i),2.);
  long double sini4 = powl(sinl(odata->i),4.);
  long double cosi = cosl(odata->i);
  long double cosi2 = powl(cosl(odata->i),2.);
  long double cosi4 = powl(cosl(odata->i),4.);
  long double sino = sinl(odata->omega);
  long double coso = cosl(odata->omega);
  long double sin2o = sinl(2.*odata->omega);
  long double cos2o = cosl(2.*odata->omega);
  //  long double sino2 = powl(sinl(odata->omega),2.);
  long double coso2 = powl(cosl(odata->omega),2.);
  long double e2 = powl(odata->e,2.);
  long double ome2 = 1.-e2;
  long double some2 = sqrtl(ome2);


  // semimajor axis: 
  a = odata->a;  // a doesn't change (\dot{a}=0)


#ifndef ORBIT_TRANS_ELEM
  // eccentricity:
  e = odata->e + odata->dt*odata->n*( -3./32.*J22*Rep4*sini2*(14.-15.*sini2)*odata->e*ome2*sin2o 
				      -0.375*J3*Rep3*sini*(4.-5.*sini2)*ome2*coso
				      -15./32.*J4*Rep4*sini2*(6.-7.*sini2)*odata->e*ome2*sin2o );
#endif


  // inclination:
  i = odata->i + odata->dt*odata->n*( +3./64.*J22*Rep4*sinl(2.*odata->i)*(14.-15.*sini2)*e2*sin2o
				      +0.375*J3*Rep3*cosi*(4.-5.*sini2)*odata->e*coso
				      +15./64.*J4*Rep4*sinl(2.*odata->i)*(6.-7.*sini2)*e2*sin2o );


  // right ascension of ascending node:
  Omega = odata->Omega + odata->dt*odata->n*( -1.5*J2*Rep2*cosi 
					      -1.5*J22*Rep4*cosi*(2.25+1.5*some2-sini2*(2.5+2.25*some2)+e2/16.*(4.+5.*sini2) - e2*0.125*(7.-15.*sini2)*cos2o)
					      -0.375*J3*Rep3*(4.-5.*sini2)*odata->e/tanl(odata->i)*sino
					      +15./16.*J4*Rep4*cosi*((4.-7.*sini2)*(1.+1.5*e2)-(3.-7.*sini2)*e2*cos2o) );

				     
#ifndef ORBIT_TRANS_ELEM
  // argument of perigee:
  omega = odata->omega;
  // contributions from J2
  omega += odata->dt*odata->n*J2*Rep2*( 0.75*(4.-5.*sini2) );
  // contributions from J2^2:
  //  omega += odata->dt*odata->n*J22*Rep4*( 3./128.*(-10.-224.*cos2o*e2*cosi2+270.*cos2o*e2*cosi4-64.*cos2o*cosi2+10.*cos2o*e2+430.*cosi4+60.*cos2o*cosi4+126.*e2*cosi2-45.*e2*cosi4-192.*some2*cosi2-25.*e2+360.*some2*cosi4-36.*cosi2+24.*some2+4.*cos2o) );
  omega += odata->dt*odata->n*J22*Rep4*3./16.*(48.-103.*sini2+215./4.*sini4+(7.-4.5*sini2-45./8.*sini4)*e2+6.*some2*(1.-1.5*sini2)*(4.-5.*sini2)-0.25*(2.*(14.-15.*sini2)*sini2-(28.-158.*sini2+135.*sini4)*e2)*cos2o);
  // contributions from J3:
  omega += odata->dt*odata->n*J3*Rep3*( 3./8.*(
					       (4.-5.*sini2)*(sini2-e2*cosi2)/(odata->e*sini)  // "bad guy"
					       +2.*odata->e*sini*(13.-15.*sini2)
					       )*sino );
  // contributions from J4:
  //  omega += odata->dt*odata->n*J4*Rep4*( 15./128.*(-12.+144.*cosi2+4.*cos2o+126.*e2*cosi2-189.*e2*cosi4-32.*cos2o*cosi2+28.*cos2o*cosi4-196.*cosi4+10.*cos2o*e2-9.*e2-112.*cos2o*e2*cosi2+126.*cos2o*e2*cosi4) );
  omega += odata->dt*odata->n*J4*Rep4*15./32.*(16.-62.*sini2+49.*sini4+0.75*(24.-84.*sini2+63.*sini4)*e2+(sini2*(6.-7.*sini2)-0.5*(12.-70.*sini2+63.*sini4)*e2)*cos2o);
#endif


#ifdef ORBIT_TRANS_ELEM
  // orbital elements for small eccentricities
  // ep ( \dot{e} without J3 contributions)
  long double ep = -3./32.*J22*Rep4*sini2*(14.-15.*sini2)*odata->e*ome2*sin2o 
    -15./32.*J4*Rep4*sini2*(6.-7.*sini2)*odata->e*ome2*sin2o;
  // op ( \dot{omega} without J3 contributions)
  long double op = 0.75*J2*Rep2*(4.-5.*sini2)
    +3./16.*J22*Rep4*(48.-103.*sini2+215./4.*sini4+(7.-4.5*sini2-45./8.*sini4)*e2+6.*some2*(1.-1.5*sini2)*(4.-5.*sini2)-0.25*(2.*(14.-15.*sini2)*sini2-(28.-158.*sini2+135.*sini4)*e2)*cos2o) 
    -15./32.*J4*Rep4*(16.-62.*sini2+49.*sini4+0.75*(24.-84.*sini2+63.*sini4)*e2+(sini2*(6.-7.*sini2)-0.5*(12.-70.*sini2+63.*sini4)*e2)*cos2o);

  // h = e * cos(omega)
  h = odata->e*coso;
  // hp = ep * cos(omega) - e * sin(omega) * op
  h += odata->dt*odata->n*(ep*coso-odata->e*sino*op);
  // J3 contribution
  h += odata->dt*odata->n*( -3./8.*J3*Rep3*(-1.-4.*e2+5.*coso2*e2+35.*e2*cosi2-41.*coso2*e2*cosi2+40.*coso2*e2*cosi4-5.*cosi4-35.*e2*cosi4+6.*cosi2)/sini );

  // k = e * sin(omega)
  k = odata->e*sino;
  // kp = ep * sin(omega) + e * cos(omega) * op
  k += odata->dt*odata->n*(ep*sino+odata->e*coso*op);
  // J3 contribution
  k += odata->dt*odata->n*( -3./8.*coso*sino*(5.+40.*cosi4-41.*cosi2)*J3*Rep3*e2/sini );
#endif

#ifndef ORBIT_TRANS_ELEM
  // mean anomaly:
  M = odata->M;
  M += odata->dt*odata->n*((1.+1.5*J2*Rep2*(1-1.5*sini2)*some2)
			   +1.5*J22*Rep4*( powl(1.-1.5*sini2,2.)*ome2+(1.25*(1.-2.5*sini2+13./8.*sini4)+0.625*(1.-sini2+0.625*sini4)*e2+1./16.*sini2*(14.-15.*sini2)*(1.-2.5*e2)*cos2o)*some2 + 0.25*(3.*(3.-7.5*sini2+47./8.*sini4+(1.5-5.*sini2+117./6.*sini4)*e2-0.125*(1.+5.*sini2-101./8.*sini4)*powl(odata->e,4.))+e2/8.*sini2*(70.-123.*sini2+(56.-66.*sini2)*e2)*cos2o+27./128.*powl(odata->e,4.)*sini4*cosl(4.*odata->omega))/some2 )
			   -3./8.*J3*Rep3*sini*(4.-5.*sini2)*(1.-4.*e2)*some2*sino/odata->e
			   -45./128.*J4*Rep4*((8.-40.*sini2+35.*sini4)*e2*some2+2./3.*sini2*(6.-7.*sini2)*(2.-5.*e2)*some2*cos2o)
			   );
#endif

#ifdef ORBIT_TRANS_ELEM
  // u = omega + M
  u = odata->omega + odata->M;
                  // contribution of d\omega/dt
  u += odata->dt*odata->n*(+0.75*J2*Rep2*(4.-5.*sini2)
			   +3./16.*J22*Rep4*(48.-103.*sini2+215./4.*sini4+(7.-4.5*sini2-45./8.*sini4)*e2+6.*some2*(1.-1.5*sini2)*(4.-5.*sini2)-0.25*(2.*(14.-15.*sini2)*sini2-(28.-158.*sini2+135.*sini4)*e2)*cos2o)
			   -15./32.*J4*Rep4*(16.-62.*sini2+49.*sini4+0.75*(24.-84.*sini2+63.*sini4)*e2+(sini2*(6.-7.*sini2)-0.5*(12.-70.*sini2+63.*sini4)*e2)*cos2o) 
		  // contribution of dM/dt
			   +(1.+1.5*J2*Rep2*(1.-1.5*sini2)*some2)
			   +1.5*J22*Rep4*( powl(1.-1.5*sini2,2.)*ome2+(1.25*(1.-2.5*sini2+13./8.*sini4)+0.625*(1.-sini2+0.625*sini4)*e2+1./16.*sini2*(14.-15.*sini2)*(1.-2.5*e2)*cos2o)*some2 + 0.25*(3.*(3.-7.5*sini2+47./8.*sini4+(1.5-5.*sini2+117./6.*sini4)*e2-0.125*(1.+5.*sini2-101./8.*sini4)*powl(odata->e,4.))+e2/8.*sini2*(70.-123.*sini2+(56.-66.*sini2)*e2)*cos2o+27./128.*powl(odata->e,4.)*sini4*cosl(4.*odata->omega))/some2 )
			   -45./128.*J4*Rep4*((8.-40.*sini2+35.*sini4)*e2*some2+2./3.*sini2*(6.-7.*sini2)*(2.-5.*e2)*some2*cos2o) 
		  // contribution of J_3 term of (d\omega + dM)/dt
			   +J3*Rep3*(-3./8.*sino*(35.*cosi4*e2+5.*cosi4-5.*cosi4*some2+20.*cosi4*some2*e2-6.*cosi2+6.*cosi2*some2-35.*e2*cosi2-24.*cosi2*some2*e2+1.-some2+4.*some2*e2+4.*e2)/(odata->e*sini))
			   );
#endif
  

#ifdef ORBIT_TRANS_ELEM
  // recalculate Keplerian orbital elements
  e = sqrtl(powl(h,2.)+powl(k,2.));
  omega = atan2l(k, h);
  M = u - omega;
#endif


  // keep all angles in a range [0°:180°] and [0°:360°] respectively:
  if (i < 0.) {
    i += M_PI;
  } else if (i > M_PI) {
    i -= M_PI;
  }
  if (Omega < 0.) {
    Omega += 2.*M_PI;
  } else if (Omega > 2.*M_PI) {
    Omega -= 2.*M_PI;
  }
  if (omega < 0.) {
    omega += 2.*M_PI;
  } else if (omega > 2.*M_PI) {
    omega -= 2.*M_PI;
  }
  if (M < 0.) {
    M += 2.*M_PI;
  } else if (M > 2.*M_PI) {
    M -= 2.*M_PI;
  }

  
  // save perturbed Keplerian orbital elements to data structure:
  odata->a = a;
  odata->e = e;
  odata->i = i;
  odata->Omega = Omega;
  odata->omega = omega;
  odata->M = M;



////////////////////////////////////////////////////
// calculate and add short-period perturbation terms


  // calculate required data:
  // frequently used values:
  e2 = powl(e,2.);
  long double e3 = powl(e,3.);
  sini2 = powl(sinl(i),2.);
  ome2 = 1.-e2;
  some2 = sqrtl(ome2);

  // approximation of eccentric anomaly:
  E = M + (e+e3/8.)*sinl(M) + 0.5*e2*sinl(2.*M) + 0.375*e3*sin(3.*M);
  // true anomaly:
  f = 2.*atanl( sqrtl((1.+e)/(1.-e))*tanl(0.5*E) );

  // remove multiples of 2 \pi (because 'f-M' is used in the following formulas):
  while (f-M < -M_PI) { 
    f+=2.*M_PI; 
  }
  while (f-M > M_PI) {
    f-=2.*M_PI;
  }

  // calculate ellipse data:
  r = a*(1.-e2)/(1.+e*cosl(f));
  p = a*(1.-e2);

  // further frequently used values, depending on p:
  Rep2 = powl(R_e/p,2.);







  // a - semimajor axis:
  a += J2*a*powl(R_e/a,2.)*( powl(a/r,3.)*((1.-1.5*sini2)+1.5*sini2*cos(2.*omega+2.*f))-(1.-1.5*sini2)/pow(some2,3.) );

  // e - eccentricity:
  e += /** 0.5*J2*Rep2*(1.-1.5*sini2)*( 1./e*(1.+1.5*e2-pow(some2,-3.))+3.*(1.+e2/4.)*cos(f)+1.5*e*cos(2.*f)+e2/4.*cos(3.*f)) */
    +0.375*J2*Rep2*sini2*((1.+2.75*e2)*cos(2.*omega+f)+e2/4.*cos(2.*omega-f)+5.*e*cos(2.*(omega+f))+(7.+4.25*e2)/3.*cos(2.*omega+3.*f)+1.5*e*cos(2.*omega+4.*f)+e2/4.*cos(2.*omega+5.*f)+1.5*e*cos(2.*omega));

  // i - inclination:
  i += 0.375*J2*Rep2*sinl(2.*i)*(e*cosl(2.*omega+f)+cosl(2.*(omega+f))+e/3.*cosl(2.*omega+3.*f));

  // Omega - right ascension of ascending node:
  Omega += -1.5*J2*Rep2*cosl(i)*(f-M
				 +e*sinl(f)-0.5*e*sinl(2.*omega+f)-0.5*sinl(2.*(omega+f))-e*sinl(2.*omega+3.*f)/6.);


  // omega - argument of perigee:
  omega += 0.75*J2*Rep2*(4.-5.*sini2)*(f-M+e*sin(f))
    /* +1.5*J2*Rep2*(1.-1.5*sini2)*((1.-0.25*e2)/e*sin(f)+0.5*sin(2.*f)+e/12.*sin(3.*f))
       -1.5*J2*Rep2*(1./e*(0.25*sini2+0.5*e2*(1.-15./8.*sini2))*sin(2.*omega+f)+e2/16.*sini2*sin(2.*omega-f)-0.5*(1.-2.5*sini2)*sin(2.*(omega+f))-1./e*(7./12.*sini2-e2/6.*(1.-2.375*sini2))*sin(2.*omega+3.*f)-0.375*sini2*sin(2.*omega+4.*f)-1./16.*e*sini2*sin(2.*omega+5.*f)) */
        -9./16.*J2*Rep2*sini2*sin(2.*omega);

  // M - mean anomaly:
  M += /* -1.5*J2*Rep2*some2/e*((1.-1.5*sini2)*((1.-0.25*e2)*sin(f)+0.5*e*sin(2.*f)+e2/12.*sin(3.*f))+0.5*sini2*(-0.5*(1.+1.25*e2)*sin(2.*omega+f)-e2/8.*sin(2.*omega-f)+7./6.*(1.+e2/28.)*sin(2.*omega+3.*f)+0.75*e*sin(2.*omega+4.*f)+e2/8.*sin(2.*omega+5.*f))) */
    +9./16.*J2*Rep2*some2*sini2*sin(2.*omega);


  


#ifdef ORBIT_TRANS_ELEM

  sino = sinl(omega);
  coso = cosl(omega);
  cosi2 = pow(cosl(odata->i),2.);
  

  // h = e * cos(omega)
  h = e*cosl(omega);
  // hp = ep * cos(omega) - e * sin(omega) * op
  h += J2*Rep2*(
		18.*e3 *sino*sinl(f)*cosi2+6.*coso*e3 * cosl(3.* f) *cosi2+ 8.*coso*some2- 8.*coso- 36.*e2*sino* sinl(2.*omega+ 2.* f) + 24.* e *sino*sinl(f)- 0.11e2 *e3 *sino* sinl(2.*omega+ 3.* f) - 0.21e2 *e3 *sino* sinl(2.*omega+ f) - 0.28e2 * e *sino* sinl(2.*omega+ 3.* f) + 12.* e *sino* sinl(2.*omega+ f) - 12.*coso*e2- 8.*coso*some2*e2- 24.*coso*some2*cosi2+ 24.*coso*cosi2+ 24.*coso*some2*cosi2*e2+72. *coso* cosl(f) * e *cosi2+ 18.*coso* cosl(f) *e3 *cosi2+ 0.28e2 * e *sino* sinl(2.*omega+ 3.* f) *cosi2+ 0.45e2 *e3 *sino* sinl(2.*omega+ f) *cosi2- 36.*e2*sino* sinl(2.* f) *cosi2-6.*e3 *sino* sinl(3.* f) *cosi2- 12.* e *sino* sinl(2.*omega+ f) *cosi2- 18.*e2*sino* sinl(2.*omega+4. * f) - 3.*e3 *sino* sinl(2.*omega+5. * f) - 3.*e2*e2*sino* sinl(-2.*omega+ f) + 36.*coso*e2*cosi2+ 0.19e2 *e3 *sino* sinl(2.*omega+ 3.* f) *cosi2+ 18.*e2*sino* sinl(2.*omega+4. * f) *cosi2+ 3.*e3 *sino* sinl(2.*omega+5. * f) *cosi2+ 60. *e2*sino* sinl(2.*omega+ 2.* f) *cosi2+ 3.*e2*e2*sino* sinl(-2.*omega+ f) *cosi2+ 36.*coso*e2* cosl(2.* f) *cosi2-72. * e *sino*sinl(f)*cosi2- 24.*coso* cosl(f) * e -6.*coso* cosl(f) *e3 - 12.*coso*e2* cosl(2.* f) - 2.*coso*e3 * cosl(3.* f) + 12.*e2*sino* sinl(2.* f) + 2.*e3 *sino* sinl(3.* f) -6.*e3 *sino* sinl(f)) / e / 32.;


  // k = e * sin(omega)
  k = e*sinl(omega);
  // kp = ep * sin(omega) + e * cos(omega) * op
  k += J2*Rep2*(
		18.*e3 *sino*sinl(f)*cosi2+6.*coso*e3 * cosl(3.* f) *cosi2+ 8.*coso*some2- 8.*coso- 36.*e2*sino* sinl(2.*omega+ 2.* f) + 24.* e *sino*sinl(f)- 0.11e2 *e3 *sino* sinl(2.*omega+ 3.* f) - 0.21e2 *e3 *sino* sinl(2.*omega+ f) - 0.28e2 * e *sino* sinl(2.*omega+ 3.* f) + 12.* e *sino* sinl(2.*omega+ f) - 12.*coso*e2- 8.*coso*some2*e2- 24.*coso*some2*cosi2+ 24.*coso*cosi2+ 24.*coso*some2*cosi2*e2+72. *coso* cosl(f) * e *cosi2+ 18.*coso* cosl(f) *e3 *cosi2+ 0.28e2 * e *sino* sinl(2.*omega+ 3.* f) *cosi2+ 0.45e2 *e3 *sino* sinl(2.*omega+ f) *cosi2- 36.*e2*sino* sinl(2.* f) *cosi2-6.*e3 *sino* sinl(3.* f) *cosi2- 12.* e *sino* sinl(2.*omega+ f) *cosi2- 18.*e2*sino* sinl(2.*omega+4. * f) - 3.*e3 *sino* sinl(2.*omega+5. * f) - 3.*e2*e2*sino* sinl(-2.*omega+ f) + 36.*coso*e2*cosi2+ 0.19e2 *e3 *sino* sinl(2.*omega+ 3.* f) *cosi2+ 18.*e2*sino* sinl(2.*omega+4. * f) *cosi2+ 3.*e3 *sino* sinl(2.*omega+5. * f) *cosi2+ 60. *e2*sino* sinl(2.*omega+ 2.* f) *cosi2+ 3.*e2*e2*sino* sinl(-2.*omega+ f) *cosi2+ 36.*coso*e2* cosl(2.* f) *cosi2-72. * e *sino*sinl(f)*cosi2- 24.*coso* cosl(f) * e -6.*coso* cosl(f) *e3 - 12.*coso*e2* cosl(2.* f) - 2.*coso*e3 * cosl(3.* f) + 12.*e2*sino* sinl(2.* f) + 2.*e3 *sino* sinl(3.* f) -6.*e3 *sino* sinl(f)) / e / 32.;



  // u = omega + M
  u = omega + M;

  // contributions of 1/e terms in Delta omega and Delta M
  u += J2*Rep2*(
		-(24.*sinl(f)+ 12.*sinl(2.*omega+ f) - 0.28e2 *sinl(2.*omega+ 3.* f) + 12.* e * sinl(2.* f) + 2.*e2* sinl(3.* f) - 18.* e *sinl(2.*omega+4. * f) - 3.*e2*sinl(2.*omega+5. * f) + 0.19e2 *sinl(2.*omega+ 3.* f) *e2*cosi2+ 18.*sinl(f)*e2*cosi2+ 12.*some2*sinl(2.*omega+ f) *cosi2-6.*e2* sinl(3.* f) *cosi2- 12.*some2* e * sinl(2.* f) - 2.*some2*e2* sinl(3.* f) +6.*some2*e2* sinl(3.* f) *cosi2- 3.*some2*e2* sinl(-2.*omega+ f) *cosi2- 18.*some2* e *sinl(2.*omega+4. * f) *cosi2+ 15.*some2*sinl(2.*omega+ f) *e2*cosi2- 18.*some2*sinl(f)*e2*cosi2-6.*sinl(f)*e2- 36.*sinl(2.*omega+ 2.* f) * e - 0.21e2 *sinl(2.*omega+ f) *e2- 0.11e2 *sinl(2.*omega+ 3.* f) *e2- 24.*some2*sinl(f)+ 3.*e2*sinl(2.*omega+5. * f) *cosi2+ 18.* e *sinl(2.*omega+4. * f) *cosi2+ 60. *sinl(2.*omega+ 2.* f) * e *cosi2+ 3.*powl(e,3.)* sinl(-2.*omega+ f) *cosi2+ 0.28e2 *sinl(2.*omega+ 3.* f) *cosi2-72. *sinl(f)*cosi2- 12.*sinl(2.*omega+ f) *cosi2+6.*some2*sinl(f)*e2- 3.*powl(e,3.)* sinl(-2.*omega+ f) + 3.*some2*e2* sinl(-2.*omega+ f) + 18.*some2* e *sinl(2.*omega+4. * f) + 3.*some2*e2*sinl(2.*omega+5. * f) - 15.*some2*sinl(2.*omega+ f) *e2-some2*sinl(2.*omega+ 3.* f) *e2- 12.*some2*sinl(2.*omega+ f) + 0.28e2 *some2*sinl(2.*omega+ 3.* f) +72. *some2*sinl(f)*cosi2- 0.28e2 *some2*sinl(2.*omega+ 3.* f) *cosi2- 36.* e * sinl(2.* f) *cosi2+ 0.45e2 *sinl(2.*omega+ f) *e2*cosi2+some2*sinl(2.*omega+ 3.* f) *e2*cosi2- 3.*some2*e2*sinl(2.*omega+5. * f) *cosi2+ 36.*some2* e * sinl(2.* f) *cosi2)) / e / 32.;


  // recalculate Keplerian orbital elements
  e = sqrtl(powl(h,2.)+powl(k,2.));
  omega = atan2l(k, h);
  M = u - omega;
#endif







/////////////////////////////////////////////////////////////////////
// calculate the position and velocity  of the satellite on the orbit

//printf("%Lf %Lf %Lf\n", oposition->t, omega*180./M_PI, M*180./M_PI);

  // get the eccentric anomaly from the mean anomaly:
  E = kepler_equation(M,e);
  // true anomaly:
  f = 2.*atanl( sqrtl((1.+e)/(1.-e))*tanl(0.5*E) );
  // parameter of the ellipse:
  p = a*(1.-powl(e,2.));

  // calculation acceleration constants:
  long double sinO = sinl(Omega);
  long double cosO = cosl(Omega);
  cosi = cosl(i);
  sini = sinl(i);

  // now we know the argument of latitude
  // and can calculate its cosine and sine
  long double cosu = cosl(omega+f);
  long double sinu = sinl(omega+f);

  // now we can determine the position and velocity of the satellite (Flury p. 38)
  long double rl = p/(1.+e*cosl(f));
  long double vr = sqrtl(mu/p)*e*sinl(f);
  long double vf = sqrtl(mu*p)/rl;
  
  // position
  oposition->r.x = rl*(cosO*cosu-sinO*cosi*sinu);
  oposition->r.y = rl*(sinO*cosu+cosO*cosi*sinu);
  oposition->r.z = rl*sini*sinu;

  // velocity
  oposition->v.x = vr*(cosO*cosu-sinO*cosi*sinu) - vf*(cosO*sinu+sinO*cosu*cosi);
  oposition->v.y = vr*(sinO*cosu+cosO*cosi*sinu) - vf*(sinO*sinu-cosO*cosu*cosi);
  oposition->v.z = vr*sini*sinu + vf*cosu*sini;

/*
  // calculate the position and velocity according to Flury p.20:
  sino = sinl(omega);
  coso = cosl(omega);
  sinO = sinl(Omega);
  cosO = cosl(Omega);
  sini = sinl(i);
  cosi = cosl(i);
  long double C1 = a*(-e+cosl(E));
  long double C2 = a*sqrtl(1.-powl(e,2.))*sinl(E);
  long double C3 = -sqrtl(mu*a)/r*sinl(E);
  long double C4 = sqrtl(mu*p)/r*cosl(E);
  struct vector ex = { coso*cosO-sino*sinO*cosi,  coso*sinO+sino*cosO*cosi, sino*sini};
  struct vector ey = {-sino*cosO-coso*sinO*cosi, -sino*sinO+coso*cosO*cosi, coso*sini};

  oposition->r.x = C1*ex.x + C2*ey.x;
  oposition->r.y = C1*ex.y + C2*ey.y;
  oposition->r.z = C1*ex.z + C2*ey.z;

  oposition->v.x = C3*ex.x + C4*ey.x;
  oposition->v.y = C3*ex.y + C4*ey.y;
  oposition->v.z = C3*ex.z + C4*ey.z;


  // calculate the position according to Steiner:
  //  oposition->r.x = p/(1.+e*cos(f)) * (cos(f)*(cos(omega)*cos(Omega)-sin(omega)*sin(Omega)*cos(i)) + sin(f)*(-sin(omega)*cos(Omega)-cos(omega)*sin(Omega)*cos(i)));
  //  oposition->r.y = p/(1.+e*cos(f)) * (cos(f)*(cos(omega)*sin(Omega)+sin(omega)*cos(Omega)*cos(i)) + sin(f)*(-sin(omega)*sin(Omega)+cos(omega)*cos(Omega)*cos(i)));
  //  oposition->r.z = p/(1.+e*cos(f)) * (cos(f)*sin(omega)*sin(i) + sin(f)*cos(omega)*sin(i));
*/

}








//////////////////////////////////////////////////////////////////////////////////////////////
// This function performs a timestep, i.e. it calculates the next position of the satellite //
// on its orbit, considering higher order perturbation terms J2, J2^2, J3, J4.              //
// The calculation used follows the approach of Lyddane.                                    //
//////////////////////////////////////////////////////////////////////////////////////////////
void orbit_step_J234t(struct orbit_data *odata, struct orbit_position *oposition) 
{
  const long double J22 = powl(J2,2.);

  // orbital elements:
  long double a;        // semimajor axis
  long double e;        // eccentricity
  long double i;        // inclination
  long double Omega;    // right ascension of ascending node
  long double omega;    // argument of perigee
  long double M;        // mean anomaly

  // elements for low-eccentricity orbits
  long double h,k,u;
  
  long double p;                   // parameter of ellipse
  long double E;                   // eccentric anomaly
  long double f;                   // true anomaly
  long double r;                   // lenght of position vector |r|

  // frequently used values for the calculation:
  long double Rep2, Rep3, Rep4;    // (R_e/p)^2, (R_e/p)^3, (R_e/p)^4
  long double sini, cosi;          // sin(i), cos(i)
  long double sini2;               // sin^2(i)
  long double sini4;               // sin^4(i)
  long double cosi2;               // cos^2(i)
  long double cosi4;               // cos^4(i)
  long double sino, coso;          // sin(omega), cos(omega)
  long double sin2o, cos2o;        // sin(2*omega), cos(2*omega)
  long double sino2, coso2;        // sin^2(omega), cos^2(omega)
  long double sinO, cosO;          // sin(Omega), cos(Omega)
  long double e2, e3;              // e^2, e^3
  long double ome2;                // 1-e^2
  long double some2;               // sqrt(1-e^2)



//////////////////////////////////////////////////////////
// change of mean orbital elements (secular perturbations)
  // required data:
  p = odata->a*(1.-powl(odata->e,2.));

  // frequently used values:
  Rep2 = powl(R_e/p,2.);
  Rep3 = powl(R_e/p,3.);
  Rep4 = powl(R_e/p,4.);
  sini = sinl(odata->i);
  sini2 = powl(sinl(odata->i),2.);
  sini4 = powl(sinl(odata->i),4.);
  cosi2 = powl(cosl(odata->i),2.);
  cosi4 = powl(cosl(odata->i),4.);
  sino = sinl(odata->omega);
  coso = cosl(odata->omega);
  sin2o = sinl(2.*odata->omega);
  cos2o = cosl(2.*odata->omega);
  sino2 = powl(sinl(odata->omega),2.);
  coso2 = powl(cosl(odata->omega),2.);
  e2 = powl(odata->e,2.);
  ome2 = 1.-e2;
  some2 = sqrtl(ome2);


  // semimajor axis: 
  long double dat = odata->dt*0.;  // a doesn't change (\dot{a}=0)
  a = odata->a + dat;
  
  // inclination:
  long double dit = odata->dt*( -3./64.*odata->n*J22*Rep4*sinl(2*odata->i)*(14.-15.*sini2)*e2*sin2o
				+0.375*odata->n*J3*Rep3*cosl(odata->i)*(4.-5.*sini2)*odata->e*cosl(odata->omega)
				+15./64.*odata->n*J4*powl(R_e/p,4.)*sinl(2.*odata->i)*(6.-7.*sini2)*e2*sin2o );
  i = odata->i + dit;

  // right ascension of ascending node:
  long double dOmegat = odata->dt*( -1.5*odata->n*J2*Rep2*cosl(odata->i) 
			-1.5*odata->n*J22*Rep4*cosl(odata->i)*(2.25+1.5*some2-sini2*(2.5+2.25*some2)+e2/16.*(4.+5.*sini2) - e2*0.125*(7.-15.*sini2)*cos2o)
			+0.375*odata->n*J3*Rep3*(4.-5.*sini2)*odata->e/tanl(odata->i)*sinl(odata->omega)
			+15./16.*odata->n*J4*Rep4*cos(odata->i)*((4.-7.*sini2)*(1.+1.5*e2)-(3.-7.*sini2)*e2*cos2o) );
  Omega = odata->Omega + dOmegat;

  // orbital elements for small eccentricity
  // h = e * cos(omega)
  h = odata->e*cosl(odata->omega);
  // hp = ep * cos(omega) - e * sin(omega) * op
  h += odata->dt*((-3./32.*odata->n*J22*Rep4*sini2*(14.-15.*sini2)*odata->e*ome2*sin2o 
		   -15./32.*odata->n*J4*Rep4*sini2*(6.-7.*sini2)*odata->e*ome2*sin2o)*cosl(odata->omega));
  h += odata->dt*(-odata->e*sinl(odata->omega)*(0.75*odata->n*J2*Rep2*(4.-5.*sini2)));
  h += odata->dt*(-odata->e*sinl(odata->omega)*(+3./16.*odata->n*J22*Rep4*(48.-103.*sini2+215./4.*sini4+(7.-4.5*sini2-45./8.*sini4)*e2+6.*some2*(1.-1.5*sini2)*(4.-5.*sini2)-0.25*(2.*(14.-15.*sini2)*sini2-(28.-158.*sini2+135.*sini4)*e2)*cos2o) 
						-15./32.*odata->n*J4*Rep4*(16.-62.*sini2+49.*sini4+0.75*(24.-84.*sini2+63.*sini4)*e2+(sini2*(6.-7.*sini2)-0.5*(12.-70.*sini2+63.*sini4)*e2)*cos2o) ) );
  // old (e=0): h += odata->dt*(+0.375*odata->n*J3*Rep3/sinl(odata->i)*(5.*cosi4 +1.-6.*cosi2) );
  h += odata->dt*(-3./8.*odata->n*J3*Rep3*(-1.-4.*e2+5.*coso2*e2+35.*e2*cosi2-41.*coso2*e2*cosi2+40.*coso2*e2*cosi4-5.*cosi4-35.*e2*cosi4+6.*cosi2)/sini );

  // k = e * sin(omega)
  k = odata->e*sinl(odata->omega);
  // kp = ep * sin(omega) + e * cos(omega) * op
  k += odata->dt*((-3./32.*odata->n*J22*Rep4*sini2*(14.-15.*sini2)*odata->e*ome2*sin2o 
		   -15./32.*odata->n*J4*Rep4*sini2*(6.-7.*sini2)*odata->e*ome2*sin2o)*sinl(odata->omega));
  k += odata->dt*(+odata->e*cosl(odata->omega)*(0.75*odata->n*J2*Rep2*(4.-5.*sini2)));
  k += odata->dt*(+odata->e*cosl(odata->omega)*(3./16.*odata->n*J22*Rep4*(48.-103.*sini2+215./4.*sini4+(7.-4.5*sini2-45./8.*sini4)*e2+6.*some2*(1.-1.5*sini2)*(4.-5.*sini2)-0.25*(2.*(14.-15.*sini2)*sini2-(28.-158.*sini2+135.*sini4)*e2)*cos2o) 
						-15./32.*odata->n*J4*Rep4*(16.-62.*sini2+49.*sini4+0.75*(24.-84.*sini2+63.*sini4)*e2+(sini2*(6.-7.*sini2)-0.5*(12.-70.*sini2+63.*sini4)*e2)*cos2o)) );
  k += odata->dt*( -3./8.*coso*sino*(5.+40.*cosi4-41.*cosi2)*odata->n*J3*Rep3*e2/sini );

  // u = omega + M
  u = odata->omega + odata->M;
                  // contribution of d\omega/dt
  u += odata->dt*(+0.75*odata->n*J2*Rep2*(4.-5.*sini2)
		  +3./16.*odata->n*J22*Rep4*(48.-103.*sini2+215./4.*sini4+(7.-4.5*sini2-45./8.*sini4)*e2+6.*some2*(1.-1.5*sini2)*(4.-5.*sini2)-0.25*(2.*(14.-15.*sini2)*sini2-(28.-158.*sini2+135.*sini4)*e2)*cos2o)
		  -15./32.*odata->n*J4*Rep4*(16.-62.*sini2+49.*sini4+0.75*(24.-84.*sini2+63.*sini4)*e2+(sini2*(6.-7.*sini2)-0.5*(12.-70.*sini2+63.*sini4)*e2)*cos2o) 
		  // contribution of dM/dt
		  +odata->n*(1.+1.5*J2*Rep2*(1.-1.5*sini2)*some2)
		  +1.5*odata->n*J22*Rep4*( powl(1.-1.5*sini2,2.)*ome2+(1.25*(1.-2.5*sini2+13./8.*sini4)+0.625*(1.-sini2+0.625*sini4)*e2+1./16.*sini2*(14.-15.*sini2)*(1.-2.5*e2)*cos2o)*some2 + 0.25*(3.*(3.-7.5*sini2+47./8.*sini4+(1.5-5.*sini2+117./6.*sini4)*e2-0.125*(1.+5.*sini2-101./8.*sini4)*powl(odata->e,4.))+e2/8.*sini2*(70.-123.*sini2+(56.-66.*sini2)*e2)*cos2o+27./128.*powl(odata->e,4.)*sini4*cosl(4.*odata->omega))/some2 )
		  +45./128.*odata->n*J4*Rep4*((8.-40.*sini2+35.*sini4)*e2*some2+2./3.*sini2*(6.-7.*sini2)*(2.-5.*e2)*some2*cos2o) 
		  // contribution of J_3 term of (d\omega + dM)/dt
		  +odata->n*J3*Rep3*(-3./8.*sino*(35.*cosi4*e2+5.*cosi4-5.*cosi4*some2+20.*cosi4*some2*e2-6.*cosi2+6.*cosi2*some2-35.*e2*cosi2-24.*cosi2*some2*e2+1.-some2+4.*some2*e2+4.*e2)/(odata->e*sinl(odata->i)))
		  );

  

  // recalculate Keplerian orbital elements
  e = sqrtl(powl(h,2.)+powl(k,2.));
  omega = atan2l(k, h);
  M = u - omega;


  // keep all angles in a range [0°:180°] and [0°:360°] respectively:
  if (i < 0.) {
    i += M_PI;
  } else if (i > M_PI) {
    i -= M_PI;
  }
  if (Omega < 0.) {
    Omega += 2.*M_PI;
  } else if (Omega > 2.*M_PI) {
    Omega -= 2.*M_PI;
  }
  if (omega < 0.) {
    omega += 2.*M_PI;
  } else if (omega > 2.*M_PI) {
    omega -= 2.*M_PI;
  }
  if (M < 0.) {
    M += 2.*M_PI;
  } else if (M > 2.*M_PI) {
    M -= 2.*M_PI;
  }

  
  // save perturbed Keplerian orbital elements to data structure:
  odata->a = a;
  odata->e = e;
  odata->i = i;
  odata->Omega = Omega;
  odata->omega = omega;
  odata->M = M;


  

////////////////////////////////////////////////////
// calculate and add short-period perturbation terms


  // calculate required data:
  // frequently used values:
  e2 = powl(e,2.);
  e3 = powl(e,3.);
  sini2 = powl(sinl(i),2.);
  ome2 = 1.-e2;
  some2 = sqrtl(ome2);

  // approximation of eccentric anomaly:
  E = M + (e+e3/8.)*sinl(M) + 0.5*e2*sinl(2.*M) + 0.375*e3*sin(3.*M);
  // true anomaly:
  f = 2.*atanl( sqrtl((1.+e)/(1.-e))*tanl(0.5*E) );

  // remove multiples of 2 \pi (because 'f-M' is used in the following formulas):
  if (f-M < 0.) { 
    M-=2.*M_PI; 
  }

  // calculate ellipse data:
  r = a*(1.-e2)/(1.+e*cosl(f));
  p = a*(1.-e2);

  // frequently used values, depending on p:
  Rep2 = powl(R_e/p,2.);


  // a - semimajor axis:
  a += J2*a*powl(R_e/a,2.)*( powl(a/r,3.)*((1.-1.5*sini2)+1.5*sini2*cos(2.*omega+2.*f))-(1.-1.5*sini2)/pow(some2,3.) );

  /*  // e - eccentricity:
  e += 0.5*J2*Rep2*(1.-1.5*sini2)*( 1./e*(1.+1.5*e2-pow(some2,-3.))+3.*(1.+e2/4.)*cos(f)+1.5*e*cos(2.*f)+e2/4.*cos(3.*f) )
    +0.375*J2*Rep2*sini2*((1.+2.75*e2)*cos(2.*omega+f)+e2/4.*cos(2.*omega-f)+5.*e*cos(2.*(omega+f))+(7.+4.25*e2)/3.*cos(2.*omega+3.*f)+1.5*e*cos(2.*omega+4.*f)+e2/4.*cos(2.*omega+5.*f)+1.5*e*cos(2.*omega));
  */

  // i - inclination:
  i += 0.375*J2*Rep2*sini2*(e*cosl(2.*omega+f)+cosl(2.*(omega+f))+1.5*cosl(2.*omega+3.*f));

  // Omega - right ascension of ascending node:
  Omega += -1.5*J2*Rep2*cosl(i)*(f-M+e*sinl(f)-0.5*e*sinl(2.*omega+f)-0.5*sinl(2.*(omega+f))-e*sinl(2.*omega+3.*f)/6.);


  /*  // omega - argument of perigee:
  omega += 0.75*J2*Rep2*(4.-5.*sini2)*(f-M+e*sin(f))
    +1.5*J2*Rep2*(1.-1.5*sini2)*((1.-0.25*e2)/e*sin(f)+0.5*sin(2.*f)+e/12.*sin(3.*f))
    -1.5*J2*Rep2*(1./e*(0.25*sini2+0.5*e2*(1.-15./8.*sini2))*sin(2.*omega+f)+e2/16.*sini2*sin(2.*omega-f)-0.5*(1.-2.5*sini2)*sin(2.*(omega+f))-1./e*(7./12.*sini2-e2/6.*(1.-2.375*sini2))*sin(2.*omega+3.*f)-0.375*sini2*sin(2.*omega+4.*f)-1./16.*e*sini2*sin(2.*omega+5.*f))
    -9./16.*J2*Rep2*sini2*sin(2.*omega);

  // M - mean anomaly:
  M += -1.5*J2*Rep2*some2/e*((1.-1.5*sini2)*((1.-0.25*e2)*sin(f)+0.5*e*sin(2.*f)+e2/12.*sin(3.*f))+0.5*sini2*(-0.5*(1.+1.25*e2)*sin(2.*omega+f)-e2/8.*sin(2.*omega-f)+7./6.*(1.+e2/28.)*sin(2.*omega+3.*f)+0.75*e*sin(2.*omega+4.*f)+e2/8.*sin(2.*omega+5.*f)))
    +9./16.*J2*Rep2*some2*sini2*sin(2.*omega);
  */

  // perturbation terms for small eccentricity
  // h = e * cos(omega)
  h += ( 0.375*J2*Rep2*sini2*((1.+2.75*e2)*cosl(2.*omega+f)+e2/4.*cosl(2.*omega-f)+5.*e*cosl(2.*(omega+f))+(7.+4.25*e2)/3.*cosl(2.*omega+3.*f)+1.5*e*cosl(2.*omega+4.*f)+e2/4.*cosl(2.*omega+5.*f)+1.5*e*cosl(2.*omega)) )*cosl(omega)
    - e*sinl(omega)*( 0.75*J2*Rep2*(4.-5.*sini2)*(f-M+e*sinl(f)) -9./16.*J2*Rep2*sini2*sinl(2.*omega) )
    +Rep2*( -7./64.*J2*cosl(-3.*omega+2.*i-3.*f)-7./64.*J2*cosl(3.*omega+2.*i+3.*f)+3./64.*J2*cosl(-3.*omega-f+2.*i)+3./64.*J2*cosl(3.*omega+f+2.*i)+15./32.*J2*cosl(omega+f)+7./32.*J2*cosl(3.*f+3.*omega)-7./32.*J2*cosl(omega+3.*f)+33./64.*J2*cosl(2.*i-omega-f)+33./64.*J2*cosl(2.*i+omega+f)-3./32.*J2*cosl(f+3.*omega)+7./64.*J2*cosl(-3.*f-omega+2*i)+7./64.*J2*cosl(3.*f+omega+2.*i) );

  k += ( 0.375*J2*Rep2*sini2*((1.+2.75*e2)*cosl(2.*omega+f)+e2/4.*cosl(2.*omega-f)+5.*e*cosl(2.*(omega+f))+(7.+4.25*e2)/3.*cosl(2.*omega+3.*f)+1.5*e*cosl(2.*omega+4.*f)+e2/4.*cosl(2.*omega+5.*f)+1.5*e*cosl(2.*omega)) )*sinl(omega)
    + e*cosl(omega)*( 0.75*J2*Rep2*(4.-5.*sini2)*(f-M+e*sinl(f)) -9./16.*J2*Rep2*sini2*sinl(2.*omega) )
    +Rep2*( 7./32.*J2*sinl(3.*f+3.*omega)-3./64.*J2*sinl(-3.*omega-f+2.*i)+3./64.*J2*sinl(3.*omega+f+2.*i)-39./64.*J2*sinl(2.*i-omega-f)+39./64.*J2*sinl(2.*i+omega+f)-3./32.*J2*sinl(f+3.*omega)+9./32.*J2*sinl(omega+f)+7./32.*J2*sinl(omega+3.*f)+7./64.*J2*sinl(-3.*omega+2.*i-3.*f)-7./64.*J2*sinl(3.*omega+2.*i+3.*f)+7./64.*J2*sinl(-3.*f-omega+2.*i)-7./64.*J2*sinl(3.*f+omega+2.*i) );

  u += 0.75*J2*Rep2*(4.-5.*sini2)*(f-M+e*sinl(f)) 
           -9./16.*J2*Rep2*sini2*sinl(2.*omega) 
    +9./16.*J2*Rep2*some2*sini2*sinl(2.*omega)
    +Rep2*( 3./16.*J2*sinl(2.*omega+2.*f)-15./32.*J2*sinl(2.*f+2.*omega+2.*i)+15./32.*J2*sinl(-2.*f-2.*omega+2.*i) );



  // recalculate Keplerian orbital elements
  e = sqrtl(powl(h,2.0)+powl(k,2.0));
  omega = atan2l(k,h);
  M = u - omega;


/////////////////////////////////////////////////////////////////////
// calculate the position and velocity  of the satellite on the orbit

  // get the eccentric anomaly from the mean anomaly:
  E = kepler_equation(M,e);
  // true anomaly:
  f = 2.*atanl( sqrtl((1.+e)/(1.-e))*tanl(0.5*E) );
  // other orbital elements:
  p = a*(1.-powl(e,2.));

  // calculation acceleration constants:
  sinO = sinl(Omega);
  cosO = cosl(Omega);
  cosi = cosl(i);
  sini = sinl(i);

  // now we know the argument of latitude
  // and can calculate its cosine and sine
  long double cosu = cosl(omega+f);
  long double sinu = sinl(omega+f);


  // now we can determine the position and velocity of the satellite (Flury p. 38)
  long double rl = p/(1.+e*cosl(f));
  long double vr = sqrtl(mu/p)*e*sinl(f);
  long double vf = sqrtl(mu*p)/rl;
  
  // position
  oposition->r.x = rl*(cosO*cosu-sinO*cosi*sinu);
  oposition->r.y = rl*(sinO*cosu+cosO*cosi*sinu);
  oposition->r.z = rl*sini*sinu;

  // velocity
  oposition->v.x = vr*(cosO*cosu-sinO*cosi*sinu) - vf*(cosO*sinu+sinO*cosu*cosi);
  oposition->v.y = vr*(sinO*cosu+cosO*cosi*sinu) - vf*(sinO*sinu-cosO*cosu*cosi);
  oposition->v.z = vr*sini*sinu + vf*cosu*sini;

}








/////////////////////////////////////////////////////////////
// Solves the Kepler equation.                             //
// Calculates the eccentric anomaly from the mean anomaly. //
/////////////////////////////////////////////////////////////
long double kepler_equation(const long double M, const long double e) {
  long double E = M;
  long double dE;

  // Newton algorithm to solve the Kepler equation
  do {
    dE = (M - E + e*sinl(E))/(1.-e*cosl(E));
    E += dE;
  } while (fabsl(dE) > 0.00000001);

  return(E);
}





