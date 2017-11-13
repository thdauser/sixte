//
// Various integrators for the TES differential equations 
//

#include "tessim.h"

// Stochastic integrator
//
// This function evolves an autonomous SDE system
// dX = a(X)dt + b_1(X)dW_1 + ... + b_m(X)dW_m
// one step of size delta_t from a given state Y, 
// using an explicit strong Runge-Kutta algorithm of order 1.5.
// The new state of the system is then stored in Y.

int sde_step(double a(double X[], int k, void *params), double b(double X[], int k, int j, void *params), int dim, int m, double Y[], const double delta_t, const gsl_rng *r, void *params) {
	
	int p = 10; // the higher p the better the approximations of the stochastic integrals
	
	// allocate memory for the sampling points
	double *Y_plus = malloc(sizeof(double) * dim);
	if (Y_plus == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double *Y_minus = malloc(sizeof(double) * dim);
	if (Y_minus == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double *Phi_plus = malloc(sizeof(double) * dim);
	if (Phi_plus == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double *Phi_minus = malloc(sizeof(double) * dim);
	if (Phi_minus == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	// allocate memory for Wiener increments, auxiliary variables and matrices to store values of double and triple integrals 
	double *dW = malloc(sizeof(double) * m); 
	if (dW == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double *sol = malloc(sizeof(double) * dim);
	if (sol == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double *xi = malloc(sizeof(double) * m); 
	if (xi == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double **zeta = malloc(sizeof(double*) * m);
	if (zeta == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	for (int i = 0; i < m; i++) {
		zeta[i] = malloc(sizeof(double) * p);
		if (zeta[i] == NULL) {
			printf("ERROR: Out of memory\n");
			return EXIT_FAILURE;;
		}
	}
	
	double **eta = malloc(sizeof(double*) * m);
	if (eta == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	for (int i = 0; i < m; i++) {
		eta[i] = malloc(sizeof(double) * p);
		if (eta[i] == NULL) {
			printf("ERROR: Out of memory\n");
			return EXIT_FAILURE;;
		}
	}
	
	double *mu = malloc(sizeof(double) * m); 
	if (mu == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double *phi = malloc(sizeof(double) * m); 
	if (phi == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	
	double **II = malloc(sizeof(double*) * (m+1));
	if (II == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	for (int i = 0; i < m+1; i++) {
		II[i] = malloc(sizeof(double) * (m+1));
		if (II[i] == NULL) {
			printf("ERROR: Out of memory\n");
			return EXIT_FAILURE;;
		}
	}
	
	double ***III = malloc(sizeof(double**) * (m+1));
	if (III == NULL) {
		printf("ERROR: Out of memory\n");
		return EXIT_FAILURE;;
	}
	for (int i = 0; i < m+1; i++) {
		III[i] = malloc(sizeof(double*) * (m+1));
		if (III[i] == NULL) {
			printf("ERROR: Out of memory\n");
			return EXIT_FAILURE;;
		}
		for (int j = 0; j < m+1; j++) {
			III[i][j] = malloc(sizeof(double) * (m+1));
			if (III[i][j] == NULL) {
				printf("ERROR: Out of memory\n");
				return EXIT_FAILURE;;
			}
		}
	}
	
	// first initialize all random variables previously allocated
	for (int i = 0; i < m; i++) {
		xi[i] = gsl_ran_gaussian(r, 1);
		mu[i] = gsl_ran_gaussian(r, 1);
		phi[i] = gsl_ran_gaussian(r, 1);
		for (int j = 0; j < p; j++) {
			zeta[i][j] = gsl_ran_gaussian(r, 1);
			eta[i][j] = gsl_ran_gaussian(r, 1);
		}
	}
	
	// calculate the Wiener increments 
	for (int j = 0; j < m; j++) {
		dW[j] = sqrt(delta_t) * xi[j];
	}
	
	// calculate stochastic integrals II(0,0), II(j,0), II(0,j), II(j,j), III(j,j,j)
	II[0][0] = 0.5 * delta_t * delta_t; // elementary, not stochastic
	// coefficients rho_p, alpha_p. Only depend on p and are needed in following calculations
	double rho_p = 0, alpha_p = 0;
	for (int r = 1; r <= p; r++) {
		rho_p += 1./(r*r);
		alpha_p += 1./(r*r*r*r);
	}
	rho_p = 1./12 - rho_p/(2*M_PI*M_PI);
	alpha_p = M_PI*M_PI/180 - alpha_p/(2*M_PI*M_PI);
	
	for (int j = 0; j < m; j++) {
		double a_j0 = 0;
		for (int r = 1; r <= p; r++) {
			a_j0 += zeta[j][r-1]/r; 
		}
		a_j0 = - a_j0 * sqrt(2*delta_t)/M_PI - 2*sqrt(delta_t*rho_p)*mu[j];
		II[j+1][0] = 0.5 * delta_t * (sqrt(delta_t)*xi[j] + a_j0);
		II[0][j+1] = dW[j]*delta_t - II[j+1][0];
		II[j+1][j+1] = 0.5 * (dW[j]*dW[j]); //this is actually the Stratonovich integral because I need this for the triple integrals later. "Ito value" is inserted after calculation of triple integrals.
		III[j+1][j+1][j+1] = 0.5 * (dW[j]*dW[j]/3 - delta_t) * dW[j];
	}
	
	// calculate stochastic integrals II(j1,j2) and Stratonivich integrals J(j1,0,j2), J(0,j1,j2), J(j1,j2,0) (I need the Stratonivich integrals for III(j1,j2,j3))
	for (int j1 = 0; j1 < m; j1++) {
		
		// coefficients a_j1, b_j1
		double a_j1 = 0, b_j1 = 0;
		for (int r = 1; r <= p; r++) {
			a_j1 += zeta[j1][r-1]/r; 
			b_j1 += eta[j1][r-1]/(r*r);
		}
		a_j1 = - a_j1 * sqrt(2*delta_t)/M_PI - 2*sqrt(delta_t*rho_p)*mu[j1];
		b_j1 = b_j1*sqrt(delta_t/2) + sqrt(delta_t*alpha_p)*phi[j1];
		
		for (int j2 = 0; j2 < m; j2++) {
			
			// coefficients a_j2, b_j2
			double a_j2 = 0, b_j2 = 0;
			for (int r = 1; r <= p; r++) {
				a_j2 += zeta[j2][r-1]/r; 
				b_j2 += eta[j2][r-1]/(r*r);
			}
			a_j2 = - a_j2 * sqrt(2*delta_t)/M_PI - 2*sqrt(delta_t*rho_p)*mu[j2];
			b_j2 = b_j2*sqrt(delta_t/2) + sqrt(delta_t*alpha_p)*phi[j2];
			
			// calculate coefficients A_p, B_p, C_p
			double A_p = 0, B_p = 0, C_p = 0;
				for (int r = 1; r <= p; r++) {
					A_p += (zeta[j1][r-1] * eta[j2][r-1] - eta[j1][r-1] * zeta[j2][r-1])/r;
					B_p += (zeta[j1][r-1] * zeta[j2][r-1] + eta[j1][r-1] * eta[j2][r-1])/(r*r);
					for (int l = 1; l <= p; l++) {
						if (r != l) {
							C_p += (zeta[j1][r-1]*zeta[j2][l-1]/l - l*eta[j1][r-1]*eta[j2][l-1]/r) * (r/(1.*r*r-l*l));
						}
					}
				}
			A_p = A_p/(2*M_PI);
			B_p = B_p/(4*M_PI*M_PI);
			C_p = - C_p/(2*M_PI*M_PI);
			
			// now calculate II[j1][j2]
			if (j1 != j2) {
				II[j1+1][j2+1] = 0.5*delta_t*xi[j1]*xi[j2] - 0.5*sqrt(delta_t)*(xi[j1]*a_j2 - xi[j2]*a_j1) + A_p * delta_t;
			}
			
			// and Stratonivich integrals J(j1,0,j2), J(0,j1,j2), J(j1,j2,0)
			III[j1+1][0][j2+1] = delta_t*delta_t*xi[j1]*xi[j2]/6 + a_j1*II[0][j2+1]/2 + pow(delta_t,1.5)*xi[j2]*b_j1/(2*M_PI) - delta_t*delta_t*B_p 
								- pow(delta_t,1.5)*a_j2*xi[j1]/4 + pow(delta_t,1.5)*xi[j1]*b_j2/(2*M_PI);
			III[0][j1+1][j2+1] = delta_t*delta_t*xi[j1]*xi[j2]/6 - pow(delta_t,1.5)*xi[j2]*b_j1/M_PI + delta_t*delta_t*B_p - pow(delta_t,1.5)*a_j2*xi[j1]/4
								+ pow(delta_t,1.5)*xi[j1]*b_j2/(2*M_PI) + delta_t*delta_t*C_p + delta_t*delta_t*A_p/2;
			III[j1+1][j2+1][0] = delta_t*delta_t*xi[j1]*xi[j2]/2 - pow(delta_t,1.5)*(a_j2*xi[j1] - a_j1*xi[j2])/2 + delta_t*delta_t*A_p
								- III[j1+1][0][j2+1] - III[0][j1+1][j2+1];
		}
	}
											
	// calculate Stratonovich integral J(j1,j2,j3)
	for (int j1 = 0; j1 < m; j1++) {
		
		// coefficients a_j1, b_j1
		double a_j1 = 0, b_j1 = 0;
		for (int r = 1; r <= p; r++) {
			a_j1 += zeta[j1][r-1]/r; 
			b_j1 += eta[j1][r-1]/(r*r);
		}
		a_j1 = - a_j1 * sqrt(2*delta_t)/M_PI - 2*sqrt(delta_t*rho_p)*mu[j1];
		b_j1 = b_j1*sqrt(delta_t/2) + sqrt(delta_t*alpha_p)*phi[j1];
		
		for (int j2 = 0; j2 < m; j2++) {
			
			// calculate A_p
			double A_p = 0;
			for (int r = 1; r <= p; r++) {
				A_p += (zeta[j1][r-1] * eta[j2][r-1] - eta[j1][r-1] * zeta[j2][r-1])/r;
			}
			A_p = A_p/(2*M_PI);
						
			for (int j3 = 0; j3 < m; j3++) {
				if (!((j1 == j2) && (j2 == j3))) {
					//coeffiecient D_p
					double D_p = 0;
					// first calculate the sums that appear in D_p
					double sum1 = 0, sum2 = 0, sum3 = 0;
					
					for (int l = 1; l <= p; l++) {
						for (int r = 1; r <= p; r++) {
							// set eta[j][r] = zeta[j][r] = 0 for r > p
							double zeta_j3_l_r, eta_j1_l_r, eta_j3_l_r;
							if ((l+r) <= p) {
								zeta_j3_l_r = zeta[j3][l+r-1];
								eta_j1_l_r = eta[j1][l+r-1];
								eta_j3_l_r = eta[j3][l+r-1];
							} else {
								zeta_j3_l_r = 0;
								eta_j1_l_r = 0;
								eta_j3_l_r = 0;
							}
															
							sum1 += (zeta[j2][l-1]*(zeta_j3_l_r*eta[j1][r-1] - zeta[j1][r-1]*eta_j1_l_r) 
									+ eta[j2][l-1]*(zeta[j1][r-1]*zeta_j3_l_r + eta[j1][r-1]*eta_j3_l_r)) / (l*(1.*l+r));
						}
					}
					sum1 = - sum1 / (M_PI*M_PI*pow(2,(2.5)));
						
					for (int l = 2; l <= p; l++) {
						for (int r = 1; r <= (l-1); r++) {
							sum2 += (zeta[j2][l-1]*(zeta[j1][r-1]*eta[j3][l-r-1] + zeta[j3][l-r-1]*eta[j1][r-1])
									- eta[j2][l-1]*(zeta[j1][r-1]*zeta[j3][l-r-1] - eta[j1][r-1]*eta[j3][l-r])) / (r*(1.*l-r));
						}
					}
					sum2 = sum2 / (M_PI*M_PI*pow(2,(2.5)));
					
					for (int l = 1; l <= p; l++) {
						for (int r = (l+1); r <= 2*p; r++) {
							
							// set eta[j][r] = zeta[j][r] = 0 for r > p
							double eta_j1_r, zeta_j1_r, eta_j3_r_l, zeta_j3_r_l;
							if (r <= p) {
								eta_j1_r = eta[j1][r-1];
								zeta_j1_r = zeta[j1][r-1];
							} else {
								eta_j1_r = 0;
								zeta_j1_r = 0;
							}
							
							if ((r-l) <= p) {
								eta_j3_r_l = eta[j3][r-l-1];
								zeta_j3_r_l = zeta[j3][r-l-1];
							} else { 
								eta_j3_r_l = 0;
								zeta_j3_r_l = 0;
							}
							
							sum3 += (zeta[j2][l-1]*(zeta_j3_r_l*eta_j1_r - zeta_j1_r*eta_j3_r_l)
									+ eta[j2][l-1]*(zeta_j1_r*zeta_j3_r_l + eta_j1_r*eta_j3_r_l)) / (r*(1.*r-l));
						}
					}
					sum3 = sum3 / (M_PI*M_PI*pow(2,(2.5)));
					
					// calculate D_p
					D_p = sum1 + sum2 + sum3;
					
					// calculate coefficients B_j1_j3, C_j2_j1
					double B_j1_j3 = 0, C_j2_j1 = 0;
					for (int r = 1; r <= p; r++) {
						B_j1_j3 += (zeta[j1][r-1] * zeta[j3][r-1] + eta[j1][r-1] * eta[j3][r-1])/(r*r);
						for (int l = 1; l <= p; l++) {
							if (r != l) {
								C_j2_j1 += (zeta[j2][r-1]*zeta[j1][l-1]/l - l*eta[j2][r-1]*eta[j1][l-1]/r) * (r/(1.*r*r-l*l));
							}
						}
					}
					B_j1_j3 = B_j1_j3/(4*M_PI*M_PI);
					C_j2_j1 = - C_j2_j1/(2*M_PI*M_PI);
					
					
					// now we can calculate the Stratonovic integral JJJ(j1,j2,j3)
					III[j1+1][j2+1][j3+1] = xi[j1]*III[0][j2+1][j3+1]/sqrt(delta_t) + 0.5*a_j1*II[j2+1][j3+1] + delta_t*b_j1*xi[j2]*xi[j3]/(2*M_PI)
											- pow(delta_t,1.5)*xi[j2]*B_j1_j3 + pow(delta_t,1.5)*xi[j3]*(0.5*A_p - C_j2_j1) + pow(delta_t,1.5)*D_p;
				}
			}
		}
	}
	
	// now that the Stratonovich integrals are calculated insert correct values for II(j,j)
	for (int j = 0; j < m; j++) {
		II[j+1][j+1] = 0.5 * (dW[j]*dW[j] - delta_t); 
	}
	
	// calculate triple Ito integrals III(j1,j2,j3)
	for (int j1 = 1; j1 <= m; j1++) {
		for (int j2 = 1; j2 <= m; j2++) {
			for (int j3 = 1; j3 <= m; j3++) {
				if (!((j1 == j2) && (j2 == j3))) { // do this only if j1, j2, j3 are not all equal (because I(j,j,j) is already calculated)
					int Ind_j1_j2 = 0, Ind_j2_j3 = 0; //indicator functions
					if (j1 == j2) {
						Ind_j1_j2 = 1;
					}
					if (j2 == j3) {
						Ind_j2_j3 = 1;
					}
					III[j1][j2][j3] = III[j1][j2][j3]  - 0.5*(Ind_j1_j2 * II[0][j3] + Ind_j2_j3 * II[j1][0]);
				}
			}
		}
	}
			
	
	// MAIN ALGORITHM: calculate Y[0], Y[1], ..., Y[dim-1]
	for (int k = 0; k < dim ; k++) { 
		double sum1, sum2, sum3, sum4; // sums that have to be calculated in the scheme
		sum1 = sum2 = sum3 = sum4 = 0;
		
		// first sum that appears in the scheme
		for (int j = 1; j <= m; j++) {
			sum1 += b(Y,k,j,params) * dW[j-1];
		}
		
		// second and third sum
		for (int j1 = 1; j1 <= m; j1++) {
			
			// calculate sampling points Y_plus and Y_minus
			for (int i = 0; i < dim; i++) {
				Y_plus[i] = Y[i] + (a(Y,i,params)*delta_t)/m + b(Y,i,j1,params)*sqrt(delta_t);
				Y_minus[i] = Y[i] + (a(Y,i,params)*delta_t)/m - b(Y,i,j1,params)*sqrt(delta_t);
			}
				
			for (int j2 = 0; j2 <= m; j2++) {
				if (j2 == 0) { // because b(k,0) = a(k)
					sum2 += (a(Y_plus,k,params) - a(Y_minus,k,params)) * II[j1][j2];
					sum3 += (a(Y_plus,k,params) - 2*a(Y,k,params) + a(Y_minus,k,params)) * II[0][j2];
				} else {
					sum2 += (b(Y_plus,k,j2,params) - b(Y_minus,k,j2,params)) * II[j1][j2];
					sum3 += (b(Y_plus,k,j2,params) - 2*b(Y,k,j2,params) + b(Y_minus,k,j2,params)) * II[0][j2];
				}
			}
		}
		
		// fourth sum
		for (int j1 = 1; j1 <= m; j1++) {
			
			// calculate sampling points Y_plus and Y_minus
			for (int i = 0; i < dim; i++) {
				Y_plus[i] = Y[i] + (a(Y,i,params)*delta_t)/m + b(Y,i,j1,params)*sqrt(delta_t);
				Y_minus[i] = Y[i] + (a(Y,i,params)*delta_t)/m - b(Y,i,j1,params)*sqrt(delta_t);
			}
				
			for (int j2 = 1; j2 <= m; j2++) {
				
				// calculate sampling points Phi_plus and Phi_minus
				for (int i = 0; i < dim; i++) {
					Phi_plus[i] = Y_plus[i] + b(Y_plus,i,j2,params)*sqrt(delta_t);
					Phi_minus[i] = Y_plus[i] - b(Y_plus,i,j2,params)*sqrt(delta_t);
				}
				
				for (int j3 = 1; j3 <= m; j3++) {
					sum4 += (b(Phi_plus,k,j3,params) - b(Phi_minus,k,j3,params) - b(Y_plus,k,j3,params) + b(Y_minus,k,j3,params)) * III[j1][j2][j3];
				}
			}
		}			
		
		sol[k] = Y[k] + a(Y,k,params)*delta_t + sum1 + sum2/(2*sqrt(delta_t)) + sum3/(2*delta_t) + sum4/(2*delta_t);
	}
	
	for (int k = 0; k < dim; k++) {
		Y[k] = sol[k];
	}
	
	// deallocate memory
	free(Y_plus);
	free(Y_minus);
	free(Phi_plus);
	free(Phi_minus);
	free(dW);
	free(sol);
	free(xi);
	for (int i = 0; i < m; i++) {
		free(zeta[i]);
	}
	free(zeta);
	for (int i = 0; i < m; i++) {
		free(eta[i]);
	}
	free(eta);
	free(mu);
	free(phi);
	for (int i = 0; i < m+1; i++) {
		free(II[i]);
	}
	free(II);
	for (int i = 0; i < m+1; i++) {
		for (int j = 0; j < m+1; j++) {
			free(III[i][j]);
		}
		free(III[i]);
	}
	free(III);
	
	return GSL_SUCCESS;
}
