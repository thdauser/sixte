#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

//#include "sixt.h"
#include "vignetting.h"

#define VIGN_FILENAME "data/dummy_vign.fits"



static Vignetting* vign_load(int* status){
  
	Vignetting* vi  = newVignetting(VIGN_FILENAME, status);

	if (vi==NULL){
		printf(" *** error: trying to load %s, but did not work\n",VIGN_FILENAME);
		*status = EXIT_FAILURE;
	}
	assert_int_equal(*status,EXIT_SUCCESS);

	return vi;
}

void test_vign_load(){

	int status = EXIT_SUCCESS;

	Vignetting* vi = vign_load(&status);

	assert_int_equal(vi!=NULL,1);

	return ;
}

void test_vign_check_dimensions(){

	int status = EXIT_SUCCESS;

	Vignetting* vi = vign_load(&status);

	assert_int_equal(vi->nenergies,2);
	assert_int_equal(vi->ntheta,3);
	assert_int_equal(vi->nphi,1);


	return ;
}


void test_print_values(){

	int status = EXIT_SUCCESS;

	float ref[2][3] = {
			{1.000000e+00, 1.666667e-01, 3.225806e-02},
			{1.000000e+00, 3.846154e-02, 6.622517e-03}
	};

	Vignetting* vi = vign_load(&status);

	for (int ii=0; ii<vi->nenergies; ii++){
		for (int jj=0; jj<vi->ntheta; jj++){
			printf(" E=%f th=%f : vign=%e\n",
					vi->energy[ii],vi->theta[jj]*180/M_PI*60,
					vi->vignet[ii][jj][0]);

			// multiply by 1000 as test is integet
			assert_int_equal(vi->vignet[ii][jj][0]*1e3,ref[ii][jj]*1e3);
		}
	}


	return ;
}

void test_get_vign_factor(){

	int status = EXIT_SUCCESS;

/*	float ref[2][3] = {
			{1.000000e+00, 1.666667e-01, 3.225806e-02},
			{1.000000e+00, 3.846154e-02, 6.622517e-03}
	}; */

	float ref_en = 0.5*(1.666667e-01+3.846154e-02);
	float ref_th = 0.5*(1.666667e-01+3.225806e-02);

	Vignetting* vi = vign_load(&status);
	float phi = 0.0;
	float theta = 5.0/60./180.*M_PI;
	float energy = 3.0;

	printf("fac_en %f \n", get_Vignetting_Factor(vi,  energy, theta, phi));
	float fac_en = get_Vignetting_Factor(vi,  energy, theta, phi);

	assert_int_equal(fac_en*1e4,ref_en*1e4);

	theta = 0.5*(5.0+30.0)/60./180.*M_PI;
	energy = 1.0;

	printf("fac_th %f \n", get_Vignetting_Factor(vi,  energy, theta, phi));
	float fac_th = get_Vignetting_Factor(vi,  energy, theta, phi);
	assert_int_equal(fac_th*1e4,ref_th*1e4);

	return;
}

int main(void)
{
  
  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_vign_load),
	cmocka_unit_test(test_vign_check_dimensions),
	cmocka_unit_test(test_print_values),
	cmocka_unit_test(test_get_vign_factor)
  };

  cmocka_set_message_output(CM_OUTPUT_TAP);

  return cmocka_run_group_tests_name("Default",tests,NULL,NULL);
}
