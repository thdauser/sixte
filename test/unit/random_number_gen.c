#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

//#include "sixt.h"
#include "rndgen.h"

const int default_seed = 0;

unsigned long ref_pseudo_rng[5] = {  // for Seed=0 only!!
		 2357136044,
	     2546248239,
	     3071714933,
	     3626093760,
	     2588848963
};


static int init_rndgen(int seed){
	int status = EXIT_SUCCESS;
	sixt_init_rng(seed, &status);
	assert_int_equal(sixt_use_pseudo_rng(),0); // make sure we are not using the pseudo rng
	return status;
}

static int init_pseudo(){
	int status = EXIT_SUCCESS;
	setenv("SIXTE_USE_PSEUDO_RNG","1",1);
	sixt_init_rng(default_seed, &status);
	unsetenv("SIXTE_USE_PSEUDO_RNG");
	return status;
}

void test_rndgen_init(){
	int status = init_rndgen(default_seed);

	assert_int_equal(status,EXIT_SUCCESS);
	assert_int_equal(sixt_rng_is_initialized(),1);
	assert_int_equal(sixt_use_pseudo_rng(),0); // not the pseudo rng
	sixt_destroy_rng();
}

void test_rndgen_init_pseudo(){
	int status = init_pseudo();

	assert_int_equal(status,EXIT_SUCCESS);
	assert_int_equal(sixt_rng_is_initialized(),1);
	assert_int_equal(sixt_use_pseudo_rng(),1);  // should be the pseudo rng

	sixt_destroy_rng();
}

void test_rndgen_exec(){
	int status = init_rndgen(default_seed);

	assert_int_equal(sixt_use_pseudo_rng(),0); // make sure we are not using the pseudo rng

	const int nrand = 10000;
	double mean = 0.0;

	double val;
	for (int ii=0; ii<nrand; ii++){
		val = sixt_get_random_number(&status);
		assert_in_range(val,0,1);
		mean +=val;
	}
	mean /= nrand;
	assert_in_range(mean,0.4,0.6);  // be very conservative here

	assert_int_equal(status,EXIT_SUCCESS);
	sixt_destroy_rng();
}

void test_random_seed(){

	int status = init_rndgen(default_seed);
	double val1 = sixt_get_random_number(&status);
	sixt_destroy_rng();

	status = init_rndgen(default_seed+1);
	double val2 = sixt_get_random_number(&status);
	sixt_destroy_rng();

	double prec = 1e-12;
	assert_false(fabs(val1-val2) < prec );

}

void test_pseudo_reproducability(){

	const int n = 5;
	unsigned long val[5];

	int status = init_pseudo(default_seed);
	for (int ii=0; ii<n; ii++){
		val[ii] = (unsigned long) (sixt_get_random_number(&status)*4294967296.0); // multiply by 2^32 for int
	}
	sixt_destroy_rng();


	printf(" ** test_pseudo_reproducability**  \n");
	for (int ii=0; ii<n; ii++){
		printf("     %lu == %lu  ?\n",val[ii],ref_pseudo_rng[ii]);
		assert_int_equal(val[ii],ref_pseudo_rng[ii]);
	}

	// just make sure calling it again will produce the same results (should always be given)
	status = init_pseudo(default_seed);
	unsigned long val2 = (unsigned long) (sixt_get_random_number(&status)*4294967296.0); // multiply by 2^32 for int
	sixt_destroy_rng();
	assert_int_equal(val2,ref_pseudo_rng[0]);

}

void test_rndgen_pseudo_exec(){
	int status = init_pseudo();

	const int nrand = 100;
	double mean = 0.0;

	double val;
	for (int ii=0; ii<nrand; ii++){
		val = sixt_get_random_number(&status);
		assert_in_range(val,0,1);
		mean +=val;
//		printf("%i - %e \n",ii,val);
	}
	mean /= nrand;


	assert_int_equal(status,EXIT_SUCCESS);
	sixt_destroy_rng();
}



int main(void)
{

  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_rndgen_init),
    cmocka_unit_test(test_rndgen_init_pseudo),
    cmocka_unit_test(test_rndgen_exec),
    cmocka_unit_test(test_rndgen_pseudo_exec),
    cmocka_unit_test(test_random_seed),
    cmocka_unit_test(test_pseudo_reproducability)

  };

  cmocka_set_message_output(CM_OUTPUT_TAP);

  return cmocka_run_group_tests_name("Default",tests,NULL,NULL);
}
