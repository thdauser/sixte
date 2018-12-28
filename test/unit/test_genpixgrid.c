#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

#include "sixt.h"

#include "geninst.h"
#include "genpixgrid.h"


const char fname_def_xml[999] = "data/default_inst.xml";
int default_seed=0;

static GenInst* get_mock_geninst(int* const status){

	GenInst* inst = loadGenInst(fname_def_xml, default_seed, status);

	assert_int_equal(inst->det->pixgrid->xwidth,9);
	assert_int_equal(inst->det->pixgrid->ywidth,9);


	return (inst);
}

/**static void get_mock_genpixgrid(int* const status){

	GenPixGrid* grid = newGenPixGrid(status);

	grid->xwidth = 5;
	grid->ywidth = 5;
} **/

static void test_single_point(double x, double y, GenPixGrid* pixgrid){
	double xp= 0.0;
	double yp=0.0;
    int xi = 0;
    int yi = 0;

    getGenDetAffectedPixel(pixgrid,
				    x,y,&xi,&yi, &xp, &yp);

	printf(" %i %i %e %e %e %e \n",xi,yi,x,y,xp,yp);

}

void test_genpixgrid(){
	int status = EXIT_SUCCESS;

	GenInst* inst = get_mock_geninst(&status);


    double delt = inst->det->pixgrid->xdelt;
 //   const double x = 0.0;
//    const double y = 0.0;
    test_single_point(0.0, 0.0, inst->det->pixgrid);
    test_single_point(0.0+2*delt, 0.0+2*delt, inst->det->pixgrid);
    test_single_point(0.0-2*delt, 0.0-2*delt, inst->det->pixgrid);

    // is this grid correct (going from 0 to xwidth-1)???
}


int main(void)
{

  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_genpixgrid)

  };

  cmocka_set_message_output(CM_OUTPUT_TAP);

  return cmocka_run_group_tests_name("Default",tests,NULL,NULL);
}
