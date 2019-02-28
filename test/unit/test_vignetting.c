#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

//#include "sixt.h"
#include "vignetting.h"

const char* = "data/dummy_vign.fits";



static int test_vign_load(int seed){
  
  return status;
}


int main(void)
{
  
  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test_vign_load),
  };

  cmocka_set_message_output(CM_OUTPUT_TAP);

  return cmocka_run_group_tests_name("Default",tests,NULL,NULL);
}
