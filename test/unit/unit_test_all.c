#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

void test1()
{
  int a = 3;
  int b = 3;
  assert_int_equal(a, b);
}

void test2()
{
  int a = 4;
  int b = 4;
  assert_int_equal(a, b);
}


int main(void)
{

  const struct CMUnitTest tests[] = {
    cmocka_unit_test(test1),
    cmocka_unit_test(test2)
  };

  cmocka_set_message_output(CM_OUTPUT_TAP);

  return cmocka_run_group_tests_name("Default",tests,NULL,NULL);
}
