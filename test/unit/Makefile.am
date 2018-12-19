AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/progressbar/include

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/build-aux/tap-driver.sh

# Try to do a proper Test setup with cmocka
check_PROGRAMS = unit_test_all random_number_gen
TESTS = unit_test_all random_number_gen

unit_test_all_LDFLAGS = -lcmocka
random_number_gen_LDFLAGS = -lcmocka


random_number_gen_LDADD =@top_builddir@/libsixt/libsixt.la
