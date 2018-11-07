AM_CPPFLAGS =-I@top_srcdir@/libsixt
AM_CPPFLAGS+=-I@top_srcdir@/extlib/progressbar/include

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/build-aux/tap-driver.sh

# Try to do a proper Test setup with cmocka
check_PROGRAMS = unit_test_all
TESTS = unit_test_all

unit_test_all_LDFLAGS = -lcmocka