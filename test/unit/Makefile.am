AM_CFLAGS =-I@top_srcdir@/libsixt
AM_CFLAGS+=-I@top_srcdir@/extlib/progressbar/include
AM_CPPFLAGS=-DDATADIR="\"$(srcdir)/data/\""

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/build-aux/tap-driver.sh

# Try to do a proper Test setup with cmocka
check_PROGRAMS = unit_test_all random_number_gen test_genpixgrid test_vignetting
TESTS = unit_test_all random_number_gen test_genpixgrid test_vignetting

unit_test_all_LDFLAGS = -lcmocka
random_number_gen_LDFLAGS = -lcmocka
test_genpixgrid_LDFLAGS = -lcmocka
test_vignetting_LDFLAGS = -lcmocka -lhdio


random_number_gen_LDADD =@top_builddir@/libsixt/libsixt.la
test_genpixgrid_LDADD =@top_builddir@/libsixt/libsixt.la
test_vignetting_LDADD =@top_builddir@/libsixt/libsixt.la

EXTRA_DIST = data
