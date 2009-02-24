AM_CPPFLAGS = -I${HEADAS}/include -I@top_srcdir@/include
#AM_CFLAGS   = -Wall
AM_CFLAGS   = -Wall -Werror

bin_PROGRAMS = binary_stream

binary_stream_SOURCES = binary_stream.c
binary_stream_LDFLAGS  = @libcfitsio@ -lnsl
binary_stream_LDFLAGS += @libape@ -lreadline
binary_stream_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
binary_stream_LDADD = @top_srcdir@/libraries/libdetectors/libdetectors.a -lm
