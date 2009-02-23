bin_PROGRAMS = binary_stream

binary_stream_SOURCES = binary_stream.c
binary_stream_LDFLAGS  = @libcfitsio@ -lnsl
binary_stream_LDFLAGS += @libape@ -lreadline
binary_stream_LDFLAGS += @libhdinit@ @libhdutils@ @libhdio@
binary_stream_LDADD = -lm
