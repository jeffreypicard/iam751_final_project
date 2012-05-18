#
# Makefile.am
#

bin_PROGRAMS = fft_jp_test \
fft_jp_mpi_test \
fftw_test \
fftw_2d_test \
fft_jp_2d_test \
mpi_alltoall_transpose_test

lib_LIBRARIES = libfft_jp.a

libfft_jp_a_SOURCES = fft_jp.c \
fft_jp.h

LDADD = libfft_jp.a \
-lfftw3 \
-lm

#TESTS = 

#EXTRA_DIST =
