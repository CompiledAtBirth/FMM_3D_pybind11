CC = g++
LD = g++

FFTW_ROOT = /usr/local/
FFTW_LIB = $(FFTW_ROOT)/lib
FFTW_INCLUDE = $(FFTW_ROOT)/include


CFLAGS  = -c -Wall -O3 -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I ./PBBFMM3D/include/ -I/usr/include -I$(FFTW_INCLUDE) -fopenmp

LDPATH = -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -L/ -L/usr/lib -I/usr/include -I ./PBBFMM3D/include/

LDFLAGS =  -L$(FFTW_LIB) -lfftw3 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group \
-lpthread -lm -ldl -fopenmp

SOURCES =  ./PBBFMM3D/src/kernel_Types.cpp ./PBBFMM3D/src/H2_3D_Tree.cpp

SOURCEA = ./test_PBBFMM3D/PBBFMM_shot_test.cpp

OBJECTA=$(SOURCES:.cpp=.o) $(SOURCEA:.cpp=.o)

EXECUTABLEA= ./exec/test_PBBFMM

test_PBBFMM: $(SOURCES) $(SOURCEA) $(EXECUTABLEA)
$(EXECUTABLEA): $(OBJECTA)
	$(CC)  $(OBJECTA) $(LDPATH) $(LDFLAGS)   -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

