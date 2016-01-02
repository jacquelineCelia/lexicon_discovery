# To build a gibbs sampler for the phone-learning project

#CC = g++
CC = /usr/users/chiaying/gcc-4.7/output/bin/g++-4.7
CFLAGS = -c -fopenmp -DQUADPREC -Wall -MMD -O6 -Wall -ffast-math -fno-finite-math-only -finline-functions -fomit-frame-pointer -fstrict-aliasing
#CC = icpc 
SOURCES = py-cfg.cc py-cky.cc interface.cc sampler.cc cluster.cc gmm.cc mixture.cc datum.cc segment.cc bound.cc config.cc toolkit.cc sym.cc rvg.cc gammadist.c mt19937ar.c xtree.cc     
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE = adaptor

INCLIDEFLAGS = -I/usr/users/chiaying/boost_1_53_0/ 

ifeq ($(INTEL_TARGET_ARCH), ia32)
MKL_LINKS=-Wl,--start-group -lmkl_intel -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread
else
MKL_LINKS=-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread
endif

MKL_FLAGS=-I$(MKLROOT)/include -L$(MKLROOT)/lib/$(INTEL_ARCH) $(MKL_LINKS)
IPP_PATHS=-I$(IPPROOT)/include -L$(IPPROOT)/lib/$(INTEL_ARCH)

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(MKL_LINKS) 

.cc.o:
	$(CC) $(CFLAGS) $< -o $@ $(MKL_LINKS) 

clean:
	rm -rf *.o 

# this command tells GNU make to look for dependencies in *.d files
-include $(patsubst %.l,%.d,$(patsubst %.c,%.d,$(SOURCES:%.cc=%.d)))

