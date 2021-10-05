CC = g++ 
LFLAGS = -lm -ldl

OMP = 1

##### Option for Illumina 8-level binning 
##B=0 disabled
##B=1 enabled

B = 0
BIN = -DB=$(B)

##### Compression mode
##M=0 max
##M=1 mean error
##M=2 constant replacement
##M=3 avg

M = 2
MODE = -DM=$(M) 

DEFINES = $(BIN) $(MODE) -DOMP=$(OMP)

OMP_LIB = -fopenmp

CPPFLAGS = -Wall -ansi -O3 -std=c++11 $(DEFINES) $(OMP_LIB) $(LFLAGS)
#-pedantic -g -O2 

####

LIBOBJ = \
	../external/rankbv/rankbv.o\
	../external/malloc_count/malloc_count.o
####

bfq_ext_obs = bfq_ext.o decode.o ${LIBOBJ}
bfq_ext: $(bfq_ext_obs)
	$(CC) $(OMP_LIB) -o bfq_ext $(bfq_ext_obs) $(LFLAGS)

clean:
	rm -f core *.o *~ bfq_ext

depend:
	$(CC) $(OMP_LIB) -MM *.cpp > dependencies.mk

include dependencies.mk
