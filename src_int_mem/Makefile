VLIB= -g -O0

MY_CXX_FLAGS= -std=c++11 #-Wall -Wextra -DNDEBUG
MY_CXX_OPT_FLAGS= -O3 -m64 
MY_CXX=g++

##### Compression mode
##M=0 max
##M=1 mean error
##M=2 constant replacement
##M=3 avg

M=2
MODE = -DM=$(M) 

##### Option for Illumina 8-level binning 
##B=0 disabled
##B=1 enabled
B=0
BIN = -DB=$(B)
LFLAGS = -lm -ldl

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(LFLAGS) $(MODE) $(BIN)

####

LIBOBJ = \
	../external/rankbv/rankbv.o\
	../external/malloc_count/malloc_count.o


all: main 

%.o:%.cpp
	$(MY_CXX) -c -o $@ $<  $(CXX_FLAGS) 

main: bfq_int.cpp ${LIBOBJ} 
	$(MY_CXX) bfq_int.cpp -o bfq_int ${LIBOBJ} $(CXX_FLAGS) 

clean:
	rm -f bfq_int 

valgrind: 
	$(MY_CXX) bfq_int.cpp -o bfq_int $(VLIB)

