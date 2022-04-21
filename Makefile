##### Option for Illumina 8-level binning 
##B=0 disabled
##B=1 enabled
B = 0

##### Compression mode
##M=0 max
##M=1 mean error
##M=2 constant replacement
##M=3 avg
M = 2

all: 
	make -C src_ext_mem B=$(B) M=$(M)
	make -C src_int_mem B=$(B) M=$(M)
	make -C external/gsufsort/ TERMINATOR=0 DNA=1 
	make -C external/egap/
	make -C external/libbsc/
	bash external/install-spring.sh

clean:
	make clean -C src_ext_mem 
	make clean -C src_int_mem 
	make clean -C external/gsufsort/
	make clean -C external/egap/
	make clean -C external/libbsc/
	rm -f external/rankbv/*.o
	rm -f external/malloc_count/*.o
