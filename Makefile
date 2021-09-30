all: 
	make -C src_ext_mem
	make -C src_int_mem
	make -C external/gsufsort/ TERMINATOR=0 DNA=1 
	make -C external/egap/
	make -C external/libbsc/

clean:
	make clean -C src_ext_mem 
	make clean -C src_int_mem 
	make clean -C external/gsufsort/
	make clean -C external/egap/
	make clean -C external/libbsc/
	rm -f external/rankbv/*.o
	rm -f external/malloc_count/*.o
