
#debug:
#    cc -Wall -Wextra -g $(SDL2_CFLAGS) $(SDL2_LIBS) -o nanopond nanopond.c -lpthread
gui:
	cc -Ofast -o nanopond nanopond.c
	nvcc -o gpupond gpupond.cu
	@echo "Timing nanopond:"
	@/usr/bin/time -f "Elapsed time: %E" ./nanopond > c1
	@echo "Timing gpupond:"
	@/usr/bin/time -f "Elapsed time: %E" ./gpupond > c2
	diff -u c1 c2
	rm nanopond mod_nanopond c1 c2
clean:
	rm -f *.o nanopond *.dSYM
