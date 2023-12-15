
#debug:
#    cc -Wall -Wextra -g $(SDL2_CFLAGS) $(SDL2_LIBS) -o nanopond nanopond.c -lpthread
gui:
	cc -Ofast -o nanopond nanopond.c
	cc -Ofast -o nanopond_walt_GPU nanopond_walt_GPU.c
	@echo "Timing nanopond:"
	@/usr/bin/time -f "Elapsed time: %E" ./nanopond > c1
	@echo "Timing gpupond:"
	@/usr/bin/time -f "Elapsed time: %E" ./nanopond_walt_GPU > c2
	diff -u c1 c2
	rm nanopond nanopond_walt_GPU  c1 c2
clean:
	rm -f *.o nanopond *.dSYM
