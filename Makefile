SDL2_CFLAGS = `pkgconf --cflags sdl2`
SDL2_LIBS = `pkgconf --libs sdl2`

#debug:
#    cc -Wall -Wextra -g $(SDL2_CFLAGS) $(SDL2_LIBS) -o nanopond nanopond.c -lpthread
gui:
	cc -Wall -Wextra -Ofast $(SDL2_CFLAGS) $(SDL2_LIBS) -o nanopond nanopond.c -lpthread
	cc -Wall -Wextra -Ofast $(SDL2_CFLAGS) $(SDL2_LIBS) -o mod_nanopond mod_nanopond.c -lpthread
	@echo "Timing nanopond:"
	@/usr/bin/time -f "Elapsed time: %E" ./nanopond > c1
	@echo "Timing mod_nanopond:"
	@/usr/bin/time -f "Elapsed time: %E" ./mod_nanopond > c2
	diff -u c1 c2
	rm nanopond mod_nanopond c1 c2
clean:
	rm -f *.o nanopond *.dSYM
