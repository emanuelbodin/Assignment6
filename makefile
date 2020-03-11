CFLAGS= -g -pg -Wall -Ofast

galsim: galsim.o
	gcc -o galsim galsim.o -lm -fopenmp

galsim.o: galsim.c 
	gcc $(CFLAGS) -c galsim.c -fopenmp

clean:
	rm -f ./galsim *.o
