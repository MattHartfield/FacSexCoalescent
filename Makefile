GCC_OPTS=-Wall	
FacSexCoalescent: FacSexCoalescent.c
	gcc ${GCC_OPTS} -O3 -o FacSexCoalescent -lm -lgsl -lgslcblas -I/usr/local/include -L/usr/local/lib FacSexCoalescent.c

clean:
	rm -f FacSexCoalescent
