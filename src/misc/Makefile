
LIBS    = -lm -fopenmp -lfreeimage -lcfitsio -lhdf5 
# --------------------------------------------

all:   kmean

kmean: kabsDImaxtree2.c
	gcc -o $@ $^  $(LIBS)


clean:
	rm kmean
