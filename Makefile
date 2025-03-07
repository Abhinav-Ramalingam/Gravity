all: galsim

galsim: galsim.c
	gcc -o galsim galsim.c -lm

galsimo: galsim.c
	gcc -O3 -march=native -ffast-math -o galsim galsim.c -lm

run: galsim.c
	./galsim 3000 ./input_data/ellipse_N_03000.gal 100 1e-5 0
compare:
	./cgf 3000 results.gal ./ref_output_data/ellipse_N_03000_after100steps.gal 
clean: 
	rm -f galsim
