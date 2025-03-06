all: galsim

galsim: galsim.c
	gcc -o galsim galsim.c -lm
run: galsim.c
	./galsim 3000 ./input_data/ellipse_N_00100.gal 100 1e-5 0
compare:
	./cgf 3000 results.gal ./ref_output_data/ellipse_N_03000_after100steps.gal 
clean: 
	rm -f galsim
