PHYS351 Assignment 5
Mathilde Papillon
(Completed on Mac OS)


To compile the program, run the following command in Mac Terminal:

	gcc main.c num_integ.c -o exec


This will make a Unix executable named ‘exec’. Execute the source code:

	./exec 

This will output a kepler.dat file. To plot its contents, execute:

	python plots.py 

To only test the integration tool with a test function, uncomment the second half in num_integ.c (it has a main in it) and only compile num_integ.c. Executing it will output a num_integ_test.dat

Both the data files and Kepler plot can be found in Results sub-directory.