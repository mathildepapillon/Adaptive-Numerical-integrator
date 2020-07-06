A tool for numerical integration. It takes in any function and approximates it to a polynomial, divides the curve into intervals, integrates over each interval, then approximates its own error. It continues to subdivide the intervals until the error respects the given tolerance. I provide an example application with  the  Kepler problem, a special case of 2 Newtonian bodies in motion.

To compile the program, run the following command in Mac Terminal:

	gcc main.c num_integ.c -o exec


This will make a Unix executable named ‘exec’. Execute the source code:

	./exec 

This will output a kepler.dat file. To plot its contents, execute:

	python plots.py 

To only test the integration tool with a test function, uncomment the second half in num_integ.c (it has a main in it) and only compile num_integ.c. Executing it will output a num_integ_test.dat
