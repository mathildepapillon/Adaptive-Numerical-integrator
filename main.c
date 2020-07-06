//main.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "num_integ.h"

double Kepler(double r); // integrand

int main()
{
    FILE *output;
    double K, mass, pphi, ene, r_up, r_down, dr, r, phi, x, y, tol, limit;
    int nmax, count;

	//paramter values
    K = 1.0/137.0;
    mass = 0.511e+6;
    pphi = 2.0;
    ene = -1.0;

	//calculate bound values
	//upper bound
    r_up= K + pow(K*K-2*fabs(ene)*pphi*pphi/mass,0.5);; 
    r_up /= 2*fabs(ene);

	//lower bound
    r_down= K-pow(K*K-2*fabs(ene)*pphi*pphi/mass,0.5);
    r_down /= 2*fabs(ene);

	//print bounds
    printf("r_down = %lf, r_up = %lf\n",r_down,r_up );

    tol = 1e-10; 
    limit = 10000;

    output = fopen("kepler.dat", "w+"); //for writing
    nmax = 100;
    dr = (r_up-r_down)/nmax;

	//calculate Kepler orbit
    for(int n=0; n<=nmax; n++){
        r = r_down + n*dr;
        count = 0;

        phi = Integrate(Kepler, r_down, r, tol, &count, limit); //integrating to get phi
      
	    //variables from phi
	    x = r*cos(phi); 
        y = r*sin(phi);

        fprintf(output,"%lf %lf\n",x,y);//print to output file
    }

    fclose(output);

    return 0;
}//main

//Kepler integrand
double Kepler(double r)
{
    double f, K, mass, pphi, ene, arg, root;

	//same paramters as in main
    K = 1.0/137.0;
    mass = 0.511e+6;
    pphi = 2.0;
    ene = -1.0;

	// define integrand for total energy eq (124)
    f = pphi/(mass*r*r);
    root= fabs(2.0/mass*(ene-pow(pphi,2)/(2*mass*pow(r,2))+K/r)+1e-16);
    f/= pow(root,0.5);

    return f;

}//Kepler