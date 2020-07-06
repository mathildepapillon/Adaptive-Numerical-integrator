#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "num_integ.h"
#define error_handler(tag) { printf(tag); exit(0);}//defined error handler

double test_func(double x); //test integrand
int TestNumInteg();

//to integrate in range h*down to h*up
double NumIntegrate (double *f, double h, int down, int up, int nmax){
    int tmp; 
    double sum;

	//print error if bounds outside of range (0, nmax)
    if(up > nmax || down > nmax) error_handler("bound > nmax ") 
    if(up < 0 || down < 0) error_handler("bound < zero ")

    if(down > up)//integral takes minus sign
	{
        tmp = up;
        up = down;
        down = tmp;
        h *= -1;
    }

    if(up == down)
	{ 
        //no integration domain
		return 0; 
    }
    else if(up == down +1)
	{   
		//single interval
        sum = f[down] + f[down +1];
        sum *= h/2.0;
    }
    else if(up == down +2)
	{   
		// Simpson's Rule
        sum = f[down] + 4*f[down +1] + f[down +2];
        sum *= h/3.0;
    }
    else if(up == down +3)
	{ // Simpson's 3/8 rule
        sum = 3*f[down] + 9*f[down +1] + 9*f[down +2] + 3*f[down +3];
        sum *= h/8.0;
    }
    else if(up == down +4)
	{ 
		// Boole's rule
        sum = 14*f[down] + 64*f[down +1] + 24*f[down +2] + 64*f[down +3] +14*f[down +4];
        sum *= h/45.0;
    }
    else if(up == down +5)
	{   
		// 6 point rule
        sum = 19*f[down] + 75*f[down +1] + 50*f[down +2] + 50*f[down +3] +75*f[down +4]+19*f[down +5];
        sum *= 5*h/45.0;
    }
    else 
	{   
		// 7+ point rule
        sum = (3/8.0)*f[down] + (7/6.0)*f[down +1] + (23/24.0)*f[down +2] + (23/24.0)*f[down +(nmax -2)] + (7/6.0)*f[down +(nmax-1)]+ (3/8.0)*f[down + nmax];
        for(int i = 3; i <= (nmax -3); i++){
            sum += f[down +i];
        }
        sum *= h; 
    }
    return sum;
}//NumIntegrate

//wrapper for NumIntegrate
double NumIntegrateV (double (*func)(double), double x_down, double x_up, int num_intervals)

{
    double f,h;
    double *vf;
    vf = (double *) malloc(sizeof(double)*(num_intervals+1)); //allocate memory
    h = (x_up - x_down)/(num_intervals); // inverval size
    for(int n = 0 ; n <= num_intervals; n++)
	{
		//popultate vf
        vf[n] = func(x_down + n*h);
    }
	//evaluate integral
    f = NumIntegrate(vf, h, 0, num_intervals, num_intervals);

    free(vf);
    return f;
}//NumIntegrateV


double AdaptiveIntegrate(double (*func)(double), double x_down, double h, double *f_prev, double ans_prev, double tol, int *count)
{
    double h_half, ans_l=0, ans_r=0, ans_now, err;
    double fl[3], fr[3]; 
    double w[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0}; //simpson rule weights

	//half interval
    h_half = h/2.0;

	//set up left sid
    fl[0] = f_prev[0];
    fl[2] = f_prev[1];
    fl[1] = func(x_down + h_half);

	//set up right side
    fr[0] = f_prev[1];
    fr[2] = f_prev[2];
    fr[1] = func(x_down + 3*h_half);

    *count += 2;//adavnce count

    for (int i = 0; i < 3; i++)
	{   
		//apply Simpon's rule for each side
        ans_l += (h/2.0)*w[i]*fl[i];
        ans_r += (h/2.0)*w[i]*fr[i];
    }

    ans_now = ans_l + ans_r;//total answer
    err = (ans_now - ans_prev)/15.0; // error
	
	//if error is larger than tolerance then repeat
    if(fabs(err) > tol)
	{
        ans_l = AdaptiveIntegrate(func, x_down, h_half, fl, ans_l, tol, count);
        ans_r = AdaptiveIntegrate(func, x_down+h, h_half, fr, ans_r, tol, count);
        ans_now = ans_l + ans_r;
    }
    else
	{
        ans_now += err; //adjust answer with error
    }

    return ans_now;
}//AdaptiveIntegrate

double Integrate(double (*integrand)(double), double x_down, double x_up, double tol, int *count, int limit) {
    double w[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0}; // weights for Simpsons

    double f_prev[3];
    double ans_prev, f;
    double h = (x_up - x_down)/2;

	//set up f_prev
    f_prev[0]=integrand(x_down); 
    f_prev[1]=integrand(x_down + 2);
    f_prev[2]=integrand(x_down + 4);

	//apply Simpsons rule for each side
    for (int i = 0; i < 3; i++)
	{ 
        ans_prev += (h)*w[i]*f_prev[i];
    }

    *count += 3;//advance count

    f = AdaptiveIntegrate(integrand, x_down, h, f_prev, ans_prev, tol, count);

    if(*count > limit) error_handler("count number larger than limit ") 

    return f;
}//Integrate

//COMMENT OUT EVERYTHING BELOW FOR KEPLER TRIAL


// //to test the integration method
// int main()
// {
//     TestNumInteg();
//     return 1;
// }

// int TestNumInteg()
// {
//     double f_prev[3];
//     double w[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0}; //Simpson weights
//     double (*integrand)(double); // function pointer for the integrand

//     FILE *output;
//     double x_down, x_up, tol, err,f,g;
//     int count, limit, num_intervals;

//     integrand = test_func; //function we want to test

// 	//bounds
//     x_down=0.001;
//     x_up= 10.0;

// 	//true value
//     g = 2.726201989096137;

//     output = fopen("num_integ_test.dat", "w+"); //open file for writing

// 	//trial 1
//     num_intervals=100;
//     f= NumIntegrateV(integrand, x_down, x_up, num_intervals);
//     err=(f-g)/g; //calculate error
//     fprintf(output,"f= %lf, number of intervals = %d, Relative Error =  %lf\n",f, num_intervals, err);
    
// 	//trial 2
// 	num_intervals=1000;
//     f= NumIntegrateV(integrand, x_down, x_up, num_intervals);
//     err=(f-g)/g;
//     fprintf(output,"f= %lf, number of intervals = %d, Relative Error =  %lf\n",f, num_intervals, err);

// 	//trial 3
//     tol= 1.0e-10; //tolerance level
//     limit=10000; //limit on counts
//     count=0;
//     f = Integrate(integrand, x_down, x_up, tol, &count, limit);
//     err=(f-g)/g;
//     fprintf(output,"f= %lf, number of intervals = %d, Relative Error =  %lf\n",f, count, err);

// 	//trial 4
//     num_intervals = count;
//     f = NumIntegrateV(integrand, x_down, x_up, num_intervals);
//     err = (f-g)/g;
//     fprintf(output,"f= %lf, number of intervals = %d, Relative Error =  %lf\n",f, num_intervals, err);

//     fclose(output);

//     return 1;
// }

// double test_func(double x)//test function
// {
//     double f, R, a;
//     f = sin(1.0/x);
//     return f;
// }
