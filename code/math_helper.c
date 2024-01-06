#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

double RadialIntegrate(double (*func)(double), double low, double high, int num_intervals) {
    /*
    Perform integration of a radial function on a disk 
    from r=low to r=high using trapezoidal rule
    */
    double h = (high - low) / num_intervals;
    double result = 0.5 * ((*func)(low) * low + (*func)(high) * high);

    for (int i = 1; i < num_intervals; i++) {
        double x_i = low + i * h;
        result += (*func)(x_i) * x_i;
    }

    result *= h * 2 * M_PI;

    return result;
}

double KernelTanh(double r){
    /*
    Un-normalized tanh kernel with cutoff at 1
    */
    double steep = 5;
    if (r >= 0 && r < 1)
        return (1+tanh(steep*(r+1/2))) * (1+tanh(steep*(1/2-r))) / 4;
    else
        return 0;
}

double KernelExp(double r){
    /*
    Un-normalized exponential kernel (formula in the MIPS review paper)
    */
    if (r >= 0 && r < 1)
        return exp(-1/(1 - r*r));
    else
        return 0;
}

double KernelGauss(double r) {
    return exp(-r*r/2);
}

double Min(int num, ...) {
    /*
    Find the minimum of an arbitrary number of values. 
    The first argument is the number of values to be inputed.
    */
    va_list valist;
    int i;

    va_start(valist, num);

    double min = va_arg(valist, double);
    for (i=0; i<num-1; i++){
        min = fmin(min, va_arg(valist, double));
    }

    va_end(valist);

    return min;
}