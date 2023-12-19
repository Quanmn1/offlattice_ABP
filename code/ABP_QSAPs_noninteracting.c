#include "Options.h" // Whether to include RNG or not

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _MT
#include "mt19937-64.c" // RNG
#endif

/*
store parameters used in the loop
*/
typedef struct param {
    double dt;
    long N;
    double Lx;
    double Ly;
    double v;
    double (*kernel)(double);
    double kernel_width;
    double kernel_normalization;
    double final_time;
    double next_store_time;
    double store_time_interval;
    double noiseamp;
    FILE* param_file;
    FILE* data_file;
} param;

/*
Store parameters used only at the beginning
*/
typedef struct inputparam {
    char name[1000];
    double Dr;
    char kernel_name[10];
#ifdef _MT
    long long seed;
#endif
} inputparam;

typedef struct particle {
    double x;
    double y;
    double theta;
} particle;

#include "ABP_QSAPs_functions_noninteracting.c"

int main(int argc, char* argv[]) {
    char* command_line_output; // pointer to start of string
    param parameters; // struct
    inputparam input_parameters; // struct
    particle* particles; // pointer to start of array
    double* local_densities;
    /*
    Read command line input into parameters and input_parameters
    */
    ReadInputParameters(argc, argv, &command_line_output, &parameters, &input_parameters);
    /*
    Allocate memory to the particle array; read file name and calculate parameters
    */
    AssignValues(&parameters, input_parameters, &particles, &local_densities);
    /*
    Write the params into a file
    */
    StoreInputParameters(argc, argv, parameters, input_parameters, command_line_output);
    /*
    Generate initial positions and orientations of particles
    */
    InitialConditionsOrigin(particles, parameters);

    double t = 0;
    StorePositions(t, &parameters, particles);
    while (t < parameters.final_time) {
        /*
        Update the positions of the particles
        */
        UpdateParticles(particles, parameters);
        t += parameters.dt;
        /*
        Store positions every store_time_interval
        */
        StorePositions(t, &parameters, particles);
    }

    /*
    Free up memory
    */
    free(particles);
    free(command_line_output);
    fclose(parameters.data_file);

    return 0; // terminal recognizes 0 for success and others for error (flagged by a red dot). accessible through echo #$
}