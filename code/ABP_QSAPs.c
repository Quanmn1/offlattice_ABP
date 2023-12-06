#include "Options.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _MT
#include "mt19937-64.c"
#endif

typedef struct param {
    double dt;
    long N;
    double v;
    double final_time;
    double next_store_time;
    double store_time_interval;
    double noiseamp;
    FILE* param_file;
    FILE* data_file;
} param;

typedef struct inputparam {
    char name[1000];
    double Dr;
#ifdef _MT
    long long seed;
#endif
} inputparam;

typedef struct particle {
    double x;
    double y;
    double theta;
} particle;

#include "ABP_QSAPs_functions.c"

int main(int argc, char* argv[]) {
    printf("Starting up...");
    char* command_line_output; // pointer to start of string
    param parameters; // struct
    inputparam input_parameters; // struct
    particle* Particles; // pointer to start of array
    ReadInputParameters(argc, argv, command_line_output, &parameters, &input_parameters);
    AssignValues(&parameters, input_parameters, &Particles);
    StoreInputParameters(argc, argv, parameters, input_parameters, command_line_output);

    InitialConditions(Particles, parameters);

    double t = 0;
    while (t < parameters.final_time) {
        UpdateParticles(Particles, parameters);
        t += parameters.dt;
        StorePositions(t, parameters, Particles);
    }

    free(Particles);
    fclose(parameters.data_file);

    return 1;
}