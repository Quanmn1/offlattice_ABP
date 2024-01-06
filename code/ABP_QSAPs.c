#include "Options.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef MT
#include "mt19937-64.c" // RNG
#endif

#include "ABP_QSAPs_functions.c"

int main(int argc, char* argv[]) {
    char* command_line_output; // pointer to start of string
    param parameters; // struct
    inputparam input_parameters; // struct
    particle* particles; // pointer to start of array

#ifdef HASHING
    int** boxes; // first particle in box i, j. malloc in AssignValues
    int* neighbors; // list of neighbors
    box*** neighboring_boxes; // [i][j] is array of neighbors of box i, j. Constant.
#endif

    /*
    Read command line input into parameters and input_parameters
    */
    ReadInputParameters(argc, argv, &command_line_output, &parameters, &input_parameters);
    /*
    Allocate memory to the particle array; read file name and calculate parameters
    */
    AssignValues(&parameters, input_parameters, &particles);
    /*
    Write the params into a file
    */
    StoreInputParameters(argc, argv, parameters, input_parameters, command_line_output);

#ifdef HASHING
    ConstructBoxes(&parameters, &boxes); // Initialize all empty boxes

    /*
    Allocate neighboring particles
    */
    neighbors = (int*) malloc( (parameters.N+1)*2 * sizeof(int) );
    ConstructNeighbors(&neighbors, parameters.N);

    /*
    Allocate neighboring boxes
    */
    // initialize columns of the matrix: malloc(Nxbox * sizeof(box*)): the matrix is an array of pointers to columns
    // each column is an array of pointers to boxes
    // before malloc: cast into box***
    // for each column: malloc(Nybox * sizeof(box))
    // we do this bc we want to alloc for each of columns
    neighboring_boxes = (box***) malloc(parameters.NxBox * sizeof(box*)); // NxBox: number of indices along x
    #ifdef QSAP
    ConstructNeighboringBoxes(parameters, neighboring_boxes);
    #endif

#endif

    /*
    Generate initial positions and orientations of particles
    */
    InitialConditions(particles, parameters);
    double t = 0;
    StorePositions(t, &parameters, particles);
    while (t < parameters.final_time) {
        /*
        Update the positions of the particles
        */
        UpdateParticles(particles, parameters); // check if box is changed. if so, update boxes
        t += parameters.dt;
        /*
        Store positions every store_time_interval
        */
        StorePositions(t, &parameters, particles);
    }

    /*
    Free up memory: every array.
    If 2d array: look at each column and free each columns
    Try not to declare array inside a function that you call a lot
    */
    
    free(particles);
    free(command_line_output);
    fclose(parameters.data_file);
#ifdef OUTPUT_DENSITY
    fclose(parameters.density_file);
#endif

#ifdef HASHING
    FreeBoxes(&parameters, &boxes);
    FreeNeighbors(&neighbors);
    #ifdef QSAP
        FreeNeighboringBoxes();
    #endif
#endif

    return 0; // terminal recognizes 0 for success and others for error (flagged by a red dot). accessible through echo #$
}