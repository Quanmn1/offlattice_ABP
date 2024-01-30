#include "Options.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef MT
#include "mt19937-64.c" // RNG
#endif

#include "ABP_functions.c"

int main(int argc, char* argv[]) {
    char* command_line_output; // pointer to start of string
    param parameters; // struct
    inputparam input_parameters; // struct
    particle* particles; // pointer to start of array
    double* density_histogram; // store the accumulated density histogram in a time period

#ifdef HASHING
    long** boxes; // first particle in box i, j. malloc in AssignValues
    long* neighbors; // list of neighbors
    box*** neighboring_boxes; // [i][j] is array of neighbors of box i, j. Constant.
#endif

    /*
    Read command line input into parameters and input_parameters
    */
    ReadInputParameters(argc, argv, &command_line_output, &parameters, &input_parameters);
    /*
    Allocate memory to the particle array; read file name and calculate parameters
    */
    AssignValues(&parameters, input_parameters, &particles, &density_histogram);
    /*
    Write the params into a file
    */
    StoreInputParameters(argc, argv, parameters, input_parameters, command_line_output);

#ifdef HASHING
    // Allocate and initialize all empty boxes
    ConstructBoxes(parameters, &boxes); 

    // Allocate and initialize empty neighbors
    ConstructNeighbors(&neighbors, parameters.N);

    // Allocate neighboring boxes
    neighboring_boxes = malloc(parameters.NxBox * sizeof(box**));
    if (neighboring_boxes == NULL) {
        printf("Memory allocation for neighboring_boxes failed.\n");
        exit(3);
    }
    // Assign values to neighboring boxes
    ConstructNeighboringBoxes(parameters, neighboring_boxes);

#endif

    /*
    Generate initial positions and orientations of particles
    */
    InitialConditions(particles, parameters
    #ifdef HASHING
    , &boxes, &neighbors, neighboring_boxes
    #endif
    );
    double t = 0;
    StorePositions(t, parameters, particles
    #ifdef HASHING
    , &boxes, &neighbors
    #endif
    );
    parameters.next_store_time += parameters.store_time_interval;

    double step;
    double histogram_count = 0;

    while (t < parameters.final_time + EPS) {
        /*
        Update the positions of the particles
        */
        UpdateParticles(particles, parameters, &step
        #ifdef HASHING
        , &boxes, &neighbors, neighboring_boxes, t
        #endif
        );
        t += step;

        // printf("Time %lf \n", t);

        /*
        Store positions every store_time_interval
        */
        if (t + EPS > parameters.next_store_time) {
            StorePositions(t, parameters, particles
            #ifdef HASHING
            , &boxes, &neighbors
            #endif
            );
            parameters.next_store_time += parameters.store_time_interval;
        }

        #if defined(HASHING) && defined(DENSITY_HISTOGRAM)
        if (t + EPS > parameters.next_histogram_update) {
            UpdateHistogram(density_histogram, parameters, boxes, neighbors);
            histogram_count += 1;
            parameters.next_histogram_update += parameters.histogram_update_interval;
        }

        if (t + EPS > parameters.next_histogram_store) {
            StoreHistogram(t, parameters, density_histogram, histogram_count, &boxes, &neighbors);
            histogram_count = 0;
            parameters.next_histogram_store += parameters.histogram_store_interval;
        }
        #endif


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
    FreeBoxes(parameters.NxBox, &boxes);
    free(neighbors);
    FreeNeighboringBoxes(&neighboring_boxes, parameters.NxBox, parameters.NyBox);
    #ifdef DENSITY_HISTOGRAM
    fclose(parameters.histogram_file);
    #endif
    #ifdef TESTING
    fclose(parameters.boxes_file);
    #endif
#endif

    return 0; // terminal recognizes 0 for success and others for error (flagged by a red dot).
    // accessible through echo $?
}