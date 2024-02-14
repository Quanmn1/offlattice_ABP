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
    double** density_matrix; // store the density at each position

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
    AssignValues(&parameters, input_parameters, &particles, &density_histogram, &density_matrix);
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
    #ifdef INIT_SLAB
    SlabInitialConditions(particles, parameters
    #ifdef HASHING
    , &boxes, &neighbors, neighboring_boxes
    #endif
    );
    #else
    RandomInitialConditions(particles, parameters
    #ifdef HASHING
    , &boxes, &neighbors, neighboring_boxes
    #endif
    );
    #endif
    double t = 0;
    StorePositions(t, parameters, particles
    #ifdef HASHING
    , &boxes, &neighbors
    #endif
    );
    parameters.next_store_time += parameters.store_time_interval;

    double step = parameters.dt;
    double histogram_count = 0;
    double store = 0;

    int progress = 0;
    double next_report_progress = 0.0;

    // for (int i = parameters.NyBox-1; i >= 0; i--) {
    //     for (int j=0; j<parameters.NxBox;j++) {
    //         fprintf(parameters.histogram_file, "%ld \t", boxes[j][i]);
    //     }
    //     fprintf(parameters.histogram_file, "\n");
    // }
    // for (int i = 0; i < parameters.N*2; i++) {
    //     fprintf(parameters.histogram_file, "%ld \t", neighbors[i]);
    // }
    // fprintf(parameters.histogram_file, "\n");
    // fflush(parameters.histogram_file);


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
            if (t + EPS > parameters.next_histogram_store) store = 1;
            else store = 0;
            UpdateDensity(density_histogram, density_matrix, parameters, boxes, neighbors, store);
            histogram_count += 1;
            parameters.next_histogram_update += parameters.histogram_update_interval;
        }

        if (t + EPS > parameters.next_histogram_store) {
            StoreDensity(t, parameters, density_histogram, density_matrix, histogram_count, &boxes, &neighbors);
            histogram_count = 0;
            parameters.next_histogram_store += parameters.histogram_store_interval;
        }
        #endif

        if (t > next_report_progress) {
            fprintf(parameters.param_file, "Simulation %d%% done!\n", progress);
            fflush(parameters.param_file);
            next_report_progress += parameters.final_time/10;
            progress += 10;
        }


    }

    /*
    Free up memory: every array.
    If 2d array: look at each column and free each columns
    Try not to declare array inside a function that you call a lot
    */

    free(particles);
    free(command_line_output);
    fclose(parameters.param_file);
    fclose(parameters.data_file);
#ifdef HASHING
    FreeBoxes(parameters.NxBox, &boxes);
    free(neighbors);
    FreeNeighboringBoxes(&neighboring_boxes, parameters.NxBox, parameters.NyBox);
    #ifdef DENSITY_HISTOGRAM
    fclose(parameters.histogram_file);
    FreeDensity(parameters.number_of_boxes_x, &density_matrix);
    fclose(parameters.density_file);
    #endif
#endif

    printf("Simulation done!\n");
    return 0; // terminal recognizes 0 for success and others for error (flagged by a red dot).
    // accessible through echo $?
}