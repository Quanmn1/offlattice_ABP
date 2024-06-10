#include "Options.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <setjmp.h>

#ifdef MT
#include "mt19937-64.c" // RNG
#endif

// Implement one huge try-catch block to clean up memory if encounter an error.
jmp_buf ex_buf__; // define the buffer here to be usable in ABP_functions.c. 
#define TRY do{ if( !setjmp(ex_buf__) ){ // If want multiple blocks, reset the buffer in here
#define CATCH } else { // If want multiple kinds of errors, use a switch case statement.
#define ETRY } }while(0) // Terminate after one interation. This is to make the block one complete statement.
#define THROW longjmp(ex_buf__, 1) // Allow for a variable here if want multiple kinds of errors

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

TRY {
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
        fprintf(stderr, "Memory allocation for neighboring_boxes failed.\n");
        ERROR(3);
    }
    // Assign values to neighboring boxes
    ConstructNeighboringBoxes(parameters, neighboring_boxes);

#endif

    double t = 0;
    /*
    Generate initial positions and orientations of particles
    */
    if (parameters.input_file == NULL) {
    #ifdef INIT_SLAB 
        #ifdef PFAP
        SlabLatticeInitialConditions(particles, parameters
        #ifdef HASHING
        , &boxes, &neighbors, neighboring_boxes
        #endif
        );
        #else
        SlabRandomInitialConditions(particles, parameters
        #ifdef HASHING
        , &boxes, &neighbors, neighboring_boxes
        #endif
        );
        #endif
    #else
        #ifdef PFAP
        LatticeInitialConditions(particles, parameters
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
    #endif
    }
    else {
        GivenInitialConditions(parameters.input_file, particles, &parameters, &t
        #ifdef HASHING
        , &boxes, &neighbors, neighboring_boxes
        #endif
        ); // already closed input_file here
    }

    if (t + EPS > parameters.next_store_time) {
        StorePositions(t, parameters, particles
        #ifdef HASHING
        , &boxes, &neighbors
        #endif
        );
        parameters.next_store_time += parameters.store_time_interval;
    }

    double step = parameters.dt;
    double histogram_count = 0;
    double store = 0;

    int progress = 0;
    double next_report_progress = t;
    double duration = parameters.final_time - t;

    // printf("lg\n", t);

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

    fprintf(parameters.param_file, "Starting simulation!\n");
    fflush(parameters.param_file);


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
            next_report_progress += duration/20;
            progress += 5;
        }
    }
    printf("Simulation done!\n");
    THROW; // To transition into cleaning memory.
}
CATCH
{
    /*
    Free up memory: every array.
    If 2d array: look at each column and free each columns
    Try not to declare array inside a function that you call a lot
    */
    free(particles);
    free(command_line_output);
    fflush(stderr);
    fflush(parameters.param_file);
    fclose(parameters.param_file);
    fclose(parameters.data_file);
    fclose(stderr);
#ifdef HASHING
    FreeBoxes(parameters.NxBox, &boxes);
    free(neighbors);
    FreeNeighboringBoxes(&neighboring_boxes, parameters.NxBox, parameters.NyBox);
    #ifdef DENSITY_HISTOGRAM
    fclose(parameters.histogram_file);
    FreeDensity(parameters.number_of_boxes_x, &density_matrix);
    free(density_histogram);
    fclose(parameters.density_file);
    #endif
    printf("Cleanup done!\n");
#endif
}
ETRY;
    return 0; // terminal recognizes 0 for success and others for error (flagged by a red dot).
    // accessible through echo $?
}