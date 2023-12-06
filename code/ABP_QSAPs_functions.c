// Pre-compiler options
#include "./Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
This function read the command line and check it has correct number of paramters.
Then it stores these parameters
*/

void ReadInputParameters(int argc, char* argv[], char* command_line_output, param* parameters, inputparam* input_parameters){
    // normally pass by value. by using pointers, pass by reference
    // always one parameter: the program's name
    int number_of_input_parameters = 1;
    int i;
    // command_line_output is a string that contains "usage: " and the correct structure of the command line
    command_line_output = (char*)malloc(1000); // Adjust the size as needed
    if (command_line_output == NULL) {
        printf("Memory allocation for command_line_output failed.\n");
        exit(3);
    }

    command_line_output = "usage: ";
    strcat(command_line_output, argv[0]);

    strcat(command_line_output, "dt ");
    number_of_input_parameters++;

    strcat(command_line_output, "N ");
    number_of_input_parameters++;

    strcat(command_line_output, "v ");
    number_of_input_parameters++;

    strcat(command_line_output, "FileName ");
    number_of_input_parameters++;

    strcat(command_line_output, "FinalTime ");
    number_of_input_parameters++;
    
    strcat(command_line_output, "NextStoreTime ");
    number_of_input_parameters++;

    strcat(command_line_output, "StoreTimeInterval ");
    number_of_input_parameters++;

    strcat(command_line_output, "Dr ");
    number_of_input_parameters++;

#ifdef _MT
    strcat(command_line_output, "seed ");
    number_of_input_parameters++;
#endif

    if(argc != number_of_input_parameters){
        printf("%s\n", command_line_output);
        exit(2); // Exit the whole program and give a different value compared to 1
    }
    // Now we know the number of params is correct, assign values
    i = 1;
    parameters[0].dt = strtod(argv[i],NULL); i++; // parameters[1] would be next block in memory, meaningless
    parameters[0].N = (long) strtod(argv[i], NULL); i++;
    parameters[0].v = strtod(argv[i], NULL); i++;
    sprintf(input_parameters[0].name , "%s", argv[i]); i++;
    parameters[0].final_time = strtod(argv[i], NULL); i++;
    parameters[0].next_store_time = strtod(argv[i], NULL); i++;
    parameters[0].store_time_interval = strtod(argv[i], NULL); i++;
    input_parameters[0].Dr = strtod(argv[i], NULL); i++;
#ifdef _MT
    input_parameters[0].seed = (long long) strtod(argv[i], NULL); i++;
#endif
}

void AssignValues(param* parameters, inputparam input_parameters, particle** Particles){
    char filename[1000]; //string of the filename
    // Assign a length of Particles
    // pointer to the pointer (array) Particles
    Particles[0] = (particle*) malloc(sizeof(particle) * parameters[0].N); // alloc a block of mem to an array of particles

    // Assign the files
    sprintf(filename, "%s_param",input_parameters.name); // print name into filename
    parameters[0].param_file = fopen(filename, "w");

    sprintf(filename, "%s_data",input_parameters.name);
    parameters[0].data_file = fopen(filename, "w");

    parameters[0].noiseamp = sqrt(2 * parameters[0].dt * input_parameters.Dr);
#ifdef _MT
    init_genrand64(input_parameters.seed);
#endif
}

void StoreInputParameters(int argc, char* argv[], param parameters, inputparam input_parameters, char* command_line_output){
    // Write to file
    fprintf(parameters.param_file,"%s\n", command_line_output);
    int i;
    for(i=0;i<argc;i++)
        fprintf(parameters.param_file, "%s ", argv[i]);
    
    fprintf(parameters.param_file, "\n");

    fflush(parameters.param_file); // empty the buffer into the file
    
    fclose(parameters.param_file);
}

void InitialConditions(particle* Particles, param parameters){
    // not modifying Particles, just where they point to, so don't need pointers
    long i;
    for (i=0; i<parameters.N; i++){
        Particles[i].x = 2*(-.5+genrand64_real3());
        Particles[i].y = 2*(-.5+genrand64_real3());
        Particles[i].theta = 2*M_PI*genrand64_real3(); //M_PI calls Pi in C
    }
}

void UpdateParticles(particle* Particles, param parameters){    
    long i;

    for (i=0;i<parameters.N;i++){
        Particles[i].x += parameters.dt * parameters.v * cos(Particles[i].theta);
        Particles[i].y += parameters.dt * parameters.v * sin(Particles[i].theta);
        Particles[i].theta += parameters.noiseamp * gasdev();
        if (Particles[i].theta>2*M_PI)
            Particles[i].theta -= 2*M_PI;
        if (Particles[i].theta<0)
            Particles[i].theta += 2*M_PI;        
    }
}

void StorePositions(double t, param parameters, particle* Particles){
    long i;
    if (t>parameters.next_store_time){
        for (i=0;i<parameters.N;i++){
            fprintf(parameters.data_file,"%lg \t %ld \t %lg \t %lg \t %lg \t", \
                    t,i,Particles[i].x,Particles[i].y,Particles[i].theta);
        }
        fprintf(parameters.data_file, "\n");
    }
    parameters.next_store_time += parameters.store_time_interval;
    
    fflush(parameters.data_file);
}