// Pre-compiler options
#include "./Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "math_helper.c"

void ReadInputParameters(int argc, char* argv[], char** command_line_output, param* parameters, inputparam* input_parameters){
    /*
    This function read the command line and check it has correct number of paramters.
    Then it stores these parameters
    */

    // normally pass by value. by using pointers, pass by reference
    // always one parameter: the program's name
    int number_of_input_parameters = 1;
    int i;
    // command_line_output is a string that contains "usage: " and the correct structure of the command line
    *command_line_output = (char*)malloc(1000);
    if (*command_line_output == NULL) {
        printf("Memory allocation for command_line_output failed.\n");
        exit(3);
    }

    strcpy(*command_line_output, "usage: ");
    strcat(*command_line_output, argv[0]);

    strcat(*command_line_output, " dt ");
    number_of_input_parameters++;

    strcat(*command_line_output, "N ");
    number_of_input_parameters++;

    strcat(*command_line_output, "Lx ");
    number_of_input_parameters++;

    strcat(*command_line_output, "Ly ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelName ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelWidth ");
    number_of_input_parameters++;

    strcat(*command_line_output, "DensityGridSpacing ");
    number_of_input_parameters++;

    strcat(*command_line_output, "FileName ");
    number_of_input_parameters++;

    strcat(*command_line_output, "FinalTime ");
    number_of_input_parameters++;
    
    strcat(*command_line_output, "NextStoreTime ");
    number_of_input_parameters++;

    strcat(*command_line_output, "StoreTimeInterval ");
    number_of_input_parameters++;

    strcat(*command_line_output, "Dr ");
    number_of_input_parameters++;

#ifdef _MT
    strcat(*command_line_output, "seed ");
    number_of_input_parameters++;
#endif

    if(argc != number_of_input_parameters){
        printf("%s\n", *command_line_output);
        exit(1); // Exit the whole program and give a different value compared to 1
    }
    // Now we know the number of params is correct, assign values
    i = 1;
    parameters[0].dt = strtod(argv[i],NULL); i++; // parameters[1] would be next block in memory, meaningless
    parameters[0].N = (long) strtod(argv[i], NULL); i++;
    parameters[0].Lx = strtod(argv[i], NULL); i++;
    parameters[0].Ly = strtod(argv[i], NULL); i++;
    parameters[0].v = strtod(argv[i], NULL); i++;
    sprintf(input_parameters[0].kernel_name , "%s", argv[i]); i++;
    parameters[0].kernel_width = strtod(argv[i], NULL); i++;
    parameters[0].density_grid_spacing = strtod(argv[i], NULL); i++;
    sprintf(input_parameters[0].name , "%s", argv[i]); i++;
    parameters[0].final_time = strtod(argv[i], NULL); i++;
    parameters[0].next_store_time = strtod(argv[i], NULL); i++;
    parameters[0].store_time_interval = strtod(argv[i], NULL); i++;
    input_parameters[0].Dr = strtod(argv[i], NULL); i++;
#ifdef _MT
    input_parameters[0].seed = (long long) strtod(argv[i], NULL); i++;
#endif
}

void AssignValues(param* parameters, inputparam input_parameters, particle** particles, double** local_densities, int store_density){
    /*
    This function allocates memory to the particles array, and calculate some parameters
    */

    char filename[1006]; // string of the filename
    // alloc a block of mem to the pointer (array) particles
    particles[0] = (particle*) malloc(sizeof(particle) * parameters[0].N);
    // local_densities[0] = (double*) malloc(sizeof())

    // Assign the files
    sprintf(filename, "%s_param",input_parameters.name); // print name into filename
    parameters[0].param_file = fopen(filename, "w");

    sprintf(filename, "%s_data",input_parameters.name);
    parameters[0].data_file = fopen(filename, "w");

    // the density grid is centered at (0,0), and the half_number_of_points are the number of points to one side
    parameters[0].half_number_of_points_x = floor(parameters[0].Lx/2 / parameters[0].density_grid_spacing);
    parameters[0].half_number_of_points_y = floor(parameters[0].Ly/2 / parameters[0].density_grid_spacing);

    if (store_density==1){
        sprintf(filename, "%s_density",input_parameters.name);
        parameters[0].density_file = fopen(filename, "w");
        fprintf(parameters[0].density_file, "%d \t %d \t %lg \n", parameters[0].half_number_of_points_x,
                parameters[0].half_number_of_points_y, parameters[0].density_grid_spacing);
    }

    parameters[0].noiseamp = sqrt(2 * parameters[0].dt * input_parameters.Dr);

    // default: 1000 intervals in integration
    long num_intervals = 1000;
    if (strcmp(input_parameters.kernel_name, "exp") == 0) {
        parameters[0].kernel = KernelExp;
        parameters[0].kernel_normalization = RadialIntegrate(parameters[0].kernel, 0, 1, num_intervals);
    }
    else if (strcmp(input_parameters.kernel_name, "tanh") == 0) {
        parameters[0].kernel = KernelTanh;
        parameters[0].kernel_normalization = RadialIntegrate(parameters[0].kernel, 0, 1, num_intervals);
    }
    else {
        printf("Kernel not supported! Supported kernel: tanh, exp");
        exit(2);
    }

#ifdef _MT
    init_genrand64(input_parameters.seed);
#endif
}

void StoreInputParameters(int argc, char* argv[], param parameters, inputparam input_parameters, char* command_line_output){
    /*
    This function writes the parameters to a file of name name_param
    */
    fprintf(parameters.param_file,"%s\n", command_line_output);
    int i;
    for(i=0;i<argc;i++)
        fprintf(parameters.param_file, "%s ", argv[i]);
    
    fprintf(parameters.param_file, "\n");

    fflush(parameters.param_file); // empty the buffer into the file
    
    fclose(parameters.param_file);
}

void InitialConditions(particle* particles, param parameters){
    /*
    This function initializes the positions and angles of each particle 
    to a uniformly random number
    */

    // not modifying particles, just where they point to, so don't need pointers
    long i;
    for (i=0; i<parameters.N; i++){
        particles[i].x = parameters.Lx*(-.5+genrand64_real3()); // From -Lx/2 to Lx/2
        particles[i].y = parameters.Ly*(-.5+genrand64_real3()); // From -Ly/2 to Ly/2
        particles[i].theta = 2*M_PI*genrand64_real3(); //M_2_PI calls 2*Pi in C
    }
}

void InitialConditionsOrigin(particle* particles, param parameters){
    /*
    This function initializes the positions and angles of each particle 
    to a uniformly random number
    */

    // not modifying particles, just where they point to, so don't need pointers
    long i;
    for (i=0; i<parameters.N; i++){
        particles[i].x = 0;
        particles[i].y = 0;
        particles[i].theta = 2*M_PI*genrand64_real3(); //M_2_PI calls 2*Pi in C
    }
}

double Kernel(double x, double y, param parameters){
    /*
    Give the value of a normalized kernel at r1 and r2
    */
    double width = parameters.kernel_width;
    double r = sqrt(x*x + y*y);
    return (*parameters.kernel)(r/width) / parameters.kernel_normalization / (width*width);
}

double Density(double x, double y, particle* particles, param parameters) {
    /*
    Calculate density at point (x,y) using the kernel specified in parameters
    */
    long j;
    double rho = 0;
    for (j=0;j<parameters.N;j++){
        double relative_x = Min(3, fabs(x-particles[j].x),
                        fabs(x-particles[j].x+parameters.Lx),
                        fabs(x-particles[j].x-parameters.Lx));
        double relative_y = Min(3, fabs(y-particles[j].y),
                        fabs(y-particles[j].y+parameters.Ly),
                        fabs(y-particles[j].y-parameters.Ly));
        rho += Kernel(relative_x, relative_y, parameters);
    }

    return rho;
}

double PositionDependentSpeed(double x, double y) {
    if (x>=0 && y>=0) return 0.5;
    else if (x>0 && y<0) return 0.3;
    else if (x<0 && y>0) return 0.2;
    else return 0.1;
}

void UpdateParticles(particle* particles, param parameters){   
    /*
    This function updates the positions of particles
    */ 
    long i;
    for (i=0;i<parameters.N;i++){
        // Non-interacting
        double v = parameters.v;
        // position dependent
        // double v = PositionDependentSpeed(particles[i].x, particles[i].y);
        particles[i].x += parameters.dt * v * cos(particles[i].theta);
        particles[i].y += parameters.dt * v * sin(particles[i].theta);
        particles[i].theta += parameters.noiseamp * gasdev();
        if (particles[i].theta>2*M_PI)
            particles[i].theta -= 2*M_PI;
        if (particles[i].theta<0)
            particles[i].theta += 2*M_PI;

        if (particles[i].x>parameters.Lx/2)
            particles[i].x -= parameters.Lx;
        if (particles[i].x<-parameters.Lx/2)
            particles[i].x += parameters.Lx;

        if (particles[i].y>parameters.Ly/2)
            particles[i].y -= parameters.Ly;
        if (particles[i].y<-parameters.Ly/2)
            particles[i].y += parameters.Ly;

    }
}

void MeanSquare(param parameters, particle* particles, double* x2, double* y2){
    long i;
    for (i=0;i<parameters.N;i++) {
        *x2 += particles[i].x * particles[i].x / parameters.N;
        *y2 += particles[i].y * particles[i].y / parameters.N;
    }
}

void StorePositions(double t, param *parameters, particle* particles, int store_density){
    /*
    This function store the positions and angles of each particles 
    every store_time_interval starting from the inputted next_store_time.
    If the last argument is 1, also store the density at regular grid points.
    */
    long i;
    if (t >= parameters[0].next_store_time){
        double x2 = 0, y2 = 0;
        MeanSquare(parameters[0], particles, &x2, &y2);
        fprintf(parameters[0].data_file,"%lg \t %lg \t", x2, y2);
        for (i=0;i<parameters[0].N;i++){
            fprintf(parameters[0].data_file,"%lg \t %ld \t %lg \t %lg \t %lg \t",
                    t,i,particles[i].x,particles[i].y,particles[i].theta);
        }
        fprintf(parameters[0].data_file, "\n");

        if (store_density==1){
            // Calculate density grid.
            int i, j;
            for (i=-parameters[0].half_number_of_points_x;i<=parameters[0].half_number_of_points_x;i++) {
                for (j=-parameters[0].half_number_of_points_y; j<=parameters[0].half_number_of_points_y;j++) {
                    double rho = Density(i*parameters[0].density_grid_spacing, 
                                        j*parameters[0].density_grid_spacing, particles, parameters[0]);
                    fprintf(parameters[0].density_file, "%lg \t", rho);
                }
                fprintf(parameters[0].density_file, "\n");
            }
            fprintf(parameters[0].density_file, "\n");
            fflush(parameters[0].density_file);
        }

        parameters[0].next_store_time += parameters[0].store_time_interval;
    }
    
    fflush(parameters[0].data_file);
}