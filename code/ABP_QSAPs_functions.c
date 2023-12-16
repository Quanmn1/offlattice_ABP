// Pre-compiler options
#include "./Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


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

    result *= h * M_2_PI;

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

    strcat(*command_line_output, "rho_m ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_min ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_max ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelName ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelWidth ");
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
    parameters[0].rho_m = strtod(argv[i], NULL); i++;
    parameters[0].v_min = strtod(argv[i], NULL); i++;
    parameters[0].v_max = strtod(argv[i], NULL); i++;
    sprintf(input_parameters[0].kernel_name , "%s", argv[i]); i++;
    parameters[0].kernel_width = strtod(argv[i], NULL); i++;
    sprintf(input_parameters[0].name , "%s", argv[i]); i++;
    parameters[0].final_time = strtod(argv[i], NULL); i++;
    parameters[0].next_store_time = strtod(argv[i], NULL); i++;
    parameters[0].store_time_interval = strtod(argv[i], NULL); i++;
    input_parameters[0].Dr = strtod(argv[i], NULL); i++;
#ifdef _MT
    input_parameters[0].seed = (long long) strtod(argv[i], NULL); i++;
#endif
}

void AssignValues(param* parameters, inputparam input_parameters, particle** particles, double** local_densities){
    /*
    This function allocates memory to the particles array, and calculate some parameters
    */

    char filename[1000]; // string of the filename
    // alloc a block of mem to the pointer (array) particles
    particles[0] = (particle*) malloc(sizeof(particle) * parameters[0].N);
    // local_densities[0] = (double*) malloc(sizeof())

    // Assign the files
    sprintf(filename, "%s_param",input_parameters.name); // print name into filename
    parameters[0].param_file = fopen(filename, "w");

    sprintf(filename, "%s_data",input_parameters.name);
    parameters[0].data_file = fopen(filename, "w");

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

double Kernel(double x1, double y1, double x2, double y2, param parameters){
    /*
    Give the value of a normalized kernel at r1 and r2
    */
    double width = parameters.kernel_width;
    double r = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    return (*parameters.kernel)(r/width) / parameters.kernel_normalization / (width*width);
}

double Speed(double rho, double rho_m, double v_min, double v_max){
    /*
    Give the density-dependent speed (formula in the MIPS review paper)
    */
    return v_max + (v_min - v_max)/2 * (1+tanh(2*rho/rho_m-2));
}

double QuorumSensingSpeed(long i, particle* particles, param parameters){
    /*
    This function calculates the density-dependent speed of particle i
    */

    // Calculate rho at the position of particle i
    long j;
    double rho = 0;
    for (j=0;j<parameters.N;j++){
        rho += Kernel(particles[i].x, particles[i].y, 
                      particles[j].x, particles[j].y,
                      parameters);
    }

    // Calculate v(rho)
    double v = Speed(rho, parameters.rho_m, parameters.v_min, parameters.v_max);

    return v;
}

void UpdateParticles(particle* particles, param parameters){   
    /*
    This function updates the positions of particles
    */ 
    long i;
    for (i=0;i<parameters.N;i++){
        // density-dependent speed
        double v = QuorumSensingSpeed(i, particles, parameters); 
        // double v = parameters.v_max; // Non-interacting
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

void StorePositions(double t, param *parameters, particle* particles){
    /*
    This function store the positions and angles of each particles 
    every store_time_interval starting from the inputted next_store_time
    */
    long i;
    if (t >= parameters[0].next_store_time){
        for (i=0;i<parameters[0].N;i++){
            fprintf(parameters[0].data_file,"%lg \t %ld \t %lg \t %lg \t %lg \t",
                    t,i,particles[i].x,particles[i].y,particles[i].theta);
        }
        fprintf(parameters[0].data_file, "\n");
    parameters[0].next_store_time += parameters[0].store_time_interval;
    }
    
    fflush(parameters[0].data_file);
}