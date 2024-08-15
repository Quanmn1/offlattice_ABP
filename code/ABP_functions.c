// Pre-compiler options
#include "Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "math_helper.c"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define EPS 1e-7
#define NORMALIZATION_EXP 0.46651239317833015
#define DENSITY_MAX 3

/*
store parameters used in the loop
*/
typedef struct param {
    double dt;
    long N;
    #if defined(INIT_SLAB)
    double rho_small;
    double rho_large;
    double liquid_fraction;
    #endif
    double Lx;
    double Ly;
#ifdef QSAP_TANH
    double rho_m;
    double v_max;
    double v_min;
#elif defined QSAP_EXP || defined QSAP_ZERO
    double rho_m;
    double v;
    double lambda;
    double phi;
#endif
#ifdef PFAP
    double epsilon;
    #if !defined(QSAP)
    double v;
    #endif
    double interaction_range_pfap;
#elif defined NONE
    double v;
#endif
    double interaction_range_qsap;
#ifdef ArbitraryKernel
    double (*kernel)(double);
    double kernel_width;
    double kernel_normalization;
#endif
    double final_time;
    double next_store_time;
    double store_time_interval;
    double Dr;
    double noiseamp;
#ifdef WALL
    double wall_size;
    double omega;
    FILE* wall_pressure_file;
#endif
#ifdef HASHING
    double box_size;
    int NxBox;
    int NyBox;
    #ifdef DENSITY_HISTOGRAM
    int density_box_size; // how many hashing boxes are in a density box
    double density_box_area;
    long max_number;
    int number_of_boxes_x;
    int number_of_boxes_y;
    FILE* histogram_file;
    double next_histogram_update;
    double histogram_update_interval;
    double next_histogram_store;
    double histogram_store_interval;
    #ifdef QSAP_ZERO
    int stopped;
    #endif
    FILE* density_file;
    #endif
    #ifdef STRESS_TENSOR
    FILE** sigmaIK_files;
    FILE* sigmaA_file;
    FILE* sigma_file;
    FILE* nematic_file;
    #endif
#endif
    FILE* param_file;
    FILE* data_file;
    FILE* input_file;
} param;

/*
Store parameters used only at the beginning
*/
typedef struct inputparam {
    char name[1000];
    char resume[10];
#ifdef ArbitraryKernel
    char kernel_name[10];
#endif
#ifdef MT
    long long seed;
#endif
} inputparam;

typedef struct particle {
    double x;
    double y;
    double theta;
    double rho;
    double v;
    double fx;
    double fy;
    double move_x;
    double move_y;
    double speed;
#ifdef HASHING
    int bi;
    int bj;
#endif
#ifdef WALL
    double fw;
#endif
} particle;

#ifdef HASHING
typedef struct box {
    int i;
    int j;
#ifdef PBC
    double epsilon_x;
    double epsilon_y;
#endif
#ifdef STRESS_TENSOR
    double line_x;
    double line_y;
#endif
} box;
#endif

void ERROR(int exit_code) {
    THROW;
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
    *command_line_output = malloc(1000);
    if (*command_line_output == NULL) {
        printf("Memory allocation for command_line_output failed.\n");
        exit(3);
    }

    strcpy(*command_line_output, "usage: ");
    strcat(*command_line_output, argv[0]);

    strcat(*command_line_output, " dt ");
    number_of_input_parameters++;

#ifdef INIT_SLAB
    strcat(*command_line_output, "rho_small ");
    number_of_input_parameters++;

    strcat(*command_line_output, "rho_large ");
    number_of_input_parameters++;

    strcat(*command_line_output, "liquid_fraction ");
    number_of_input_parameters++;
#else
    strcat(*command_line_output, "N ");
    number_of_input_parameters++;
#endif

    strcat(*command_line_output, "Lx ");
    number_of_input_parameters++;

    strcat(*command_line_output, "Ly ");
    number_of_input_parameters++;

#ifdef QSAP_TANH
    strcat(*command_line_output, "rho_m ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_min ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_max ");
    number_of_input_parameters++;

    strcat(*command_line_output, "r_max_qsap ");
    number_of_input_parameters++;

#elif defined QSAP_EXP || defined QSAP_ZERO
    strcat(*command_line_output, "rho_m ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v ");
    number_of_input_parameters++;

    strcat(*command_line_output, "lambda ");
    number_of_input_parameters++;

    strcat(*command_line_output, "phi ");
    number_of_input_parameters++;

    strcat(*command_line_output, "r_max_qsap ");
    number_of_input_parameters++;
#endif

#ifdef PFAP
    #if !defined QSAP
    strcat(*command_line_output, "v ");
    number_of_input_parameters++;
    #endif

    strcat(*command_line_output, "epsilon ");
    number_of_input_parameters++;

    strcat(*command_line_output, "r_max_pfap ");
    number_of_input_parameters++;

#elif defined NONE
    strcat(*command_line_output, "v ");
    number_of_input_parameters++;

#endif

    strcat(*command_line_output, "Dr ");
    number_of_input_parameters++;

    strcat(*command_line_output, "final_time ");
    number_of_input_parameters++;
    
#ifdef ArbitraryKernel
    strcat(*command_line_output, "KernelName ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelWidth ");
    number_of_input_parameters++;
#endif

#ifdef HASHING
    #ifdef WALL
    strcat(*command_line_output, "wall_size ");
    number_of_input_parameters++;

    strcat(*command_line_output, "omega ");
    number_of_input_parameters++;
    #endif

    #ifdef DENSITY_HISTOGRAM
    strcat(*command_line_output, "density_box_size ");
    number_of_input_parameters++;

    strcat(*command_line_output, "next_histogram_update ");
    number_of_input_parameters++;

    strcat(*command_line_output, "histogram_update_interval ");
    number_of_input_parameters++;

    strcat(*command_line_output, "next_histogram_store ");
    number_of_input_parameters++;

    strcat(*command_line_output, "histogram_store_interval ");
    number_of_input_parameters++;
    #endif

#endif   

    strcat(*command_line_output, "next_store_time ");
    number_of_input_parameters++;

    strcat(*command_line_output, "store_time_interval ");
    number_of_input_parameters++;

    strcat(*command_line_output, "file_name ");
    number_of_input_parameters++;

    strcat(*command_line_output, "resume? ");
    number_of_input_parameters++;

#ifdef MT
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
#ifdef INIT_SLAB
    parameters[0].rho_small = strtod(argv[i], NULL); i++;
    parameters[0].rho_large = strtod(argv[i], NULL); i++;
    parameters[0].liquid_fraction = strtod(argv[i], NULL); i++;
#else
    parameters[0].N = (long) strtod(argv[i], NULL); i++;
#endif
    parameters[0].Lx = strtod(argv[i], NULL); i++;
    parameters[0].Ly = strtod(argv[i], NULL); i++;
#ifdef QSAP_TANH
    parameters[0].rho_m = strtod(argv[i], NULL); i++;
    parameters[0].v_min = strtod(argv[i], NULL); i++;
    parameters[0].v_max = strtod(argv[i], NULL); i++;
    parameters[0].interaction_range_qsap = strtod(argv[i], NULL); i++;
#elif defined QSAP_EXP || defined QSAP_ZERO
    parameters[0].rho_m = strtod(argv[i], NULL); i++;
    parameters[0].v = strtod(argv[i], NULL); i++;
    parameters[0].lambda = strtod(argv[i], NULL); i++;
    parameters[0].phi = strtod(argv[i], NULL); i++;
    parameters[0].interaction_range_qsap = strtod(argv[i], NULL); i++;
#endif
#ifdef PFAP
    #if !defined QSAP
    parameters[0].v = strtod(argv[i], NULL); i++;
    #endif
    parameters[0].epsilon = strtod(argv[i], NULL); i++;
    parameters[0].interaction_range_pfap = strtod(argv[i], NULL); i++;
    #if !defined QSAP
    parameters[0].interaction_range_qsap = ceil(parameters[0].interaction_range_pfap); 
    #endif
#elif defined NONE
    parameters[0].v = strtod(argv[i], NULL); i++;
#endif
    parameters[0].Dr = strtod(argv[i], NULL); i++;
    parameters[0].final_time = strtod(argv[i], NULL); i++;

#ifdef ArbitraryKernel
    sprintf(input_parameters[0].kernel_name , "%s", argv[i]); i++;
    parameters[0].kernel_width = strtod(argv[i], NULL); i++;
#endif

#ifdef HASHING
    #ifdef WALL
    parameters[0].wall_size = strtod(argv[i], NULL); i++;
    parameters[0].omega = strtod(argv[i], NULL); i++;
    #endif
    #ifdef DENSITY_HISTOGRAM
    parameters[0].density_box_size = strtod(argv[i], NULL); i++;
    parameters[0].next_histogram_update = strtod(argv[i], NULL); i++;
    parameters[0].histogram_update_interval = strtod(argv[i], NULL); i++;
    parameters[0].next_histogram_store = strtod(argv[i], NULL); i++;
    parameters[0].histogram_store_interval = strtod(argv[i], NULL); i++;
    #endif
#endif
    parameters[0].next_store_time = strtod(argv[i], NULL); i++;
    parameters[0].store_time_interval = strtod(argv[i], NULL); i++;
    sprintf(input_parameters[0].name , "%s", argv[i]); i++;
    sprintf(input_parameters[0].resume , "%s", argv[i]); i++;
#ifdef MT
    input_parameters[0].seed = (long long) strtod(argv[i], NULL); i++;
#endif
}

void AssignValues(param* parameters, inputparam input_parameters, particle** particles, 
                    double** density_histogram, double*** density_matrix, double**** sigmaIK, double*** sigmaA, double*** nematic){
    /*
    This function allocates memory to the particles array, and calculate some parameters
    */
    #ifdef INIT_SLAB
    long N_left = (long) (parameters[0].rho_small * parameters[0].Lx * parameters[0].Ly * (1-parameters[0].liquid_fraction));
    long N_right =  (long) (parameters[0].rho_large * parameters[0].Lx * parameters[0].Ly * parameters[0].liquid_fraction);
    parameters[0].N = N_left + N_right;
    #endif

    char filename[2050]; // string of the filenames

    // Assign the files
    sprintf(filename, "%s_param",input_parameters.name); // print name into filename
    parameters[0].param_file = fopen(filename, "w");
    freopen(filename, "a", stderr);
    #ifdef QSAP_ZERO
    parameters[0].stopped = 0;
    #endif

    int resume = strcmp(input_parameters.resume,"yes");
    char write_mode[] = "w";

    if (resume == 0) {
        sprintf(filename, "%s_video/%s_last_state",input_parameters.name, input_parameters.name);
        parameters[0].input_file = fopen(filename, "r");
        if (parameters[0].input_file == NULL) {
            fprintf(stderr, "Input file not found.\n");
            fflush(stderr);
            ERROR(4);
        }
        write_mode[0] = 'a';
    }
    else {
        parameters[0].input_file = NULL;
    }

    sprintf(filename, "%s_data",input_parameters.name);
    parameters[0].data_file = fopen(filename, write_mode);

    // Calculate the noise amplitude for each time step
    parameters[0].noiseamp = sqrt(2 * parameters[0].dt * parameters[0].Dr);

    // alloc a block of mem to the pointer (array) particles
    particles[0] = malloc(sizeof(particle) * parameters[0].N);
    if (particles[0] == NULL) {
        fprintf(stderr, "Memory allocation for particles failed.\n");
        fflush(stderr);
        ERROR(3);
    }

#ifdef HASHING
    parameters[0].box_size = parameters[0].interaction_range_qsap; // qsap radius always defined, always >= pfap radius
    parameters[0].NxBox = (int) parameters[0].Lx / parameters[0].box_size; // number of spatial hashing boxes
    parameters[0].NyBox = (int) parameters[0].Ly / parameters[0].box_size;

    // Check that Lx, Ly is a multiple of box_size
    if ((fabs(parameters[0].Lx - parameters[0].NxBox*parameters[0].box_size) > EPS) || 
        (fabs(parameters[0].Ly - parameters[0].NyBox*parameters[0].box_size) > EPS)) {
        fprintf(stderr, "Length not an integer multiple of size of box.\n");
        fflush(stderr);
        ERROR(2);
    }

    #ifdef DENSITY_HISTOGRAM
    sprintf(filename, "%s_histogram",input_parameters.name);
    parameters[0].histogram_file = fopen(filename, write_mode);

    sprintf(filename, "%s_density",input_parameters.name);
    parameters[0].density_file = fopen(filename, write_mode);

    // number of density boxes (density_box_size: number of spatial hashing boxes in a density box)
    int number_of_boxes_x = (int) parameters[0].NxBox / parameters[0].density_box_size;
    int number_of_boxes_y = (int) parameters[0].NyBox / parameters[0].density_box_size;
    parameters[0].number_of_boxes_x = number_of_boxes_x; 
    parameters[0].number_of_boxes_y = number_of_boxes_y;
    if ((abs(parameters[0].NxBox - number_of_boxes_x*parameters[0].density_box_size) > EPS) || 
        (abs(parameters[0].NyBox - number_of_boxes_y*parameters[0].density_box_size) > EPS)) {
        fprintf(stderr, "System does not contain integer numbers of density boxes: %d %d.\n", number_of_boxes_x,number_of_boxes_y);
        fflush(stderr);
        ERROR(2);
    }

    parameters[0].density_box_area = parameters[0].density_box_size*parameters[0].density_box_size * \
                                     parameters[0].box_size*parameters[0].box_size;
    #ifdef PFAP
    // Almost total exclusion, so max density when close packing (and prolly achieved regardless of starting rho).  
    #ifdef QSAP
    if (parameters[0].interaction_range_pfap < EPS)
        parameters[0].max_number = parameters[0].rho_m * parameters[0].density_box_area * 10;
    else {
        // PFAP+QSAP
        double density_max = DENSITY_MAX / (parameters[0].interaction_range_pfap*parameters[0].interaction_range_pfap);
        parameters[0].max_number = ceil(density_max*parameters[0].density_box_area);
    }
    #else
    double density_max = DENSITY_MAX / (parameters[0].interaction_range_pfap*parameters[0].interaction_range_pfap);
    parameters[0].max_number = ceil(density_max*parameters[0].density_box_area);
    #endif
    #else
    parameters[0].max_number = parameters[0].rho_m * parameters[0].density_box_area * 10;
    #endif

    density_histogram[0] = malloc((parameters[0].max_number+1) * sizeof(double));
    if (density_histogram[0] == NULL) {
        fprintf(stderr, "Memory allocation for density_histogram failed.\n");
        fflush(stderr);
        ERROR(3);
    }

    density_matrix[0] = malloc(number_of_boxes_x * sizeof(double*));
    if (density_matrix[0] == NULL) {
        fprintf(stderr, "Memory allocation for density_matrix failed.\n");
        fflush(stderr);
        ERROR(3);
    }
    for (int i = 0; i < number_of_boxes_x; i++) {
        density_matrix[0][i] = malloc(number_of_boxes_y * sizeof(double));
    }
    #endif

    #ifdef WALL
    // wall pressure
    sprintf(filename, "%s_wall_pressure", input_parameters.name);
    parameters[0].wall_pressure_file = fopen(filename, write_mode);
    #endif

    #ifdef STRESS_TENSOR
    sprintf(filename, "%s_sigma", input_parameters.name);
    parameters[0].sigma_file = fopen(filename, write_mode);

    // IK pressure
    struct stat st = {0};
    char directory[1010];
    sprintf(directory, "%s_sigmaIK", input_parameters.name);
    if (stat(directory, &st) == -1) {
        mkdir(directory, 0777);
    }
    parameters[0].sigmaIK_files = malloc(4 * sizeof(FILE*));
    sprintf(filename, "%s/%s_sigmaIK_xx", directory, input_parameters.name);
    parameters[0].sigmaIK_files[0] = fopen(filename, write_mode);

    sprintf(filename, "%s/%s_sigmaIK_xy", directory, input_parameters.name);
    parameters[0].sigmaIK_files[1] = fopen(filename, write_mode);

    sprintf(filename, "%s/%s_sigmaIK_yx", directory, input_parameters.name);
    parameters[0].sigmaIK_files[2] = fopen(filename, write_mode);

    sprintf(filename, "%s/%s_sigmaIK_yy", directory, input_parameters.name);
    parameters[0].sigmaIK_files[3] = fopen(filename, write_mode);

    sigmaIK[0] = malloc(4*sizeof(double**));
    if (sigmaIK[0] == NULL) {
        fprintf(stderr, "Memory allocation for sigmaIK failed.\n");
        fflush(stderr);
        ERROR(3);
    }
    for (int i=0;i<4;i++) {
        sigmaIK[0][i] = malloc(parameters[0].NxBox * sizeof(double*));
        if (sigmaIK[0][i] == NULL) {
            fprintf(stderr, "Memory allocation for sigmaIK failed.\n");
            fflush(stderr);
            ERROR(3);
        }
        for (int j=0;j < parameters[0].NxBox; j++) {
            sigmaIK[0][i][j] = malloc(parameters[0].NyBox * sizeof(double));
            if (sigmaIK[0][i][j] == NULL) {
                fprintf(stderr, "Memory allocation for sigmaIK failed.\n");
                fflush(stderr);
                ERROR(3);
            }
        }
    }

    // active pressure
    sprintf(filename, "%s_sigmaAxx", input_parameters.name);
    parameters[0].sigmaA_file = fopen(filename, write_mode);

    sigmaA[0] = malloc(parameters[0].NxBox * sizeof(double*));
    if (sigmaA[0] == NULL) {
        fprintf(stderr, "Memory allocation for sigmaA failed.\n");
        fflush(stderr);
        ERROR(3);
    }
    for (int i=0;i < parameters[0].NxBox; i++) {
        sigmaA[0][i] = malloc(parameters[0].NyBox * sizeof(double));
        if (sigmaA[0][i] == NULL) {
            fprintf(stderr, "Memory allocation for sigmaA failed.\n");
            fflush(stderr);
            ERROR(3);
        }
    }

    // nematic tensor Qxx
    sprintf(filename, "%s_Qxx", input_parameters.name);
    parameters[0].nematic_file = fopen(filename, write_mode);

    nematic[0] = malloc(parameters[0].NxBox * sizeof(double*));
    if (sigmaA[0] == NULL) {
        fprintf(stderr, "Memory allocation for Qxx failed.\n");
        fflush(stderr);
        ERROR(3);
    }
    for (int i=0;i < parameters[0].NxBox; i++) {
        nematic[0][i] = malloc(parameters[0].NyBox * sizeof(double));
        if (nematic[0][i] == NULL) {
            fprintf(stderr, "Memory allocation for Qxx failed.\n");
            fflush(stderr);
            ERROR(3);
        }
    }
    #endif
#endif

#ifdef ArbitraryKernel
    // default: 1000 intervals in integration
    int num_intervals = 1000;
    if (strcmp(input_parameters.kernel_name, "exp") == 0) {
        parameters[0].kernel = KernelExp;
        parameters[0].kernel_normalization = RadialIntegrate(parameters[0].kernel, 0, 1, num_intervals);
    }
    else if (strcmp(input_parameters.kernel_name, "tanh") == 0) {
        parameters[0].kernel = KernelTanh;
        parameters[0].kernel_normalization = RadialIntegrate(parameters[0].kernel, 0, 1, num_intervals);
    }
    else {
        fprintf(stderr, , "Kernel not supported! Supported kernel: tanh, exp.\n");
        fflush(stderr);
        ERROR(2);
    }
#endif
#ifdef MT
    init_genrand64(input_parameters.seed);
#endif
}

void StoreInputParameters(int argc, char* argv[], param parameters, inputparam input_parameters, char* command_line_output){
    /*
    This function writes the parameters to a file of name name_param
    */
    fprintf(parameters.param_file,"%s\n", command_line_output);
    for(int i=0;i<argc;i++)
        fprintf(parameters.param_file, "%s ", argv[i]);
    
    fprintf(parameters.param_file, "\n");

    fprintf(parameters.param_file, "%ld\n", parameters.N);

    double rho = parameters.N/parameters.Lx/parameters.Ly;
    #ifdef PFAP
    double pe = parameters.v/parameters.Dr/(parameters.interaction_range_pfap*0.89);
    fprintf(parameters.param_file, "Density = %lf; Pe = %lf \n", rho, pe);
    #elif defined QSAP_TANH
    double v_ratio = parameters.v_max/parameters.v_min;
    fprintf(parameters.param_file, "Rho = %lf \t v0/v1 = %lf \n", rho, v_ratio);
    #elif defined QSAP_EXP || defined QSAP_ZERO
    fprintf(parameters.param_file, "Rho = %lf \t lambda = %lf \n", rho, parameters.lambda);
    #endif

    fflush(parameters.param_file); // empty the buffer into the file
}

double Kernel(double x, double y, double rmax_squared){
    /*
    Give the value of a normalized kernel at (x,y). Use exponential kernel.
    */
    double r_squared = x*x + y*y;
    if (r_squared < rmax_squared) {
        double Z = NORMALIZATION_EXP * rmax_squared;
        return 1/Z * exp(-rmax_squared/(rmax_squared - r_squared));
    }
    else
        return 0;
}

#ifdef ArbitraryKernel
double KernelArbitrary(double x, double y, param parameters){
    /*
    Give the value of a normalized kernel at r1 and r2.
    */
    double width = parameters.kernel_width;
    double r = sqrt(x*x + y*y);
    return (*parameters.kernel)(r/width) / parameters.kernel_normalization / (width*width); //Z=0.4665... * r0^2
}
#endif

double DensityDependentSpeed(double rho, double rho_m, double v_min, double v_max){
    /*
    Give the density-dependent speed (formula in the MIPS review paper)
    */
    return v_max + (v_min - v_max)/2 * (1+tanh(2*rho/rho_m-2));
}

double DensityDependentSpeed2(double rho, double rho_m, double v, double lambda, double phi){
    /*
    Give the density-dependent speed (formula used by Gianmarco)
    */
    return v*exp(-lambda*tanh((rho-rho_m)/phi));
}

double DensityDependentSpeed3(double rho, double rho_m, double v, double lambda, double phi) {
    /*
    Smooth.
    Give speed that is v at rho=0, reaches v(1-lambda) at finite rho_m
    */
    return v * (1 - lambda * (1 - SmoothZero(rho, rho_m, phi) ));
}

double DensityDependentSpeed4(double rho, double rho_m, double v, double lambda, double phi) {
    /*
    Linear.
    Give speed that is v at rho<rho_m-phi/2, v(1-lambda) at rho>rho_m+phi/2, linearly interpolate in between
    */    
    if (rho < rho_m - phi) return v;
    else if (rho < rho_m) return v * (1 - lambda + lambda * (rho_m-rho)/phi);
    else return v * (1 - lambda);
}

double PositionDependentSpeed(double x, double y){
    /*
    Give a quadrant-dependent speed 
    */
    if (x>=0 && y>=0) return 0.5;
    else if (x>0 && y<0) return 0.3;
    else if (x<0 && y>0) return 0.2;
    else return 0.1;
}

double ForceHarmonic(double r, double rmax, double epsilon) {
    /*
    Return the magnitude of the harmonic force at distance r
    */
    double ratio = r / rmax;
    if (ratio <= 1)
        return 2*epsilon/rmax * (1-ratio);
    else
        return 0;
}

double ForceWCA(double r, double rmax, double epsilon) {
    /*
    Return the magnitude of the WCA force at distance r
    */
    double ratio = rmax / r;
    if (ratio >= 1)
        return 12*epsilon/rmax * (pow(ratio, 13) - pow(ratio, 7));
    else
        return 0;
}

double Force(double r, double rmax, double epsilon) {
    #ifdef HARMONIC
    return ForceHarmonic(r, rmax, epsilon);
    #elif defined WCA
    return ForceWCA(r, rmax, epsilon);
    #else
    return 0;
    #endif
}

double ForceWall(double x, double wall_size, double L, double omega) {
    /*
    Return the (signed scalar) force exerted on particles by walls placed at 0 and L, with width wall_size.
    Force formula: F = - nu * omega * (x-x_wall)^(nu-1), nu=2
    Typical values:
    */
    if (x < wall_size) return 2 * omega * (wall_size - x);
    else if (x > L - wall_size) return 2 * omega * (L - wall_size - x);
    else return 0;
}

double BruteForceDensity(double x, double y, particle* particles, param parameters) {
    /*
    Calculate density at point (x,y) using the kernel specified in parameters
    */
    int j;
    double dx, dy;
    double rmax_squared = parameters.interaction_range_qsap*parameters.interaction_range_qsap;
    double rho = 0;
    for (j=0;j<parameters.N;j++){
        dx = Min(3, fabs(x-particles[j].x),
                        fabs(x-particles[j].x+parameters.Lx),
                        fabs(x-particles[j].x-parameters.Lx));
        dy = Min(3, fabs(y-particles[j].y),
                        fabs(y-particles[j].y+parameters.Ly),
                        fabs(y-particles[j].y-parameters.Ly));
        rho += Kernel(dx, dy, rmax_squared);
    }

    return rho;
}

// #ifdef PFAP
// void BruteForceForce(double* fx, double* fy, long i, particle* particles, param parameters) {
//     /*
//     Calculate force at point (x,y) using the kernel specified in parameters
//     */
//     int j;
//     double dx, dy, dr;
//     double rmax = parameters.interaction_range_pfap;
//     double epsilon = parameters.epsilon;
//     double force_magnitude;
//     double x = particles[i].x;
//     double y = particles[i].y;
//     fx[0] = 0;
//     fy[0] = 0;
//     for (j=0;j<parameters.N;j++){
//         if (j != i) {
//             dx = AbsMin(3, x-particles[j].x,
//                             x-particles[j].x+parameters.Lx,
//                             x-particles[j].x-parameters.Lx);
//             dy = AbsMin(3, y-particles[j].y,
//                             y-particles[j].y+parameters.Ly,
//                             y-particles[j].y-parameters.Ly);
//             dr = sqrt(dx*dx + dy*dy);
//             force_magnitude = Force(dr, rmax, epsilon);
//             fx[0] += force_magnitude * dx/dr;
//             fy[0] += force_magnitude * dy/dr;
//         }
//     }
// }
// #endif

#ifdef HASHING
void GetBox(int* bi, int* bj, double x, double y, double Lx, double Ly, double box_size) {
    // Get the box indices from coordinates
    bi[0] = floor( x / box_size);
    bj[0] = floor( y / box_size);
}

void ConstructBoxes(param parameters, long*** boxes) {
    /*
    Construct the array containing the first particle of each box
    */
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    int i, j;
    // Construct and initialize boxes
    boxes[0] = malloc(NxBox * sizeof(long*));
    if (boxes[0] == NULL) {
        fprintf(stderr, "Memory allocation for boxes failed.\n");
        fflush(stderr);
        ERROR(3);
    }
    for ( i = 0; i < NxBox; i++) {
        boxes[0][i] = malloc(NyBox * sizeof(long));
        for ( j = 0; j < NyBox; j++) {
            boxes[0][i][j] = -1;
        }
    }
}

void ConstructNeighbors(long** neighbors, long N) {
    // set all entries to -1: no relation yet between particles
	neighbors[0] = malloc((2*N) * sizeof(long));
    if (neighbors[0] == NULL) {
        fprintf(stderr, "Memory allocation for neighbors failed.\n");
        fflush(stderr);
        ERROR(3);
    }
	for( int i = 0 ; i < 2*N; i++)
		neighbors[0][i] = -1;
}

void FreeBoxes(long*** boxes, int NxBox) {
    for (int i = 0; i < NxBox; i++)
        free(boxes[0][i]);
    free(boxes[0]);
}

void FreeDensity(double*** density_matrix, int number_of_boxes_x) {
    for (int i = 0; i < number_of_boxes_x; i++)
        free(density_matrix[0][i]);
    free(density_matrix[0]);
}

void FreeNeighboringBoxes(box**** neighboring_boxes, int NxBox, int NyBox) {
    int bi, bj;
    for ( bi=0; bi<NxBox; bi++) {
        for ( bj=0; bj<NyBox; bj++) {
            free(neighboring_boxes[0][bi][bj]);
        }
        free(neighboring_boxes[0][bi]);
    }
    free(neighboring_boxes[0]);
}

void FreeSigma(double**** sigmaIK, int NxBox) {
    int bi, bj;
    for ( bi=0; bi<4; bi++) {
        for ( bj=0; bj<NxBox; bj++) {
            free(sigmaIK[0][bi][bj]);
        }
        free(sigmaIK[0][bi]);
    }
    free(sigmaIK[0]);
}

void AddInBox(long index_particle, int bi, int bj, long*** boxes, long** neighbors, particle* particles) {
    long k = boxes[0][bi][bj]; // Save the previous first particle of the box
    boxes[0][bi][bj] = index_particle; // particle index_particle becomes the new first of box
    neighbors[0][2*index_particle] = -1; // index_part doesn't have preceding particle
    neighbors[0][2*index_particle+1] = k; // k becomes succeeding of index_particle
    if (k != -1) neighbors[0][2*k] = index_particle; // if k is a particle, set its preceding to index_particle
    particles[index_particle].bi = bi;
    particles[index_particle].bj = bj;
}

void RemoveFromBox(long index_particle, int bi, int bj, long*** boxes, long** neighbors) {
    long next = neighbors[0][2*index_particle+1]; // Store the particle after index_particle
    long prev;
    // check if index_particle is the first one of the box
    if (boxes[0][bi][bj] == index_particle) {
        // The first particle of the box becomes the one after idx_part
        boxes[0][bi][bj] = next;
        if (next != -1) neighbors[0][2*next] = -1;
    }
    else {
        // get the previous one
        prev = neighbors[0][2*index_particle];
        if (prev == -1) {
            fprintf(stderr, "boxes and neighbors uncompatible: %d %d for particle %ld\n", bi, bj, index_particle);
            fflush(stderr);
            ERROR(5);
        }
        // The next of the previous becomes the current next
        neighbors[0][2*prev+1] = next;
        // If next is a particle, its prev becomes the current prev
        if (next != -1) neighbors[0][2*next] = prev;
    }
}

void ConstructNeighboringBoxes(param parameters, box*** neighboring_boxes) {
    /*
    Neighboring box 0: self
    Neighboring box 1-4: up and right of self
    epsilons: calculate distance through the periodic boundary
    */
    int bi, bj;
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    double Lx = parameters.Lx;
    double Ly = parameters.Ly;

    for (bi=0; bi<NxBox; bi++) {
        neighboring_boxes[bi] = malloc(NyBox*sizeof(box*));
        for (bj=0; bj<NyBox; bj++) {
            neighboring_boxes[bi][bj] = malloc(5*sizeof(box));

            neighboring_boxes[bi][bj][0].i = bi;
            neighboring_boxes[bi][bj][0].j = bj;
            neighboring_boxes[bi][bj][0].epsilon_x = 0;
            neighboring_boxes[bi][bj][0].epsilon_y = 0;
            #ifdef STRESS_TENSOR
            neighboring_boxes[bi][bj][0].line_x = (bi+1)*parameters.box_size; // doesnt matter
            neighboring_boxes[bi][bj][0].line_y = (bj+1)*parameters.box_size; // doesnt matter
            #endif

            neighboring_boxes[bi][bj][1].i = bi;
            neighboring_boxes[bi][bj][1].j = (bj+1)%NyBox;
            neighboring_boxes[bi][bj][1].epsilon_x = 0;
            neighboring_boxes[bi][bj][1].epsilon_y = (bj==NyBox-1) ? Ly : 0;
            #ifdef STRESS_TENSOR
            neighboring_boxes[bi][bj][1].line_x = (bi+1)*parameters.box_size; // doesnt matter
            neighboring_boxes[bi][bj][1].line_y = (bj+1)*parameters.box_size; 
            #endif

            neighboring_boxes[bi][bj][2].i = (bi+1)%NxBox;
            neighboring_boxes[bi][bj][2].j = (bj+1)%NyBox;
            neighboring_boxes[bi][bj][2].epsilon_x = (bi==NxBox-1) ? Lx : 0;
            neighboring_boxes[bi][bj][2].epsilon_y = (bj==NyBox-1) ? Ly : 0;
            #ifdef STRESS_TENSOR
            neighboring_boxes[bi][bj][2].line_x = (bi+1)*parameters.box_size;
            neighboring_boxes[bi][bj][2].line_y = (bj+1)*parameters.box_size;
            #endif

            neighboring_boxes[bi][bj][3].i = (bi+1)%NxBox;
            neighboring_boxes[bi][bj][3].j = bj;
            neighboring_boxes[bi][bj][3].epsilon_x = (bi==NxBox-1) ? Lx : 0;
            neighboring_boxes[bi][bj][3].epsilon_y = 0;
            #ifdef STRESS_TENSOR
            neighboring_boxes[bi][bj][3].line_x = (bi+1)*parameters.box_size; 
            neighboring_boxes[bi][bj][3].line_y = (bj+1)*parameters.box_size; // doesnt matter
            #endif

            neighboring_boxes[bi][bj][4].i = (bi+1)%NxBox;
            neighboring_boxes[bi][bj][4].j = (bj+NyBox-1)%NyBox;
            neighboring_boxes[bi][bj][4].epsilon_x = (bi==NxBox-1) ? Lx : 0;
            neighboring_boxes[bi][bj][4].epsilon_y = (bj==0) ? -Ly : 0;
            #ifdef STRESS_TENSOR
            neighboring_boxes[bi][bj][4].line_x = (bi+1)*parameters.box_size;
            neighboring_boxes[bi][bj][4].line_y = bj*parameters.box_size;
            #endif
        }
    }
}

void MeasureDensityAndForce(long* neighbors, particle* particles, long** boxes, param* parameters, box*** neighboring_boxes) {
    /*
    Measure the density and forces around each particle by looking at particles in adjacent boxes 
    and compute their contributions.
    */
    int bi, bj, nbi, nbj;
    long i, j, k;
    double epsx=0, epsy=0;
    double rmax_squared = parameters[0].interaction_range_qsap*parameters[0].interaction_range_qsap;
    double Z = NORMALIZATION_EXP * rmax_squared;
    #ifdef PFAP
    double rmax_pfap = parameters[0].interaction_range_pfap;
    double epsilon = parameters[0].epsilon;
    double force_magnitude;
    #endif
    #ifdef WALL
    double wall_size = parameters[0].wall_size;
    double Lx = parameters[0].Lx;
    double omega = parameters[0].omega;
    #endif
    double self_density = 1/Z/M_E;
    int NxBox = parameters[0].NxBox;
    int NyBox = parameters[0].NyBox;
    double xi, yi, xj, yj, dx, dy, dr;
    double density_ij;
    // Initialize density fieds to zero for all particles
    for (i=0; i<parameters[0].N; i++) {
        particles[i].rho = self_density; // include the self density
        #ifdef WALL
        particles[i].fw = ForceWall(particles[i].x, wall_size, Lx, omega);
        particles[i].fx = particles[i].fw;
        #else
        particles[i].fx = 0;
        #endif
        particles[i].fy = 0;
    }
    // Loop through all boxes
    for (bi = 0; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            // First particle of box
            i = boxes[bi][bj];
            while (i!=-1) {
                xi = particles[i].x;
                yi = particles[i].y;

                // particles in same box
                j = neighbors[2*i+1]; // consider pairs (i,j) with j after i
                while (j != -1) {
                    xj = particles[j].x;
                    yj = particles[j].y;
                    dx = xi - xj;
                    dy = yi - yj;
                    density_ij = Kernel(dx, dy, rmax_squared);
                    particles[i].rho += density_ij;
                    particles[j].rho += density_ij;

                    #ifdef PFAP
                    dr = sqrt(dx*dx + dy*dy);
                    force_magnitude = Force(dr, rmax_pfap, epsilon);
                    particles[i].fx += force_magnitude * dx/dr;
                    particles[i].fy += force_magnitude * dy/dr;
                    particles[j].fx -= force_magnitude * dx/dr;
                    particles[j].fy -= force_magnitude * dy/dr;
                    #endif
                    // go to next particle
                    j = neighbors[2*j+1];
                }

                // particles in neighboring boxes
                for (k=1; k<5; k++) {
                    nbi = neighboring_boxes[bi][bj][k].i;
                    nbj = neighboring_boxes[bi][bj][k].j;
                    #ifdef PBC
                    epsx = neighboring_boxes[bi][bj][k].epsilon_x;
                    epsy = neighboring_boxes[bi][bj][k].epsilon_y;
                    #endif
                    
                    j = boxes[nbi][nbj];
                    // Look through all particles in box nbi, nbj and compute their contributions to rho
                    while (j != -1) {
                        xj = particles[j].x;
                        yj = particles[j].y;
                        dx = xi - (xj + epsx);
                        dy = yi - (yj + epsy);
                        density_ij = Kernel(dx, dy, rmax_squared);
                        particles[i].rho += density_ij;
                        particles[j].rho += density_ij;

                        #ifdef PFAP
                        dr = sqrt(dx*dx + dy*dy);
                        force_magnitude = Force(dr, rmax_pfap, epsilon);
                        particles[i].fx += force_magnitude * dx/dr;
                        particles[i].fy += force_magnitude * dy/dr;
                        particles[j].fx -= force_magnitude * dx/dr;
                        particles[j].fy -= force_magnitude * dy/dr;
                        #endif
                        // go to next particle
                        j = neighbors[2*j+1];
                    }
                }
                i = neighbors[2*i+1];
            }
        }
    }
    #ifdef QSAP_ZERO
    parameters[0].stopped = 1;
    #endif
    for (i=0;i<parameters[0].N;i++) {
        #ifdef QSAP_TANH
        particles[i].v = DensityDependentSpeed(particles[i].rho, parameters[0].rho_m, parameters[0].v_min, parameters[0].v_max);
        #elif defined(QSAP_EXP)
        particles[i].v = DensityDependentSpeed2(particles[i].rho, parameters[0].rho_m, parameters[0].v, parameters[0].lambda, parameters[0].phi);
        #elif defined(QSAP_ZERO_LINEAR)
        particles[i].v = DensityDependentSpeed4(particles[i].rho, parameters[0].rho_m, parameters[0].v, parameters[0].lambda, parameters[0].phi);
        #elif defined(QSAP_ZERO_SMOOTH)
        particles[i].v = DensityDependentSpeed3(particles[i].rho, parameters[0].rho_m, parameters[0].v, parameters[0].lambda, parameters[0].phi);
        #endif
        particles[i].move_x = particles[i].v*cos(particles[i].theta) + particles[i].fx;
        particles[i].move_y = particles[i].v*sin(particles[i].theta) + particles[i].fy;
        particles[i].speed = particles[i].move_x * cos(particles[i].theta) + particles[i].move_y * sin(particles[i].theta);
        #ifdef QSAP_ZERO
        if ((particles[i].move_x > EPS) || (parameters[0].move_y > EPS)) parameters[0].stopped = 0;
        #endif
    }
}

#if defined STRESS_TENSOR
// should only be on if PFAP is on
// don't need to do every timestep
void MeasureSigmaIK(double*** sigmaIK, long* neighbors, particle* particles, long** boxes, param parameters, box*** neighboring_boxes) {
    /*
    Measure the average IK stress tensor in each box by looking at particles in adjacent boxes 
    and compute their contributions.
    */
    int bi, bj, nbi, nbj, third_box_i, third_box_j;
    long i, j, k;
    int si, sj, sdim;
    double r0 = parameters.box_size;
    double rmax_pfap = parameters.interaction_range_pfap;
    double epsilon = parameters.epsilon;
    double force_magnitude;
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    double xi, yi, xj, yj, dx, dy, dr;
    double frac, horizontal_frac, vertical_frac, small_frac, big_frac;
    // Initialize stress tensor to be 0
    for (sdim=0; sdim<4; sdim++) {
        for (si=0; si<NxBox; si++) {
            for (sj=0; sj<NyBox; sj++) {
                sigmaIK[sdim][si][sj] = 0;
            }
        }
    }
    // Loop through all boxes
    for (bi = 0; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            // First particle of box
            i = boxes[bi][bj];
            while (i!=-1) {
                xi = particles[i].x;
                yi = particles[i].y;

                // particles in same box
                j = neighbors[2*i+1]; // consider pairs (i,j) with j after i
                while (j != -1) {
                    xj = particles[j].x;
                    yj = particles[j].y;
                    dx = xj - xi;
                    dy = yj - yi;
                    dr = sqrt(dx*dx + dy*dy);
                    force_magnitude = Force(dr, rmax_pfap, epsilon);
                    sigmaIK[0][bi][bj] += (force_magnitude * dx/dr) * dx / (r0*r0);
                    sigmaIK[1][bi][bj] += (force_magnitude * dx/dr) * dy / (r0*r0);
                    sigmaIK[2][bi][bj] += (force_magnitude * dy/dr) * dx / (r0*r0);
                    sigmaIK[3][bi][bj] += (force_magnitude * dy/dr) * dy / (r0*r0);
                    // go to next particle
                    j = neighbors[2*j+1];
                }

                // particles in neighboring box 1 or 3
                for (k=1; k<5; k+=2) {
                    nbi = neighboring_boxes[bi][bj][k].i;
                    nbj = neighboring_boxes[bi][bj][k].j;                
                    j = boxes[nbi][nbj];
                    // Look through all particles in box nbi, nbj and compute their contributions to rho
                    while (j != -1) {
                        xj = particles[j].x + neighboring_boxes[bi][bj][k].epsilon_x;
                        yj = particles[j].y + neighboring_boxes[bi][bj][k].epsilon_y;
                        dx = xj - xi;
                        dy = yj - yi;
                        if (k==1) frac = (neighboring_boxes[bi][bj][k].line_y - yi) / dy;
                        else frac = (neighboring_boxes[bi][bj][k].line_x - xi) / dx;
                        dr = sqrt(dx*dx + dy*dy);
                        force_magnitude = Force(dr, rmax_pfap, epsilon);
                        sigmaIK[0][bi][bj] += (frac * dx) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[1][bi][bj] += (frac * dx) * (force_magnitude * dy/dr) / (r0*r0);
                        sigmaIK[2][bi][bj] += (frac * dy) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[3][bi][bj] += (frac * dy) * (force_magnitude * dy/dr) / (r0*r0);

                        sigmaIK[0][nbi][nbj] += ((1-frac) * dx) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[1][nbi][nbj] += ((1-frac) * dx) * (force_magnitude * dy/dr) / (r0*r0);
                        sigmaIK[2][nbi][nbj] += ((1-frac) * dy) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[3][nbi][nbj] += ((1-frac) * dy) * (force_magnitude * dy/dr) / (r0*r0);

                        // go to next particle
                        j = neighbors[2*j+1];
                    }
                }

                // particles in neighboring box 2, 4
                for (k=2; k<5; k+=2) {
                    nbi = neighboring_boxes[bi][bj][k].i;
                    nbj = neighboring_boxes[bi][bj][k].j;                
                    j = boxes[nbi][nbj];
                    // Look through all particles in box nbi, nbj and compute their contributions to rho
                    while (j != -1) {
                        xj = particles[j].x + neighboring_boxes[bi][bj][k].epsilon_x;
                        yj = particles[j].y + neighboring_boxes[bi][bj][k].epsilon_y;
                        dx = xj - xi;
                        dy = yj - yi;
                        vertical_frac = (neighboring_boxes[bi][bj][k].line_y - yi) / dy;
                        horizontal_frac = (neighboring_boxes[bi][bj][k].line_x - xi) / dx;
                        dr = sqrt(dx*dx + dy*dy);
                        force_magnitude = Force(dr, rmax_pfap, epsilon);
                        small_frac = fmin(vertical_frac, horizontal_frac);
                        big_frac = fmax(vertical_frac, horizontal_frac);
                        if (vertical_frac < horizontal_frac) {
                            third_box_i = bi;
                            third_box_j = nbj;
                        }
                        else {
                            third_box_i = nbi;
                            third_box_j = bj;
                        }
                        sigmaIK[0][bi][bj] += (small_frac * dx) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[1][bi][bj] += (small_frac * dx) * (force_magnitude * dy/dr) / (r0*r0);
                        sigmaIK[2][bi][bj] += (small_frac * dy) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[3][bi][bj] += (small_frac * dy) * (force_magnitude * dy/dr) / (r0*r0);

                        sigmaIK[0][nbi][nbj] += ((1-big_frac) * dx) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[1][nbi][nbj] += ((1-big_frac) * dx) * (force_magnitude * dy/dr) / (r0*r0);
                        sigmaIK[2][nbi][nbj] += ((1-big_frac) * dy) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[3][nbi][nbj] += ((1-big_frac) * dy) * (force_magnitude * dy/dr) / (r0*r0);

                        sigmaIK[0][third_box_i][third_box_j] += ((big_frac-small_frac) * dx) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[1][third_box_i][third_box_j] += ((big_frac-small_frac) * dx) * (force_magnitude * dy/dr) / (r0*r0);
                        sigmaIK[2][third_box_i][third_box_j] += ((big_frac-small_frac) * dy) * (force_magnitude * dx/dr) / (r0*r0);
                        sigmaIK[3][third_box_i][third_box_j] += ((big_frac-small_frac) * dy) * (force_magnitude * dy/dr) / (r0*r0);
                        // go to next particle
                        j = neighbors[2*j+1];
                    }
                }
                i = neighbors[2*i+1];
            }
        }
    }
}

// compute these sigmas in bigger boxes
void MeasureSigmaActive(double** sigmaA, long* neighbors, particle* particles, long** boxes, param parameters) {
    /*
    Measure the average active stress tensor in each box by looking at particles in that box 
    and compute their contributions.
    */
    int bi, bj;
    long i;
    double r0 = parameters.box_size;
    double Dr = parameters.Dr;
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    // Loop through all boxes
    for (bi = 0; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            sigmaA[bi][bj] = 0;
            // First particle of box
            i = boxes[bi][bj];
            while (i != -1) {
                sigmaA[bi][bj] += particles[i].move_x * cos(particles[i].theta) * particles[i].v / Dr / (r0*r0);
                // go to next particle
                i = neighbors[2*i+1];
            }
        }
    }
}

void MeasureNematic(double** nematic, long* neighbors, particle* particles, long** boxes, param parameters) {
    /*
    Measure the average active stress tensor in each box by looking at particles in that box 
    and compute their contributions.
    */
    int bi, bj;
    long i;
    double r0 = parameters.box_size;
    // double Dr = parameters.Dr;
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    // Loop through all boxes
    for (bi = 0; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            nematic[bi][bj] = 0;
            // First particle of box
            i = boxes[bi][bj];
            while (i != -1) {
                nematic[bi][bj] += ( cos(particles[i].theta)*cos(particles[i].theta) - 0.5 ) / (r0*r0);
                // go to next particle
                i = neighbors[2*i+1];
            }
        }
    }
}

#ifdef WALL
void MeasureWallPressure(double* pressure_left, double* pressure_right, long* neighbors, particle* particles, long** boxes, param parameters) {
    /*
    Measure the average wall pressure by looking at particles within the potential range 
    */
    int bi, bj;
    long i;
    double wall_size = parameters.wall_size;
    double Ly = parameters.Ly;
    int num_boxes_in_wall = ceil(wall_size / parameters.box_size);
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    pressure_left[0] = 0;
    pressure_right[0] = 0;
    for (bi = 0; bi<num_boxes_in_wall; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            i = boxes[bi][bj];
            while (i != -1) {
                pressure_left[0] += particles[i].fw / Ly;
                i = neighbors[2*i+1];
            }
        }
    }
    for (bi = NxBox - num_boxes_in_wall; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            i = boxes[bi][bj];
            while (i != -1) {
                pressure_right[0] -= particles[i].fw / Ly;
                i = neighbors[2*i+1];
            }
        }
    }
}
#endif
#endif

#endif

double Distance2(long i, long j, particle* particles, param parameters) {
    /*
    Compute the squared distance between i and j
    */
    double dx, dy, dr2;
    double x = particles[i].x, y = particles[i].y;
    dx = Min(3, fabs(x-particles[j].x),
                    fabs(x-particles[j].x+parameters.Lx),
                    fabs(x-particles[j].x-parameters.Lx));
    dy = Min(3, fabs(y-particles[j].y),
                    fabs(y-particles[j].y+parameters.Ly),
                    fabs(y-particles[j].y-parameters.Ly));
    dr2 = dx*dx + dy*dy;
    return dr2;
}


void GivenInitialConditions(FILE* input, particle* particles, param* parameters, double *t
                        #ifdef HASHING
                        , long*** boxes, long** neighbors
                        #endif
                        ){
    /*
    Read initial conditions from input file. 
    Other parameters of this simul are the same as the params that generate the input file.
    */
    long i = 0;
    // double x, y, theta;
    double Lx = parameters[0].Lx;
    double Ly = parameters[0].Ly;
    #ifdef HASHING
    double box_size = parameters[0].box_size;
    int bi, bj;
    #endif
    char line[100];
    char *token;

    fprintf(parameters[0].param_file, "Initialization starting now!\n");
    fflush(parameters[0].param_file);

    fscanf(input,"%lg\n",t);
    parameters[0].next_store_time = floor(t[0]) + parameters[0].store_time_interval;
    parameters[0].next_histogram_update = floor(t[0]) + parameters[0].histogram_update_interval;
    parameters[0].next_histogram_store = floor(t[0]) + parameters[0].histogram_store_interval;

    while (fgets(line, sizeof(line), input)) {
        line[strcspn(line, "\n")] = 0;

        /* 
        don't know how many columns there are, but the first three are x,y,theta
        each column should be separated by either spaces, tabs, or both
        */
        token = strtok(line, " \t");
        particles[i].x = strtod(token, NULL);

        token = strtok(NULL, " \t");
        particles[i].y = strtod(token, NULL);

        token = strtok(NULL, " \t");
        particles[i].theta = strtod(token, NULL);

        // For debugging
        // fprintf(parameters[0].param_file, "particle %ld at %lg \t %lg\n", i, particles[i].x, particles[i].y);
        // fflush(parameters[0].param_file);

        if (particles[i].theta>=2*M_PI) particles[i].theta -= 2*M_PI;
        if (particles[i].theta<0) particles[i].theta += 2*M_PI;

        if (particles[i].x>=Lx) particles[i].x -= Lx;
        if (particles[i].x<0) particles[i].x += Lx;

        if (particles[i].y>=Ly) particles[i].y -= Ly;
        if (particles[i].y<0) particles[i].y += Ly;
        
        particles[i].fx = 0;
        particles[i].fy = 0;
        #if defined NONE || defined PFAP
        particles[i].v = parameters[0].v;
        #endif
        #ifdef HASHING
        // Store the box and neighbor information of particles
        GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        AddInBox(i, bi, bj, boxes, neighbors, particles);
        #endif
        i++;
    }
    if (i != parameters[0].N) {
        fprintf(stderr, "Input file has %ld lines of data instead of %ld\n", i, parameters[0].N);
        fflush(stderr);
        ERROR(4);
    }

    fclose(parameters[0].input_file);
    fprintf(parameters[0].param_file, "Initialization success!\n");
    fflush(parameters[0].param_file);
}

// void OriginInitialConditions(particle* particles, param parameters
//                         #ifdef HASHING
//                         , long*** boxes, long** neighbors
//                         #endif
//                         ){
//     /*
//     Put particles at random position.
//     If PFAP: reject a particle's initial position if too close to previous particles.
//     Can be used for hard core interaction at low and moderate densities.
//     */
//     long i;
//     double Lx = parameters.Lx;
//     double Ly = parameters.Ly;
//     #ifdef HASHING
//     double box_size = parameters.box_size;
//     int bi, bj;
//     #endif

//     int N = parameters.N;
//     for (i=0; i<N; i++){
//         particles[i].x = 0;
//         particles[i].y = 0;
//         particles[i].theta = 2*M_PI*genrand64_real3();
//         particles[i].fx = 0;
//         particles[i].fy = 0;
//         #if defined NONE || defined PFAP
//         particles[i].v = parameters.v;
//         #endif
//         #ifdef HASHING
//         // Store the box and neighbor information of particles
//         GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
//         AddInBox(i, bi, bj, boxes, neighbors, particles);
//         #endif
//     }
// }

// void ManualInitialConditions(particle* particles, param parameters
//                         #ifdef HASHING
//                         , long*** boxes, long** neighbors
//                         #endif
//                         ){
//     double Lx = parameters.Lx;
//     double Ly = parameters.Ly;
//     long N = parameters.N;
//     long i;
//     #ifdef HASHING
//     double box_size = parameters.box_size;
//     int bi, bj;
//     #endif
//     particles[0].x = 0.4;
//     particles[0].y = 0.2;

//     particles[1].x = 2.9;
//     particles[1].y = 2.9;

//     particles[2].x = 1.9;
//     particles[2].y = 2.2;

//     particles[3].x = 3.9;
//     particles[3].y = 0.2;

//     particles[4].x = 0.4;
//     particles[4].y = 3.7;

//     particles[5].x = 1.9;
//     particles[5].y = 1.7;

//     particles[6].x = 2.4;
//     particles[6].y = 1.7;

//     particles[7].x = 3.9;
//     particles[7].y = 3.7;

//     particles[8].x = 2.4;
//     particles[8].y = 2.2;

//     particles[9].x = 2.7;
//     particles[9].y = 2.9;
    
//     for (i=0; i<N; i++) {
//         particles[i].theta = 2*M_PI*genrand64_real3();
//         particles[i].fx = 0;
//         particles[i].fy = 0;
//         #if defined NONE || defined PFAP
//         particles[i].v = parameters.v;
//         #endif
//         #ifdef HASHING
//         // Store the box and neighbor information of particles
//         GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
//         AddInBox(i, bi, bj, boxes, neighbors, particles);
//         #endif
//     }
//     fprintf(parameters.param_file, "Done initializing boxes, start initializing densities!\n");
// }

void RandomInitialConditions(particle* particles, param parameters
                        #ifdef HASHING
                        , long*** boxes, long** neighbors
                        #endif
                        ){
    /*
    Put particles at random position.
    If PFAP: reject a particle's initial position if too close to previous particles.
    Can be used for hard core interaction at low and moderate densities.
    */
    long i;
    double Lx = parameters.Lx;
    double Ly = parameters.Ly;
    double mindist2, dr2;
    #ifdef HASHING
    double box_size = parameters.box_size;
    int bi, bj;
    #endif

    int N = parameters.N;
    for (i=0; i<N; i++){
        #ifdef PFAP
        mindist2 = -1;
        while (mindist2 < 0.9*parameters.interaction_range_pfap) {
            mindist2 = Lx*Ly;
            particles[i].x = Lx*(genrand64_real3()); // From 0 to Lx
            particles[i].y = Ly*(genrand64_real3()); // From 0 to Ly
            for (long j=1; j<i; j++) {
                dr2 = Distance2(i, j, particles, parameters);
                if (mindist2 > dr2) mindist2 = dr2;
            }
        }
        #else
        particles[i].x = Lx*(genrand64_real3()); // From 0 to Lx
        particles[i].y = Ly*(genrand64_real3()); // From 0 to Ly
        #endif

        particles[i].theta = 2*M_PI*genrand64_real3();
        particles[i].fx = 0;
        particles[i].fy = 0;
        #if defined NONE || defined PFAP
        particles[i].v = parameters.v;
        #endif
        // if (i%100 == 0)
        // fprintf(parameters.param_file, "Done initializing particle!\n");
        #ifdef HASHING
        // Store the box and neighbor information of particles
        GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        AddInBox(i, bi, bj, boxes, neighbors, particles);
        #endif
    }
    fprintf(parameters.param_file, "Done initializing boxes, start initializing densities!\n");
}

#ifdef PFAP
void LatticeInitialConditions(particle* particles, param parameters
                        #ifdef HASHING
                        , long*** boxes, long** neighbors
                        #endif
                        ){
    /*
    Create a triangular lattice with minimum distance and randomly place particles.
    */
    #ifndef WALL
    double Lx = parameters.Lx;
    #else
    // when have walls, the prescribed number of particles is true, and the density might be different.
    // the hope is that wall accumulation decreases the density and compensates for the increase because system size is smaller.
    double wall_size = parameters.wall_size;
    double Lx = parameters.Lx - 2 * wall_size;
    #endif
    double Ly = parameters.Ly;
    #ifdef HASHING
    double box_size = parameters.box_size;
    int bi, bj;
    #endif

    double rho = parameters.N / Lx / Ly * parameters.interaction_range_pfap * parameters.interaction_range_pfap;
    double lattice_size = fmin(sqrt(2/sqrt(3)/rho) - 0.01, 1);
    // horizontal lattice spacing: should keep finite when r_f is small
    double a = fmax(0.1, lattice_size * parameters.interaction_range_pfap); 
    if (a < EPS) a = 0.1 * parameters.interaction_range_qsap;
    double h = a*cos(M_PI/6.); //height between two layers. Equal to a*cos(pi/6)
    long N = parameters.N;
    int Nx = floor(Lx / a ); //Number of sites in x direction
    int Ny = floor(Ly / h );  //Number of sites in y direction
    long Nsites = Nx*Ny; //Total number of available sites in the gas phase
    int* lattice; //Lattice of gas phase
    long nb; //Number of particles placed so far
    long i,j,k;
    int bool;

    //We create the corresponding lattices, which are currently empty
    lattice = calloc(Nsites,sizeof(int));
    if (lattice==NULL) {
        fprintf(stderr, "Memory allocation for lattice failed.\n");
        fflush(stderr);
        ERROR(3);
    }

    //if there are enough spaces for the N particles
    if(N<=Nsites){
        for(i=0;i<N;i++){
            bool=1;
            // We pull a site at random and place the particle there if the
            // site is not occupied
            while(bool==1){
                k=genrand64_int64()%Nsites;
                if(lattice[k]==0){
                    lattice[k]=1;
                    bool=0;
                }
            }
        }
    }
    else {
        fprintf(stderr, "Not enough sites: %ld > %ld\n", N, Nsites);
        fflush(stderr);
        ERROR(4);
    }

    // We now place particles in particles, given the arrays
    //nb is the label of the particle to be placed. 
    nb=0;
    //We start with the particles in NsitesLiq. We go through all the
    //sites of the lattice LiquidPhase.
    for(i=0;i<Nsites;i++){
        //If it is occupied, place the particle.
        if(lattice[i]==1){
            // printf("nb = %ld N = %ld\n", nb, parameters.N);
            //we know that i = k + Nx * j
            k = i % Nx;
            j = (i-k) / Nx;
            //Compute the corresponding x. Depending on whether the line is
            //even or odd, there is an a/2 offset
            #ifndef WALL
            particles[nb].x = k*a + a * .5 * (j%2) ;
            #else
            particles[nb].x = wall_size + k*a + a * .5 * (j%2) ;
            #endif
            particles[nb].y = j*h;
            particles[nb].theta = 2*M_PI*genrand64_real2();
            particles[nb].fx = 0;
            particles[nb].fy = 0;
            #if defined NONE || defined PFAP
            particles[nb].v = parameters.v;
            #endif

            //Add particle in the good box.
            GetBox(&bi, &bj, particles[nb].x, particles[nb].y, Lx, Ly, box_size);
            AddInBox(nb, bi, bj, boxes, neighbors, particles);

            nb++;
        }
    }

    if(nb!=parameters.N) {
        fprintf(stderr, "Wrong number of particles placed\n");
        fflush(stderr);
        ERROR(4);
    }
    free(lattice);
    fflush(stderr);
}
#endif

#ifdef INIT_SLAB
void SlabRandomInitialConditions(particle* particles, param parameters
                        #ifdef HASHING
                        , long*** boxes, long** neighbors
                        #endif
                        ){
    /*
    Create two slabs of particles with different densities.
    Don't use with PFAPs.
    */

    long i;
    double Lx = parameters.Lx;
    double Ly = parameters.Ly;
    double liquid_fraction = parameters.liquid_fraction;
    #ifdef HASHING
    double box_size = parameters.box_size;
    int bi, bj;
    #endif

    long Ngas = (long) (parameters.rho_small * Lx * Ly * (1-liquid_fraction));
    long Nliquid = (long) (parameters.rho_large * Lx * Ly * liquid_fraction);

    for (i=0; i<Ngas; i++){
        particles[i].x = Lx*(1-liquid_fraction)*(-.5+genrand64_real3());
        if (particles[i].x < 0) particles[i].x += Lx;
        particles[i].y = Ly*genrand64_real3(); // From 0 to Ly

        particles[i].theta = 2*M_PI*genrand64_real3();
        particles[i].fx = 0;
        particles[i].fy = 0;
        #if defined NONE || defined PFAP
        particles[i].v = parameters.v;
        #endif
        #ifdef HASHING
        // Store the box and neighbor information of particles
        GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        AddInBox(i, bi, bj, boxes, neighbors, particles);
        #endif
    }

    for (i=Ngas; i<Ngas+Nliquid; i++){
        particles[i].x = Lx*liquid_fraction*genrand64_real3() + Lx*(0.5-liquid_fraction/2);
        particles[i].y = Ly*genrand64_real3(); // From 0 to Ly

        particles[i].theta = 2*M_PI*genrand64_real3();
        particles[i].fx = 0;
        particles[i].fy = 0;
        #if defined NONE || defined PFAP
        particles[i].v = parameters.v;
        #endif
        #ifdef HASHING
        // Store the box and neighbor information of particles
        GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        AddInBox(i, bi, bj, boxes, neighbors, particles);
        #endif
    }
}

void SlabLatticeInitialConditions(particle* particles, param parameters
                        #ifdef HASHING
                        , long*** boxes, long** neighbors
                        #endif
                        ){
    /*
    Create two triangular lattices of particles in two slabs with different densities.
    */

    #ifndef WALL
    double Lx = parameters.Lx;
    #else
    // when have walls, still use the prescribed density, will have fewer particles than if we have no walls
    double wall_size = parameters.wall_size;
    double Lx = parameters.Lx - 2 * wall_size;
    #endif
    double Ly = parameters.Ly;
    double liquid_fraction = parameters.liquid_fraction;
    #ifdef HASHING
    double box_size = parameters.box_size;
    int bi, bj;
    #endif


    double lattice_size = fmin(sqrt(2/sqrt(3)/parameters.rho_large) - 0.01, 1);
    // horizontal lattice spacing: should keep finite when r_f is small
    double a = fmax(0.3, lattice_size * parameters.interaction_range_pfap); 
    if (a < EPS) a = 0.1 * parameters.interaction_range_qsap;

    // double a = 0.9*parameters.interaction_range_pfap; //horizontal lattice spacing
    double h = a*cos(M_PI/6.); //height between two layers. Equal to a*cos(pi/6)
    long Ngas = (long) (parameters.rho_small * Lx * Ly * (1-liquid_fraction));
    long Nliquid = (long) (parameters.rho_large * Lx * Ly * liquid_fraction);
    int NxL = floor(Lx * liquid_fraction     / a ); //Number of sites in x direction in liquid phase
    int NxG = floor(Lx * (1-liquid_fraction) / a ); //Number of sites in x direction in gas phase
    int Ny = floor(Ly / h );  //Number of sites in y direction
    long NsitesGas = NxG*Ny; //Total number of available sites in the gas phase
    long NsitesLiq = NxL*Ny; //Total number of available sites in the liquid phase
    int* GasPhase; //Lattice of gas phase
    int* LiquidPhase;//Lattice of liquid phase
    long nb; //Number of particles placed so far
    long i,j,k;
    int bool;
  
    //We create the corresponding lattices, which are currently empty
    GasPhase    = calloc(NsitesGas,sizeof(int));
    LiquidPhase = calloc(NsitesLiq,sizeof(int));
    if (GasPhase==NULL || LiquidPhase==NULL) {
        fprintf(stderr, "Memory allocation for GasPhase or LiquidPhase failed.\n");
        fflush(stderr);
        ERROR(3);
    }
    //if there are enough spaces for the Ngas particles
    if(Ngas<=NsitesGas){
        for(i=0;i<Ngas;i++){
            bool=1;
            // We pull a site at random and place the particle there if the
            // site is not occupied
            while(bool==1){
                k=genrand64_int64()%NsitesGas;
                if(GasPhase[k]==0){
                    GasPhase[k]=1;
                    bool=0;
                }
            }
        }
    }
    else {
        fprintf(stderr, "Not enough gas sites: %ld > %ld\n", Ngas, NsitesGas);
        fflush(stderr);
        ERROR(4);
    }

    //if there are enough sites for the Nliquid particles
    if(Nliquid<=NsitesLiq){
        for(i=0;i<Nliquid;i++){
            bool=1;
            // We pull a site at random and place the particle there if the
            // site is not occupied
            while(bool==1){
                k=genrand64_int64()%NsitesLiq;
                // printf("particle %ld site %ld\n", i, k);
                if(LiquidPhase[k]==0){
                    LiquidPhase[k]=1;
                    bool=0;
                }
            }
            // printf("particle %ld done\n", i);
        }
    }
    else {
        fprintf(stderr, "Not enough liquid sites: %ld > %ld\n", Nliquid, NsitesLiq);
        fflush(stderr);
        ERROR(4);
    }
    
    fprintf(parameters.param_file,"Actual gas density is %lf\n", Ngas/Lx/Ly/(1-liquid_fraction));
    fprintf(parameters.param_file,"Actual liquid density is %lf\n", Nliquid/Lx/Ly/liquid_fraction);
    
    // We now place particles in particles, given the arrays
    //nb is the label of the particle to be placed. 
    nb=0;
    //We start with the particles in NsitesLiq. We go through all the
    //sites of the lattice LiquidPhase.
    for(i=0;i<NsitesLiq;i++){
        //If it is occupied, place the particle.
        if(LiquidPhase[i]==1){
            // printf("nb = %ld N = %ld\n", nb, parameters.N);
            //we know that i = k + NxL * j
            k = i % NxL;
            j = (i-k) / NxL;
            //Compute the corresponding x. The left end of the liquid phase
            //is at Lx (1-liquid_fraction)/2. Depending on whether the line is
            //even or odd, there is an a/2 offset
            #ifndef WALL
            particles[nb].x = (1-liquid_fraction)/2 * Lx + k*a + a * .5 * (j%2) ;
            #else
            particles[nb].x = wall_size + (1-liquid_fraction)/2 * Lx + k*a + a * .5 * (j%2) ;
            #endif
            particles[nb].y = j*h;
            particles[nb].theta = 2*M_PI*genrand64_real2();
            particles[nb].fx = 0;
            particles[nb].fy = 0;
            #if defined NONE || defined PFAP
            particles[nb].v = parameters.v;
            #endif

            //Add particle in the good box.
            GetBox(&bi, &bj, particles[nb].x, particles[nb].y, Lx, Ly, box_size);
            AddInBox(nb, bi, bj, boxes, neighbors, particles);

            nb++;
        }
    }

    //Let's do the same with the gas phase. x starts at (1+liquid
    //fraction)Lx but we have to use periodic boundary conditions
    for(i=0;i<NsitesGas;i++){
        if(GasPhase[i]==1){
            k = i % NxG;
            j = (i-k) / NxG;
            #ifndef WALL
            particles[nb].x = (1+liquid_fraction)/2 * Lx + k*a + a * .5 * (j%2);
            if (particles[nb].x >= Lx) particles[nb].x -= Lx;
            #else
            particles[nb].x = wall_size + (1+liquid_fraction)/2 * Lx + k*a + a * .5 * (j%2);
            if (particles[nb].x >= Lx+wall_size) particles[nb].x -= Lx;
            #endif
            particles[nb].y = j*h;
            particles[nb].theta = 2*M_PI*genrand64_real2();
            particles[nb].fx = 0;
            particles[nb].fy = 0;
            #if defined NONE || defined PFAP
            particles[nb].v = parameters.v;
            #endif

            //Add particle in the good box.
            GetBox(&bi, &bj, particles[nb].x, particles[nb].y, Lx, Ly, box_size);
            AddInBox(nb, bi, bj, boxes, neighbors, particles);

            nb++;
        }
    }

    fprintf(parameters.param_file,"Total number of particles placed %ld\n",nb);
    fprintf(parameters.param_file,"Planned number of particles %ld or %ld\n",Ngas+Nliquid,parameters.N);
    if(nb!=parameters.N) {
        fprintf(stderr, "Wrong number of particles placed\n");
        fflush(stderr);
        ERROR(4);
    }
    free(GasPhase);
    free(LiquidPhase);

    fflush(parameters.param_file);
    fflush(stderr);
}
#endif

void UpdateParticles(particle* particles, param* parameters, double* step
    #ifdef HASHING 
    ,long*** boxes, long** neighbors, box*** neighboring_boxes, double t
    #endif
    ) {   
    /*
    This function updates the positions of particles
    */ 
    long i;    
    double Lx = parameters[0].Lx;
    double Ly = parameters[0].Ly;
    double noiseamp = parameters[0].noiseamp;
    #ifdef HASHING
    double box_size = parameters[0].box_size;
    int old_bi, old_bj, new_bi, new_bj;
    #endif

    // Adaptive time step
    double fmax = 0;
    for(i = 0 ; i < parameters[0].N ; i++) {
        // printf("Fx %lf \n", particles[i].fx);
        // printf("Fy %lf \n", particles[i].fy);
        // printf("v %lf \n", particles[i].v);
        // printf("rho %lf \n", particles[i].rho);    
        if (fabs(particles[i].move_x) > fmax) fmax = fabs(particles[i].move_x);
        if (fabs(particles[i].move_y) > fmax) fmax = fabs(particles[i].move_y);
    }
    // printf("fmax %lf \n", fmax);
    #ifdef PFAP
    double ratio = parameters[0].interaction_range_pfap / (10*fmax*parameters[0].dt); // pfap radius always smaller than qsap
    #else // if only QSAP defined
    double ratio = parameters[0].interaction_range_qsap / (10*fmax*parameters[0].dt);
    #endif
    step[0] = (ratio > 1 || ratio < EPS) ? parameters[0].dt : ratio * parameters[0].dt; 
    // printf("t %lf \n", t);
    // Update all the quantities dt dependent that you use in the update
    noiseamp = parameters[0].noiseamp * sqrt(step[0]/parameters[0].dt); 
    
    for (i=0;i<parameters[0].N;i++){
        #if defined(QSAP) && !defined(HASHING)
        particles[i].rho = BruteForceDensity(particles[i].x, particles[i].y, particles, parameters);
        particles[i].v = DensityDependentSpeed(particles[i].rho, parameters.rho_m, parameters.v_min, parameters.v_max);
        #elif defined POSITION_DEPENDENT_SPEED
        particles[i].v = PositionDependentSpeed(particles[i].x, particles[i].y);
        #endif
        // printf("x = %lf, y = %lf\n", particles[i].x, particles[i].y);
        // error: particles at 0 and 100

        particles[i].x += step[0] * particles[i].move_x;
        particles[i].y += step[0] * particles[i].move_y;
        particles[i].theta += noiseamp * gasdev();

        if (particles[i].theta>=2*M_PI) particles[i].theta -= 2*M_PI;
        if (particles[i].theta<0) particles[i].theta += 2*M_PI;

        if (particles[i].x>=Lx) particles[i].x -= Lx;
        if (particles[i].x<0) particles[i].x += Lx;

        if (particles[i].y>=Ly) particles[i].y -= Ly;
        if (particles[i].y<0) particles[i].y += Ly;

        #ifdef HASHING
        old_bi = particles[i].bi;
        old_bj = particles[i].bj;
        GetBox(&new_bi, &new_bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        
        if ((new_bi != old_bi) || (new_bj != old_bj)) {
            RemoveFromBox(i, old_bi, old_bj, boxes, neighbors);
            AddInBox(i, new_bi, new_bj, boxes, neighbors, particles);
        }
        #endif
    }

#ifdef HASHING
    MeasureDensityAndForce(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
#endif

    #ifdef TESTING_DENSITY
    double rho;
    for (i=0;i<parameters.N;i++){
        rho = BruteForceDensity(particles[i].x, particles[i].y, particles, parameters);
        if (fabs(rho - particles[i].rho) > EPS) {
            printf("Rho: %lf \t %lf \n", rho, particles[i].rho);
            printf("%lf \t %d \n", t, i);
        }
    }
    #endif

    #ifdef TESTING_FORCE
    double fx, fy;
    for (i=0;i<parameters.N;i++){
        fx = 0;
        fy = 0;
        BruteForceForce(&fx, &fy, i, particles, parameters);
        if ((fabs(fx - particles[i].fx) > EPS) || 
            (fabs(fy - particles[i].fy) > EPS)) {
            printf("Forces: %lf \t %lf \t %lf \t %lf \n", fx, particles[i].fx, fy, particles[i].fy);
            printf("%lf \t %lf \t %d \n", t, step[0], i);
        }
    }
    #endif

}


void StorePositions(double t, param parameters, particle* particles) {
    /*
    This function store the positions and angles of each particles 
    every store_time_interval starting from the inputted next_store_time.
    */
    int i;
    fprintf(parameters.data_file, "%lf\n", t);

    for (i=0;i<parameters.N;i++){
        fprintf(parameters.data_file,"%lg\t%lg\t%lg\t%lg\t%lg\n",
                particles[i].x,particles[i].y,particles[i].theta,particles[i].rho,particles[i].speed);
    }
    fprintf(parameters.data_file, "\n");
    fflush(parameters.data_file);
}

#ifdef DENSITY_HISTOGRAM
void UpdateDensity(double* density_histogram, double** density_matrix, param parameters, long** boxes, long* neighbors, int store) {
    int i, j, i1, j1;
    long n0;
    long particle;
    for(i = 0; i < parameters.number_of_boxes_x; i++) {
        for(j = 0 ; j < parameters.number_of_boxes_y; j++) {
            n0 = 0;
            for(i1 = i*parameters.density_box_size ; i1 < (i+1)*parameters.density_box_size ; i1++) {
                for(j1 = j*parameters.density_box_size ; j1 < (j+1)*parameters.density_box_size; j1++){
                    particle = boxes[i1][j1];
                    while( particle != -1 ) {
                        n0++;
                        particle = neighbors[2*particle+1];
                    }
                }
            } // Counted all particles in big box
            
            if (n0 <= parameters.max_number)  
                density_histogram[n0] += 1.;
            else {
                fprintf(stderr, "max_number %ld too small, recorded n = %ld \n", parameters.max_number, n0);
                fflush(stderr);
                ERROR(4);
            }

            if (store == 1) {
                density_matrix[i][j] = n0/parameters.density_box_area;
            }
        }
    }
}

void StoreDensity(double t, param parameters, double* density_histogram, double** density_matrix, double histogram_count) {
    int i, j;
    double density, average_count;

    // Store density histogram
    fprintf(parameters.histogram_file, "%lf\n", t);
    for (i=0; i<parameters.max_number+1; i++) {
        density = i/parameters.density_box_area;
        average_count = density_histogram[i]/histogram_count;
        fprintf(parameters.histogram_file, "%lf\t%lf\n", density, average_count);  
		density_histogram[i] = 0;
    }
    fprintf(parameters.histogram_file, "\n");
    fflush(parameters.histogram_file);

    // Store density matrix
    fprintf(parameters.density_file, "%lf\n", t);
    for(j = 0; j < parameters.number_of_boxes_y; j++) {
        for(i = 0 ; i < parameters.number_of_boxes_x; i++) {
            fprintf(parameters.density_file, "%lf\t", density_matrix[i][j]);
        }
        fprintf(parameters.density_file, "\n");
    }
    fprintf(parameters.density_file, "\n");
    fflush(parameters.density_file);
}
#endif

#ifdef STRESS_TENSOR
void StoreStressTensor(double t, param parameters, double*** sigmaIK, double** sigmaA, double** nematic) {
    /*
    Store each component of the 2D stress tensor sigmaIK into four files in the folder opened.
    Called every store_position, similar to v_eff
    */
    int i, j, k;
    for (k=0; k<4; k++) {
        fprintf(parameters.sigmaIK_files[k], "%lf\n", t);
        for(j = parameters.NyBox - 1; j >= 0; j--) {
            for(i = 0 ; i < parameters.NxBox; i++) {
                fprintf(parameters.sigmaIK_files[k], "%lf\t", sigmaIK[k][i][j]);
            }
            fprintf(parameters.sigmaIK_files[k], "\n");
        }
        fprintf(parameters.sigmaIK_files[k], "\n");
        fflush(parameters.sigmaIK_files[k]);
    }

    fprintf(parameters.sigmaA_file, "%lf\n", t);
    for(j = parameters.NyBox - 1; j >= 0; j--) {
        for(i = 0 ; i < parameters.NxBox; i++) {
            fprintf(parameters.sigmaA_file, "%lf\t", sigmaA[i][j]);
        }
        fprintf(parameters.sigmaA_file, "\n");
    }
    fprintf(parameters.sigmaA_file, "\n");
    fflush(parameters.sigmaA_file);

    fprintf(parameters.nematic_file, "%lf\n", t);
    for(j = parameters.NyBox - 1; j >= 0; j--) {
        for(i = 0 ; i < parameters.NxBox; i++) {
            fprintf(parameters.nematic_file, "%lf\t", nematic[i][j]);
        }
        fprintf(parameters.nematic_file, "\n");
    }
    fprintf(parameters.nematic_file, "\n");
    fflush(parameters.nematic_file);

    fprintf(parameters.sigma_file, "%lf\n", t);
    for(j = parameters.NyBox - 1; j >= 0; j--) {
        for(i = 0 ; i < parameters.NxBox; i++) {
            fprintf(parameters.sigma_file, "%lf\t", sigmaA[i][j] + sigmaIK[0][i][j] );
        }
        fprintf(parameters.sigma_file, "\n");
    }
    fprintf(parameters.sigma_file, "\n");
    fflush(parameters.sigma_file);
}
#endif

#ifdef WALL
void StoreWallPressure(double t, param parameters, double pressure_left, double pressure_right) {
    fprintf(parameters.wall_pressure_file,"%lg\t%lg\t%lg\n", t, pressure_left, pressure_right);
    fflush(parameters.wall_pressure_file);
}
#endif
