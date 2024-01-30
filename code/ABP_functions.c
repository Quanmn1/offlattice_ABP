// Pre-compiler options
// #include "./Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "math_helper.c"

#define EPS 1e-7
#define NORMALIZATION_EXP 0.46651239317833015

/*
store parameters used in the loop
*/
typedef struct param {
    double dt;
    int N;
    double Lx;
    double Ly;
#ifdef QSAP
    double rho_m;
    double v_max;
    double v_min;
#elif defined PFAP
    double epsilon;
    double v;
    // double mu;   // mu=1
#elif defined NONE
    double v;
#endif
    double interaction_range;
#ifdef ArbitraryKernel
    double (*kernel)(double);
    double kernel_width;
    double kernel_normalization;
#endif
    double final_time;
    double next_store_time;
    double store_time_interval;
    double noiseamp;
#ifdef HASHING
    double box_size;
    int NxBox;
    int NyBox;
    #ifdef DENSITY_HISTOGRAM
    int density_box_size; // how many hashing boxes are in a density box
    double density_box_area;
    long max_number;
    FILE* histogram_file;
    double next_histogram_update;
    double histogram_update_interval;
    double next_histogram_store;
    double histogram_store_interval;
    #endif
#endif
#ifdef OUTPUT_DENSITY
    double density_grid_spacing;
    int half_number_of_points_x;
    int half_number_of_points_y;
    FILE* density_file;
#endif
#ifdef TESTING
    FILE* boxes_file;
#endif
    FILE* param_file;
    FILE* data_file;
} param;

/*
Store parameters used only at the beginning
*/
typedef struct inputparam {
    char name[1000];
    double Dr;
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
#ifdef HASHING
    int bi;
    int bj;
#endif
} particle;

#ifdef HASHING
typedef struct box {
    int i;
    int j;
#ifdef PBC
    double epsilon_x; // For PBC
    double epsilon_y;
#endif
} box;
#endif

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

    strcat(*command_line_output, "N ");
    number_of_input_parameters++;

    strcat(*command_line_output, "Lx ");
    number_of_input_parameters++;

    strcat(*command_line_output, "Ly ");
    number_of_input_parameters++;

#ifdef QSAP
    strcat(*command_line_output, "rho_m ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_min ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_max ");
    number_of_input_parameters++;

#elif defined PFAP
    strcat(*command_line_output, "v ");
    number_of_input_parameters++;

    strcat(*command_line_output, "epsilon ");
    number_of_input_parameters++;

    // strcat(*command_line_output, "mu ");
    // number_of_input_parameters++;

#elif defined NONE
    strcat(*command_line_output, "v ");
    number_of_input_parameters++;

#endif

    strcat(*command_line_output, "Dr ");
    number_of_input_parameters++;

    strcat(*command_line_output, "r_max ");
    number_of_input_parameters++;

    strcat(*command_line_output, "final_time ");
    number_of_input_parameters++;
    
#ifdef ArbitraryKernel
    strcat(*command_line_output, "KernelName ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelWidth ");
    number_of_input_parameters++;
#endif

#ifdef OUTPUT_DENSITY
    strcat(*command_line_output, "DensityGridSpacing ");
    number_of_input_parameters++;
#endif

#ifdef HASHING
    // strcat(*command_line_output, "box_size ");
    // number_of_input_parameters++;

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
    parameters[0].N = (long) strtod(argv[i], NULL); i++;
    parameters[0].Lx = strtod(argv[i], NULL); i++;
    parameters[0].Ly = strtod(argv[i], NULL); i++;
#ifdef QSAP
    parameters[0].rho_m = strtod(argv[i], NULL); i++;
    parameters[0].v_min = strtod(argv[i], NULL); i++;
    parameters[0].v_max = strtod(argv[i], NULL); i++;
#elif defined PFAP
    parameters[0].v = strtod(argv[i], NULL); i++;
    parameters[0].epsilon = strtod(argv[i], NULL); i++;
    // parameters[0].mu = strtod(argv[i], NULL); i++;
#elif defined NONE
    parameters[0].v = strtod(argv[i], NULL); i++;
#endif
    input_parameters[0].Dr = strtod(argv[i], NULL); i++;
    parameters[0].interaction_range = strtod(argv[i], NULL); i++;
    parameters[0].final_time = strtod(argv[i], NULL); i++;

#ifdef ArbitraryKernel
    sprintf(input_parameters[0].kernel_name , "%s", argv[i]); i++;
    parameters[0].kernel_width = strtod(argv[i], NULL); i++;
#endif
#ifdef OUTPUT_DENSITY
    parameters[0].density_grid_spacing = strtod(argv[i], NULL); i++;
#endif

#ifdef HASHING
    // parameters[0].box_size = strtod(argv[i], NULL); i++;
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
#ifdef MT
    input_parameters[0].seed = (long long) strtod(argv[i], NULL); i++;
#endif
}

void AssignValues(param* parameters, inputparam input_parameters, particle** particles, double** density_histogram){
    /*
    This function allocates memory to the particles array, and calculate some parameters
    */

    char filename[1010]; // string of the filenames
    // alloc a block of mem to the pointer (array) particles
    particles[0] = malloc(sizeof(particle) * parameters[0].N);
    if (particles[0] == NULL) {
        printf("Memory allocation for particles failed.\n");
        exit(3);
    }

    // Assign the files
    sprintf(filename, "%s_param",input_parameters.name); // print name into filename
    parameters[0].param_file = fopen(filename, "w");

    sprintf(filename, "%s_data",input_parameters.name);
    parameters[0].data_file = fopen(filename, "w");

    // Calculate the noise amplitude for each time step
    parameters[0].noiseamp = sqrt(2 * parameters[0].dt * input_parameters.Dr);

#ifdef HASHING
    parameters[0].box_size = parameters[0].interaction_range;

    parameters[0].NxBox = (int) parameters[0].Lx / parameters[0].box_size;
    parameters[0].NyBox = (int) parameters[0].Ly / parameters[0].box_size;

    // Check that Lx, Ly is a multiple of box_size
    if ((fabs(parameters[0].Lx - parameters[0].NxBox*parameters[0].box_size) > EPS) || 
        (fabs(parameters[0].Ly - parameters[0].NyBox*parameters[0].box_size) > EPS)) {
        printf("Length not an integer multiple of size of box");
        exit(2);
    }

    #ifdef DENSITY_HISTOGRAM
    sprintf(filename, "%s_histogram",input_parameters.name);
    parameters[0].histogram_file = fopen(filename, "w");

    int number_of_boxes_x = (int) parameters[0].NxBox / parameters[0].density_box_size;
    int number_of_boxes_y = (int) parameters[0].NyBox / parameters[0].density_box_size;
    if ((abs(parameters[0].NxBox - number_of_boxes_x*parameters[0].density_box_size) > EPS) || 
        (abs(parameters[0].NyBox - number_of_boxes_y*parameters[0].density_box_size) > EPS)) {
        printf("System does not contain integer numbers of density boxes");
        exit(2);
    }

    parameters[0].density_box_area = parameters[0].density_box_size*parameters[0].density_box_size * \
                                     parameters[0].interaction_range*parameters[0].interaction_range;
    parameters[0].max_number = parameters[0].N / (number_of_boxes_x*number_of_boxes_y) * 10;

    density_histogram[0] = malloc(parameters[0].max_number * sizeof(double));
    if (density_histogram[0] == NULL) {
        printf("Memory allocation for density_histogram failed.\n");
        exit(3);
    }
    #endif

#endif

#ifdef OUTPUT_DENSITY
    // the density grid is centered at (0,0), and the half_number_of_points are the number of points to one side
    parameters[0].half_number_of_points_x = floor(parameters[0].Lx/2 / parameters[0].density_grid_spacing);
    parameters[0].half_number_of_points_y = floor(parameters[0].Ly/2 / parameters[0].density_grid_spacing);

    sprintf(filename, "%s_density",input_parameters.name);
    parameters[0].density_file = fopen(filename, "w");
    fprintf(parameters[0].density_file, "%d \t %d \t %lg \n", parameters[0].half_number_of_points_x,
            parameters[0].half_number_of_points_y, parameters[0].density_grid_spacing);
#endif

#ifdef TESTING
    sprintf(filename, "%s_boxes",input_parameters.name);
    parameters[0].boxes_file = fopen(filename, "w");
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
        printf("Kernel not supported! Supported kernel: tanh, exp");
        exit(1);
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

    fflush(parameters.param_file); // empty the buffer into the file
    
    fclose(parameters.param_file);
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

double PositionDependentSpeed(double x, double y){
    /*
    Give a quadrant-dependent speed 
    */
    if (x>=0 && y>=0) return 0.5;
    else if (x>0 && y<0) return 0.3;
    else if (x<0 && y>0) return 0.2;
    else return 0.1;
}

double Force(double r, double rmax, double epsilon) {
    /*
    Return the magnitude of the force at distance r
    */
    double ratio = rmax / r;
    if (ratio >= 1)
        return 12*epsilon/rmax * (pow(ratio, 13) - pow(ratio, 7));
    else
        return 0;
}

double BruteForceDensity(double x, double y, particle* particles, param parameters) {
    /*
    Calculate density at point (x,y) using the kernel specified in parameters
    */
    int j;
    double dx, dy;
    double rmax_squared = parameters.interaction_range*parameters.interaction_range;
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

void BruteForceForce(double* fx, double* fy, long i, particle* particles, param parameters) {
    /*
    Calculate force at point (x,y) using the kernel specified in parameters
    */
    int j;
    double dx, dy, dr;
    double rmax = parameters.interaction_range;
    double epsilon = parameters.epsilon;
    double force_magnitude;
    double x = particles[i].x;
    double y = particles[i].y;
    fx[0] = 0;
    fy[0] = 0;
    for (j=0;j<parameters.N;j++){
        if (j != i) {
            dx = AbsMin(3, x-particles[j].x,
                            x-particles[j].x+parameters.Lx,
                            x-particles[j].x-parameters.Lx);
            dy = AbsMin(3, y-particles[j].y,
                            y-particles[j].y+parameters.Ly,
                            y-particles[j].y-parameters.Ly);
            dr = sqrt(dx*dx + dy*dy);
            force_magnitude = Force(dr, rmax, epsilon);
            fx[0] += force_magnitude * dx/dr;
            fy[0] += force_magnitude * dy/dr;
        }
    }
}


#ifdef OUTPUT_DENSITY
void MeanSquareDisplacement(int N, particle* particles, double* x2, double* y2){
    int i;
    for (i=0;i<N;i++) {
        *x2 += particles[i].x * particles[i].x / N;
        *y2 += particles[i].y * particles[i].y / N;
    }
}
#endif

#ifdef HASHING
void GetBox(int* bi, int* bj, double x, double y, double Lx, double Ly, double box_size) {
    // Get the box indices from coordinates
    bi[0] = floor( (x+Lx/2) / box_size);
    bj[0] = floor( (y+Ly/2) / box_size);
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
        printf("Memory allocation for boxes failed.\n");
        exit(3);
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
        printf("Memory allocation for neighbors failed.\n");
        exit(3);
    }
	for( int i = 0 ; i < 2*N; i++)
		neighbors[0][i] = -1;
}

void FreeBoxes(int NxBox, long*** boxes) {
    for (int i = 0; i < NxBox; i++)
        free(boxes[0][i]);
    free(boxes[0]);
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

void AddInBox(int index_particle, int bi, int bj, long*** boxes, long** neighbors, particle* particles) {
    int k = boxes[0][bi][bj]; // Save the previous first particle of the box
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
            printf("boxes and neighbors uncompatible: %d %d", bi, bj);
            exit(4);
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

            neighboring_boxes[bi][bj][1].i = bi;
            neighboring_boxes[bi][bj][1].j = (bj+1)%NyBox;
            neighboring_boxes[bi][bj][1].epsilon_x = 0;
            neighboring_boxes[bi][bj][1].epsilon_y = (bj==NyBox-1) ? Ly : 0;

            neighboring_boxes[bi][bj][2].i = (bi+1)%NxBox;
            neighboring_boxes[bi][bj][2].j = (bj+1)%NyBox;
            neighboring_boxes[bi][bj][2].epsilon_x = (bi==NxBox-1) ? Lx : 0;
            neighboring_boxes[bi][bj][2].epsilon_y = (bj==NyBox-1) ? Ly : 0;

            neighboring_boxes[bi][bj][3].i = (bi+1)%NxBox;
            neighboring_boxes[bi][bj][3].j = bj;
            neighboring_boxes[bi][bj][3].epsilon_x = (bi==NxBox-1) ? Lx : 0;
            neighboring_boxes[bi][bj][3].epsilon_y = 0;

            neighboring_boxes[bi][bj][4].i = (bi+1)%NxBox;
            neighboring_boxes[bi][bj][4].j = (bj+NyBox-1)%NyBox;
            neighboring_boxes[bi][bj][4].epsilon_x = (bi==NxBox-1) ? Lx : 0;
            neighboring_boxes[bi][bj][4].epsilon_y = (bj==0) ? -Ly : 0;
        }
    }
}

void MeasureDensity(long* neighbors, particle* particles, long** boxes, param parameters, box*** neighboring_boxes) {
    /*
    Measure the density around each particle by looking at particles in adjacent boxes 
    and compute their contributions.
    */
    int bi, bj, nbi, nbj;
    long i, j, k;
    double epsx=0, epsy=0;
    double rmax_squared = parameters.interaction_range*parameters.interaction_range;
    double Z = NORMALIZATION_EXP * rmax_squared;
    double self_density = 1/Z/M_E;
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    double xi, yi, xj, yj, dx, dy;
    double density_ij;
    // Initialize density fieds to zero for all particles
    for (i=0; i<parameters.N; i++) {
        particles[i].rho = self_density;
    }
    // Loop through all boxes
    for (bi = 0; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            // First particle of box
            i = boxes[bi][bj];
            while (i!=-1) {
                xi = particles[i].x;
                yi = particles[i].y;

                j = neighbors[2*i+1]; // consider pairs (i,j) with j including and after i
                while (j != -1) {
                    xj = particles[j].x;
                    yj = particles[j].y;
                    dx = xi - xj;
                    dy = yi - yj;
                    density_ij = Kernel(dx, dy, rmax_squared);
                    particles[i].rho += density_ij;
                    particles[j].rho += density_ij;
                    // go to next particle
                    j = neighbors[2*j+1];
                }

                // Particle in neighboring boxes
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
                        // go to next particle
                        j = neighbors[2*j+1];
                    }
                }
                i = neighbors[2*i+1];
            }
        }
    }
    #ifdef QSAP
    for (i=0;i<parameters.N;i++) {
        particles[i].v = DensityDependentSpeed(particles[i].rho, parameters.rho_m, parameters.v_min, parameters.v_max);
    }
    #endif
}

#ifdef PFAP
void MeasureForce(long* neighbors, particle* particles, long** boxes, param parameters, box*** neighboring_boxes) {
    /*
    Measure the force each particle experiences.
    */
    int bi, bj, nbi, nbj;
    long i, j, k;
    double epsx=0, epsy=0;
    double rmax = parameters.interaction_range;
    double epsilon = parameters.epsilon;
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    double xi, yi, xj, yj, dx, dy, dr;
    double force_magnitude;
    // Initialize density fieds to zero for all particles
    for (i=0; i<parameters.N; i++) {
        particles[i].fx = 0;
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

                j = neighbors[2*i+1]; // consider pairs (i,j) with j after i
                while (j != -1) {
                    xj = particles[j].x;
                    yj = particles[j].y;
                    dx = xi - xj;
                    dy = yi - yj;
                    dr = sqrt(dx*dx + dy*dy);
                    force_magnitude = Force(dr, rmax, epsilon);
                    particles[i].fx += force_magnitude * dx/dr;
                    particles[i].fy += force_magnitude * dy/dr;
                    particles[j].fx -= force_magnitude * dx/dr;
                    particles[j].fy -= force_magnitude * dy/dr;
                    // go to next particle
                    j = neighbors[2*j+1];
                }

                // Particle in neighboring boxes
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
                        dr = sqrt(dx*dx + dy*dy);
                        force_magnitude = Force(dr, rmax, epsilon);
                        particles[i].fx += force_magnitude * dx/dr;
                        particles[i].fy += force_magnitude * dy/dr;
                        particles[j].fx -= force_magnitude * dx/dr;
                        particles[j].fy -= force_magnitude * dy/dr;
                        // go to next particle
                        j = neighbors[2*j+1];
                    }
                }
                i = neighbors[2*i+1];
            }
        }
    }
    
}
#endif
#endif

void InitialConditions(particle* particles, param parameters
                        #ifdef HASHING
                        , long*** boxes, long** neighbors, box*** neighboring_boxes
                        #endif
                        ){
    /*
    This function initializes the positions and angles of each particle 
    to a uniformly random number
    */

    // not modifying particles, just where they point to, so don't need pointers
    int i;
    double Lx = parameters.Lx;
    double Ly = parameters.Ly;
    #ifdef HASHING
    double box_size = parameters.box_size;
    int bi, bj;
    #endif
    for (i=0; i<parameters.N; i++){
        particles[i].x = Lx*(-.5+genrand64_real3()); // From -Lx/2 to Lx/2
        particles[i].y = Ly*(-.5+genrand64_real3()); // From -Ly/2 to Ly/2
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
    #ifdef HASHING
    MeasureDensity(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
    #if defined PFAP
        MeasureForce(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
    #endif
    #endif
}


void UpdateParticles(particle* particles, param parameters, double* step
    #ifdef HASHING 
    ,long*** boxes, long** neighbors, box*** neighboring_boxes, double t 
    #endif
    ) {   
    /*
    This function updates the positions of particles
    */ 
    int i;    
    double Lx = parameters.Lx;
    double Ly = parameters.Ly;
    double v;
    double noiseamp = parameters.noiseamp;
    #ifdef HASHING
    double box_size = parameters.box_size;
    int old_bi, old_bj, new_bi, new_bj;
    #endif

    // Adaptive time step
    double fmax = 0;
    for(i = 0 ; i < parameters.N ; i++) {
        // printf("Fx %lf \n", particles[i].fx);
        // printf("Fy %lf \n", particles[i].fy);
        particles[i].move_x = particles[i].v*cos(particles[i].theta) + particles[i].fx;
        particles[i].move_y = particles[i].v*sin(particles[i].theta) + particles[i].fy;    
        if (fabs(particles[i].move_x) > fmax) fmax = fabs(particles[i].move_x);
        if (fabs(particles[i].move_y) > fmax) fmax = fabs(particles[i].move_y);
    }
    double ratio = parameters.interaction_range / (10*fmax*parameters.dt);
    step[0] = (ratio > 1) ? parameters.dt : ratio * parameters.dt; 
    // Update all the quantities dt dependent that you use in the update
    noiseamp = parameters.noiseamp * sqrt(step[0]/parameters.dt); 

    for (i=0;i<parameters.N;i++){
        #if defined(QSAP) && !defined(HASHING)
        particles[i].rho = BruteForceDensity(particles[i].x, particles[i].y, particles, parameters);
        particles[i].v = DensityDependentSpeed(particles[i].rho, parameters.rho_m, parameters.v_min, parameters.v_max);
        #elif defined POSITION_DEPENDENT_SPEED
        particles[i].v = PositionDependentSpeed(particles[i].x, particles[i].y);
        #endif

        particles[i].x += step[0] * particles[i].move_x;
        particles[i].y += step[0] * particles[i].move_y;
        particles[i].theta += noiseamp * gasdev();

        if (particles[i].theta>2*M_PI) particles[i].theta -= 2*M_PI;
        if (particles[i].theta<0) particles[i].theta += 2*M_PI;

        if (particles[i].x>Lx/2) particles[i].x -= Lx;
        if (particles[i].x<-Lx/2) particles[i].x += Lx;

        if (particles[i].y>Ly/2) particles[i].y -= Ly;
        if (particles[i].y<-Ly/2) particles[i].y += Ly;

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
    MeasureDensity(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
    #if defined PFAP
        MeasureForce(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
    #endif
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

#ifdef DENSITY_HISTOGRAM
void UpdateHistogram(double* density_histogram, param parameters, long** boxes, long* neighbors) {
    int i, j, i1, j1;
    long n0;
    long particle;
    for(i = 0; i < parameters.NxBox; i+= parameters.density_box_size) {
        for(j = 0 ; j < parameters.NyBox; j+= parameters.density_box_size) {
            n0 = 0;
            for(i1 = i ; i1 < i + parameters.density_box_size ; i1++) {
                for(j1 = j ; j1 < j+ parameters.density_box_size; j1++){
                    particle = boxes[i1][j1];
                    while( particle != -1 ) {
                        n0++;
                        particle = neighbors[2*particle+1];
                    }
                }
            } // Counted all particles in big box
            
            if (n0 < parameters.max_number)  
                density_histogram[n0] += 1.;
            else {
                /* Maybe  you can relax it and adapt live the histograms*/
                printf("max_number too small, recorded n = %ld \n", n0);
                exit(4);
            }
        }
    }
}

void StoreHistogram(double t, param parameters, double* density_histogram, double histogram_count, long*** boxes, long** neighbors) {
    /*
    Store the density histogram
    */
    int i;
    double density, average_count;
    fprintf(parameters.histogram_file, "%lf\n", t);
    for (i=0; i<parameters.max_number; i++) {
        density = i/parameters.density_box_area;
        average_count = density_histogram[i]/histogram_count;
        fprintf(parameters.histogram_file, "%lf \t %lf \n", density, average_count);  
		density_histogram[i] = 0;
    }
    fprintf(parameters.histogram_file, "\n");
    fflush(parameters.histogram_file);
}
#endif

void StorePositions(double t, param parameters, particle* particles
    #ifdef HASHING
    ,long*** boxes, long** neighbors
    #endif
    ) {
    /*
    This function store the positions and angles of each particles 
    every store_time_interval starting from the inputted next_store_time.
    */
    int i;
    fprintf(parameters.data_file, "%lf\n", t);

    #ifdef OUTPUT_DENSITY
    // Print out the density profile to density file
    int i, j;
    for (i=-parameters.half_number_of_points_x;i<=parameters.half_number_of_points_x;i++) {
        for (j=-parameters.half_number_of_points_y; j<=parameters.half_number_of_points_y;j++) {
            // Maybe output density from each particle instead of calculating again on a new grid?
            double rho = Density(i*parameters.density_grid_spacing, 
                                j*parameters.density_grid_spacing, particles, parameters);
            fprintf(parameters.density_file, "%lg \t", rho);
        }
        fprintf(parameters.density_file, "\n");
    }
    fprintf(parameters.density_file, "\n");
    fflush(parameters.density_file);

    // Print out the MSD to data file
    double x2 = 0, y2 = 0;
    MeanSquareDisplacement(parameters.N, particles, &x2, &y2);
    fprintf(parameters.data_file,"%lg \t %lg \t", x2, y2); // First two entries will be x2 and y2
    #endif

    for (i=0;i<parameters.N;i++){
        fprintf(parameters.data_file,"%lg \t %lg \t %lg \t %lg \n",
                particles[i].x,particles[i].y,particles[i].theta,particles[i].rho);
    }
    fprintf(parameters.data_file, "\n");
    fflush(parameters.data_file);

    #ifdef TESTING_BOX
    // Print out boxes and neighbors to check the logic
    for (int i = parameters.NyBox-1; i >= 0; i--) {
        for (int j=0; j<parameters.NxBox;j++) {
            fprintf(parameters.boxes_file, "%ld \t", boxes[0][j][i]);
        }
        fprintf(parameters.boxes_file, "\n");
    }
    for (int i = 0; i < parameters.N*2; i++) {
        fprintf(parameters.boxes_file, "%ld \t", neighbors[0][i]);
    }
    fprintf(parameters.boxes_file, "\n");
    fflush(parameters.boxes_file);
    #endif
}