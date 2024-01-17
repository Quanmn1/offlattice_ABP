// Pre-compiler options
// #include "./Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "math_helper.c"

#define EPS 1e-7
#define NORMALIZATION_EXP 0.46651239317833015
#ifdef HASHING
    #define HASING_ARGS , long*** boxes, long** neighbors, box*** neighboring_boxes
#else
    #define HASING_ARGS 
#endif

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
#ifdef HASHING
    int bi;
    int bj;
#endif
#ifdef QSAP
    double rho;
#endif
#ifdef PFAP
    double fx;
    double fy;
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

#ifdef HASHING
    strcat(*command_line_output, "box_size ");
    number_of_input_parameters++;
#endif   

#ifdef QSAP
    strcat(*command_line_output, "rho_m ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_min ");
    number_of_input_parameters++;

    strcat(*command_line_output, "v_max ");
    number_of_input_parameters++;

#elif defined NONE
    strcat(*command_line_output, "v ");
    number_of_input_parameters++;

#endif

    strcat(*command_line_output, "interaction_range ");
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
#ifdef HASHING
    parameters[0].box_size = strtod(argv[i], NULL); i++;
#endif
#ifdef QSAP
    parameters[0].rho_m = strtod(argv[i], NULL); i++;
    parameters[0].v_min = strtod(argv[i], NULL); i++;
    parameters[0].v_max = strtod(argv[i], NULL); i++;
#elif defined NONE
    parameters[0].v = strtod(argv[i], NULL); i++;
#endif
    parameters[0].interaction_range = strtod(argv[i], NULL); i++;
#ifdef ArbitraryKernel
    sprintf(input_parameters[0].kernel_name , "%s", argv[i]); i++;
    parameters[0].kernel_width = strtod(argv[i], NULL); i++;
#endif
#ifdef OUTPUT_DENSITY
    parameters[0].density_grid_spacing = strtod(argv[i], NULL); i++;
#endif
    sprintf(input_parameters[0].name , "%s", argv[i]); i++;
    parameters[0].final_time = strtod(argv[i], NULL); i++;
    parameters[0].next_store_time = strtod(argv[i], NULL); i++;
    parameters[0].store_time_interval = strtod(argv[i], NULL); i++;
    input_parameters[0].Dr = strtod(argv[i], NULL); i++;
#ifdef MT
    input_parameters[0].seed = (long long) strtod(argv[i], NULL); i++;
#endif
}

void AssignValues(param* parameters, inputparam input_parameters, particle** particles){
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

    // Calculate the noise amplitude for each time step
    parameters[0].noiseamp = sqrt(2 * parameters[0].dt * input_parameters.Dr);

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

#ifdef HASHING
    parameters[0].NxBox = (int) parameters[0].Lx / parameters[0].box_size;
    parameters[0].NyBox = (int) parameters[0].Ly / parameters[0].box_size;
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
    Give the value of a normalized kernel at r1 and r2. Use exponential kernel.
    */
    double r_squared = x*x + y*y;
    if (r_squared < rmax_squared) {
        double Z = NORMALIZATION_EXP * rmax_squared;
        return 1/Z * exp(-rmax_squared/(rmax_squared - r_squared));
    }
    else
        return 0;
}

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
    // Check that Lx, Ly is a multiple of box_size
    if ((fabs(parameters.Lx - NxBox*parameters.box_size) > EPS) || 
        (fabs(parameters.Ly - NyBox*parameters.box_size) > EPS)) {
        printf("Length not an integer multiple of size of box");
        exit(2);
    }
    // Construct and initialize boxes
    boxes[0] = (long**) malloc(NxBox * sizeof(long*));
    for ( i = 0; i < NxBox; i++) {
        boxes[0][i] = (long*) malloc(NyBox * sizeof(long));
        for ( j = 0; j < NyBox; j++) {
            boxes[0][i][j] = -1;
        }
    }
}

void ConstructNeighbors(long** neighbors, long N) {
    // set all entries to -1: no relation yet between particles
	neighbors[0] = (long*) malloc((2*N) * sizeof(long));
	for( int i = 0 ; i < 2*N; i++)
		neighbors[0][i] = -1;
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
            exit(5);
        }
        // The next of the previous becomes the current next
        neighbors[0][2*prev+1] = next;
        // If next is a particle, its prev becomes the current prev
        if (next != -1) neighbors[0][2*next] = prev;
    }
}

void FreeBoxes(int NxBox, long*** boxes) {
    for (int i = 0; i < NxBox; i++)
        free(boxes[0][i]);
    free(boxes[0]);
}

#ifdef QSAP
void ConstructNeighboringBoxes(param parameters, box*** neighboring_boxes) {
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
            neighboring_boxes[bi][bj][4].epsilon_y = (bj==0) ? - Ly : 0;
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
    int NxBox = parameters.NxBox;
    int NyBox = parameters.NyBox;
    double xi, yi, xj, yj, dx, dy;
    double density_ij;
    // Initialize density fieds to zero for all particles
    for (i=0; i<parameters.N; i++) {
        particles[i].rho = 0;
    }
    // Loop through all boxes
    for (bi = 0; bi<NxBox; bi++) {
        for (bj = 0; bj<NyBox; bj++) {
            // First particle of box
            i = boxes[bi][bj];
            while (i!=-1) {
                xi = particles[i].x;
                yi = particles[i].y;

                // Particle in neighboring boxes
                for (k=0; k<5; k++) {
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
                        if (k > 0) particles[j].rho += density_ij;
                        // go to next particle
                        j = neighbors[2*j+1];
                    }
                }
                i = neighbors[2*i+1];
            }
        }
    }
}

double ComputeForces_withSpatialHashing() {
    return 0;
// try this by myself, might do with Gianmarco
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
        particles[i].theta = 2*M_PI*genrand64_real3(); //M_2_PI calls 2*Pi in C
        #ifdef HASHING
        // Store the box and neighbor information of particles
        GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        AddInBox(i, bi, bj, boxes, neighbors, particles);
        #endif
    }
    #ifdef HASHING
    MeasureDensity(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
    #endif
}

void InitialConditionsOrigin(particle* particles, param parameters, long*** boxes, long** neighbors){
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
        particles[i].x = 0;
        particles[i].y = 0;
        particles[i].theta = 2*M_PI*genrand64_real3();
        #ifdef HASHING
        GetBox(&bi, &bj, particles[i].x, particles[i].y, Lx, Ly, box_size);
        AddInBox(i, bi, bj, boxes, neighbors, particles);
        #endif

    }
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

double BruteForceDensity(double x, double y, particle* particles, param parameters) {
    /*
    Calculate density at point (x,y) using the kernel specified in parameters
    */
    int j;
    double relative_x, relative_y;
    double rmax_squared = parameters.interaction_range*parameters.interaction_range;
    double rho = 0;
    for (j=0;j<parameters.N;j++){
        relative_x = Min(3, fabs(x-particles[j].x),
                        fabs(x-particles[j].x+parameters.Lx),
                        fabs(x-particles[j].x-parameters.Lx));
        relative_y = Min(3, fabs(y-particles[j].y),
                        fabs(y-particles[j].y+parameters.Ly),
                        fabs(y-particles[j].y-parameters.Ly));
        rho += Kernel(relative_x, relative_y, rmax_squared);
    }

    return rho;
}

#ifdef QSAP
double DensityDependentSpeed(double rho, double rho_m, double v_min, double v_max){
    /*
    Give the density-dependent speed (formula in the MIPS review paper)
    */
    return v_max + (v_min - v_max)/2 * (1+tanh(2*rho/rho_m-2));
}

#elif defined POSITION_DEPENDENT_SPEED
double PositionDependentSpeed(double x, double y){
    /*
    Give a quadrant-dependent speed 
    */
    if (x>=0 && y>=0) return 0.5;
    else if (x>0 && y<0) return 0.3;
    else if (x<0 && y>0) return 0.2;
    else return 0.1;
}

#endif

void UpdateParticles(particle* particles, param parameters
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
    #ifdef HASHING
    double box_size = parameters.box_size;
    int old_bi, old_bj, new_bi, new_bj;
    #else 
    double rho;
    #endif

    for (i=0;i<parameters.N;i++){
        #ifdef QSAP
            #ifdef HASHING
            v = DensityDependentSpeed(particles[i].rho, parameters.rho_m, parameters.v_min, parameters.v_max);
            #else
            rho = BruteForceDensity(particles[i].x, particles[i].y, particles, parameters);
            v = DensityDependentSpeed(rho, parameters.rho_m, parameters.v_min, parameters.v_max);
            #endif
        #elif defined NONE
        v = parameters.v;
        #elif defined POSITION_DEPENDENT_SPEED
        v = PositionDependentSpeed(particles[i].x, particles[i].y);
        #endif

        #ifdef HASHING
        old_bi = particles[i].bi;
        old_bj = particles[i].bj;
        #endif

        particles[i].x += parameters.dt * v * cos(particles[i].theta);
        particles[i].y += parameters.dt * v * sin(particles[i].theta);
        particles[i].theta += parameters.noiseamp * gasdev();
        if (particles[i].theta>2*M_PI)
            particles[i].theta -= 2*M_PI;
        if (particles[i].theta<0)
            particles[i].theta += 2*M_PI;

        if (particles[i].x>Lx/2)
            particles[i].x -= Lx;
        if (particles[i].x<-Lx/2)
            particles[i].x += Lx;

        if (particles[i].y>Ly/2)
            particles[i].y -= Ly;
        if (particles[i].y<-Ly/2)
            particles[i].y += Ly;

        #ifdef HASHING
        GetBox(&new_bi, &new_bj, particles[i].x, particles[i].y, Lx, Ly, box_size);

        if ((new_bi != old_bi) || (new_bj != old_bj)) {
            RemoveFromBox(i, old_bi, old_bj, boxes, neighbors);
            AddInBox(i, new_bi, new_bj, boxes, neighbors, particles);
        }
        #endif
    }
    #ifdef HASHING
    MeasureDensity(neighbors[0], particles, boxes[0], parameters, neighboring_boxes);
    #endif

    #ifdef TESTING_DENSITY
    double rho;
    for (i=0;i<parameters.N;i++){
        rho = Density(particles[i].x, particles[i].y, particles, parameters);
        if (fabs(rho - particles[i].rho) > EPS) {
            printf("%lf \n", rho - particles[i].rho);
            printf("%lf \n", t);
            printf("%d \n", i);
        }
    }
    #endif

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

void StorePositions(double t, param *parameters, particle* particles
    #ifdef HASHING
    ,long*** boxes, long** neighbors
    #endif
    ) {
    /*
    This function store the positions and angles of each particles 
    every store_time_interval starting from the inputted next_store_time.
    */
    int i;
    if (t >= parameters[0].next_store_time){
        #ifdef OUTPUT_DENSITY
        // Print out the density profile to density file
        int i, j;
        for (i=-parameters[0].half_number_of_points_x;i<=parameters[0].half_number_of_points_x;i++) {
            for (j=-parameters[0].half_number_of_points_y; j<=parameters[0].half_number_of_points_y;j++) {
                // Maybe output density from each particle instead of calculating again on a new grid?
                double rho = Density(i*parameters[0].density_grid_spacing, 
                                    j*parameters[0].density_grid_spacing, particles, parameters[0]);
                fprintf(parameters[0].density_file, "%lg \t", rho);
            }
            fprintf(parameters[0].density_file, "\n");
        }
        fprintf(parameters[0].density_file, "\n");
        fflush(parameters[0].density_file);

        // Print out the MSD to data file
        double x2 = 0, y2 = 0;
        MeanSquareDisplacement(parameters[0].N, particles, &x2, &y2);
        fprintf(parameters[0].data_file,"%lg \t %lg \t", x2, y2); // First two entries will be x2 and y2
        #endif

        #ifdef HASHING
        for (i=0;i<parameters[0].N;i++){
            fprintf(parameters[0].data_file,"%lg \t %d \t %lg \t %lg \t %lg \t %lg \t %d \t %d \t",
                    t,i,particles[i].x,particles[i].y,particles[i].theta,particles[i].rho,particles[i].bi,particles[i].bj);
        }
        #else
        for (i=0;i<parameters[0].N;i++){
            fprintf(parameters[0].data_file,"%lg \t %d \t %lg \t %lg \t %lg \t",
                    t,i,particles[i].x,particles[i].y,particles[i].theta);
        }
        #endif
        fprintf(parameters[0].data_file, "\n");

        #ifdef TESTING_BOX
        // Print out boxes and neighbors to check the logic
        for (int i = parameters[0].NyBox-1; i >= 0; i--) {
            for (int j=0; j<parameters[0].NxBox;j++) {
                fprintf(parameters[0].boxes_file, "%ld \t", boxes[0][j][i]);
            }
            fprintf(parameters[0].boxes_file, "\n");
        }
        for (int i = 0; i < parameters[0].N*2; i++) {
            fprintf(parameters[0].boxes_file, "%ld \t", neighbors[0][i]);
        }
        fprintf(parameters[0].boxes_file, "\n");
        fflush(parameters[0].boxes_file);
        #endif
    
        parameters[0].next_store_time += parameters[0].store_time_interval;
    }
    
    fflush(parameters[0].data_file);
}