// Pre-compiler options
// #include "./Options.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "math_helper.c"

#define EPS 1e-7

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
    double (*kernel)(double);
    double kernel_width;
    double kernel_normalization;
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
    char kernel_name[10];
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
} particle;

#ifdef HASHING
typedef struct box {
    int i;
    int j;
    double epsilon_x; // For PBC
    double epsilon_y;
} box; // only need for neighboring boxes
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

    strcat(*command_line_output, "KernelName ");
    number_of_input_parameters++;

    strcat(*command_line_output, "KernelWidth ");
    number_of_input_parameters++;

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
    parameters[0].N = (int) strtod(argv[i], NULL); i++;
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
    sprintf(input_parameters[0].kernel_name , "%s", argv[i]); i++;
    parameters[0].kernel_width = strtod(argv[i], NULL); i++;
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

#ifdef MT
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

#ifdef HASHING
void GetBox(int* bi, int* bj, double x, double y, param parameters) {
    bi[0] = floor( (x+parameters.Lx/2) / parameters.box_size);
    bj[0] = floor( (y+parameters.Ly/2) / parameters.box_size);
}

void ConstructBoxes(param* parameters, int*** boxes) {
    // CHECK THAT Lx, Ly IS AN INTEGER MULTIPLE OF box_size
    if ((parameters[0].Lx - parameters[0].NxBox*parameters[0].box_size > EPS) || 
        (parameters[0].Ly - parameters[0].NyBox*parameters[0].box_size > EPS)) {
        printf("Length not an integer multiple of size of box");
        exit(2);
    }
    // Construct and initialize boxes
    boxes[0] = (int**) malloc(parameters[0].NxBox * sizeof(int*));
    for (int i = 0; i < parameters[0].NxBox; i++) {
        boxes[0][i] = (int*) malloc(parameters[0].NyBox * sizeof(int));
        for (int j=0; j<parameters[0].NyBox; j++) {
            boxes[0][i][j] = -1;
        }
    }
}

void ConstructNeighbors(int** neighbors, int N) {
    // set all entries to -1
	neighbors[0] = (int*) malloc((2*N) * sizeof(int));
	for( int i = 0 ; i < 2*N; i++)
		neighbors[0][i] = -1;
}

void AddInBox(int index_particle, int bi, int bj, int*** boxes, int** neighbors) {
    int k = boxes[0][bi][bj]; // Save the previous first particle of the box
    boxes[0][bi][bj] = index_particle; // particle index_particle becomes the new first of box
    neighbors[0][2*index_particle] = -1; // index_part doesn't have preceding particle
    neighbors[0][2*index_particle+1] = k; // k becomes succeeding of index_particle
    if (k != -1) neighbors[0][2*k] = index_particle; // if k is a particle, set its preceding to index_particle 
}

void RemoveFromBox(int index_particle, int bi, int bj, int*** boxes, int** neighbors) {
    int next = neighbors[0][2*index_particle+1]; // Store the particle after index_particle
    int prev;
    // check if index_particle is the first one of the box
    if (boxes[0][bi][bj] == index_particle) {
        // The first particle of the box becomes the one after idx_part
        boxes[0][bi][bj] = next;
        neighbors[0][2*next] = -1;
    }
    else {
        // get the previous one
        prev = neighbors[0][2*index_particle];
        // The next of the previous becomes the current next
        neighbors[0][2*prev+1] = next;
        // If next is a particle, its prev becomes the current prev
        if (next != -1) neighbors[0][2*next] = prev;
    }
}

void FreeBoxes(param* parameters, int*** boxes) {
    for (int i = 0; i < parameters[0].NxBox; i++)
        free(boxes[0][i]);
    free(boxes[0]);
}

void FreeNeighbors(int** neighbors) {
    free(neighbors[0]);
}

#ifdef QSAP
void ConstructNeighboringBoxes(param parameters, box*** neighboring_boxes) {

}

double ComputeForces_withSpatialHashing() {

}

double ComputeDensities_withSpatialHashing() {

}

void FreeNeighboringBoxes() {
    
}
#endif

#endif

void InitialConditions(particle* particles, param parameters, int*** boxes, int** neighbors){
    /*
    This function initializes the positions and angles of each particle 
    to a uniformly random number
    */

    // not modifying particles, just where they point to, so don't need pointers
    int i;
    for (i=0; i<parameters.N; i++){
        particles[i].x = parameters.Lx*(-.5+genrand64_real3()); // From -Lx/2 to Lx/2
        particles[i].y = parameters.Ly*(-.5+genrand64_real3()); // From -Ly/2 to Ly/2
        particles[i].theta = 2*M_PI*genrand64_real3(); //M_2_PI calls 2*Pi in C
        #ifdef HASHING
        int bi, bj;
        GetBox(&bi, &bj, particles[i].x, particles[i].y, parameters);
        AddInBox(i, bi, bj, boxes, neighbors);
        #endif
    }
}

void InitialConditionsOrigin(particle* particles, param parameters, int*** boxes, int** neighbors){
    /*
    This function initializes the positions and angles of each particle 
    to a uniformly random number
    */

    // not modifying particles, just where they point to, so don't need pointers
    int i;
    for (i=0; i<parameters.N; i++){
        particles[i].x = 0;
        particles[i].y = 0;
        particles[i].theta = 2*M_PI*genrand64_real3();
        #ifdef HASHING
        int bi, bj;
        GetBox(&bi, &bj, particles[i].x, particles[i].y, parameters);
        AddInBox(i, bi, bj, boxes, neighbors);
        #endif

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
    int j;
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

#ifdef QSAP
double DensityDependentSpeed(double rho, double rho_m, double v_min, double v_max){
    /*
    Give the density-dependent speed (formula in the MIPS review paper)
    */
    return v_max + (v_min - v_max)/2 * (1+tanh(2*rho/rho_m-2));
}

double QuorumSensingSpeed(int i, particle* particles, param parameters){
    /*
    This function calculates the density-dependent speed of particle i
    */

    // Calculate rho at the position of particle i
    double rho = Density(particles[i].x, particles[i].y, particles, parameters);

    // Calculate v(rho)
    double v = DensityDependentSpeed(rho, parameters.rho_m, parameters.v_min, parameters.v_max);

    return v;
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

void UpdateParticles(particle* particles, param parameters, int*** boxes, int** neighbors){   
    /*
    This function updates the positions of particles
    */ 
    int i;
    for (i=0;i<parameters.N;i++){
        #ifdef QSAP
        double v = QuorumSensingSpeed(i, particles, parameters); 
        #elif defined NONE
        double v = parameters.v;
        #elif defined POSITION_DEPENDENT_SPEED
        double v = PositionDependentSpeed(particles[i].x, particles[i].y);
        #endif
        #ifdef HASHING
        int old_bi, old_bj, new_bi, new_bj;
        GetBox(&old_bi, &old_bj, particles[i].x, particles[i].y, parameters);
        #endif
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

        #ifdef HASHING
        GetBox(&new_bi, &new_bj, particles[i].x, particles[i].y, parameters);
        
        if ((new_bi != old_bi) || (new_bj != old_bj)) {
            RemoveFromBox(i, old_bi, old_bj, boxes, neighbors);
            AddInBox(i, new_bi, new_bj, boxes, neighbors);
        }
        #endif
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

void StorePositions(double t, param *parameters, particle* particles, int*** boxes, int** neighbors){
    /*
    This function store the positions and angles of each particles 
    every store_time_interval starting from the inputted next_store_time
    */
    int i;
    if (t >= parameters[0].next_store_time){
        #ifdef OUTPUT_DENSITY
        double x2 = 0, y2 = 0;
        MeanSquareDisplacement(parameters[0].N, particles, &x2, &y2);
        fprintf(parameters[0].data_file,"%lg \t %lg \t", x2, y2); // First two entries will be x2 and y2
        #endif

        for (i=0;i<parameters[0].N;i++){
            fprintf(parameters[0].data_file,"%lg \t %d \t %lg \t %lg \t %lg \t",
                    t,i,particles[i].x,particles[i].y,particles[i].theta);
        }
        fprintf(parameters[0].data_file, "\n");

        #ifdef OUTPUT_DENSITY
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
        #endif

        #ifdef TESTING
        for (int i = parameters[0].NyBox-1; i >= 0; i--) {
            for (int j=0; j<parameters[0].NxBox;j++) {
                fprintf(parameters[0].boxes_file, "%d \t", boxes[0][j][i]);
            }
            fprintf(parameters[0].boxes_file, "\n");
        }
        for (int i = 0; i < parameters[0].N*2; i++) {
            fprintf(parameters[0].boxes_file, "%d \t", neighbors[0][i]);
        }
        fprintf(parameters[0].boxes_file, "\n");
        #endif
    
        // Print out the density profile to a file
        parameters[0].next_store_time += parameters[0].store_time_interval;
    }
    
    fflush(parameters[0].data_file);
}