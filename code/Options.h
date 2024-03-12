#define MT
#define PBC
// #define QSAP_TANH
#define QSAP_EXP
// #define PFAP
// #define NONE
// #define POSITION_DEPENDENT_SPEED
#if defined(QSAP_TANH) || defined(QSAP_EXP)
#define QSAP
#endif
#define HASHING
#define DENSITY_HISTOGRAM
// #define TESTING_DENSITY
// #define TESTING_FORCE
#define INIT_SLAB

/*
Exit codes:
Nevermine it's just 1
1: Incorrect input
2: Incompatible parameters
3: Memory allocation failed
4: Error arise in initializing
5: Error arise in simulation
*/
