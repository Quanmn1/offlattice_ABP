#define MT
#define PBC
// #define QSAP_TANH
#define QSAP_EXP
// #define QSAP_ZERO_LINEAR
// #define QSAP_ZERO_SMOOTH
// #define QSAP_ZERO_EXP
// #define QSAP_ZERO_ASYMP
#if defined(QSAP_ZERO_LINEAR) || defined(QSAP_ZERO_SMOOTH) || defined(QSAP_ZERO_EXP) || defined(QSAP_ZERO_ASYMP)
#define QSAP_ZERO
#endif
#define HARMONIC
// #define WCA
// #define NONE
// #define POSITION_DEPENDENT_SPEED
#if defined(QSAP_TANH) || defined(QSAP_EXP) || defined(QSAP_ZERO)
#define QSAP
#endif
#if defined(HARMONIC) || defined(WCA)
#define PFAP
#define STRESS_TENSOR
#endif
#define HASHING
#define DENSITY_HISTOGRAM
// #define TESTING_DENSITY
// #define TESTING_FORCE
// #define INIT_SLAB
#define INIT_CIRCLE
// #define WALL

/*
Exit codes:
Nevermind it's just 1
1: Incorrect input
2: Incompatible parameters
3: Memory allocation failed
4: Error arise in initializing
5: Error arise in simulation
*/
