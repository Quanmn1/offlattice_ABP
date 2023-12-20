#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "math_helper.c"

int main(){
    double a = RadialIntegrate(KernelExp, 0, 1, 1000);
    printf("%f", a);
}
