#include "Swarm.h"
#include "example_functions.h"



int main(int argc, char* argv[]) {
    PSO pso; 

    pso.init(
        0.4, 0.5, 0.6,
        {{-5, 5}, {-6, 6}},
        5000
    );
    pso.run(exampleFunctions::mccormick, 100, "./run/mccromic.txt"); 

    pso.init(
        0.4, 0.5, 0.6,
        {{-5, 5}, {-6, 6}},
        5000
    );
    pso.run(exampleFunctions::test, 50, "./run/test.txt"); 

    return 0; 
}