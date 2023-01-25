#include "include/vmc.hpp"

int main(int argc, char *argv[]){
    int MC_cycles = 1000;
    double step = 1.0;
    monte_carlo(MC_cycles, step);
    return 0;
}