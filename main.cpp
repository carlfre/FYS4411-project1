#include "include/vmc.hpp"
using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
    if (argc != 6){
        cout << "Usage: ./main number_of_particles number_of_dimensions order_of_MC_cycles important_sampling time_step" << endl;
        exit(1);
    }
    double step = 1.0;
    int N_particles = atoi(argv[1]);
    int N_dimensions = atoi(argv[2]);
    int MC_cycles = pow(10, atoi(argv[3]));
    double time_step = atof(argv[5]);
    bool importance_sampling;
    if (atoi(argv[4]) == 0){
        importance_sampling = false;
        cout << "No important sampling" << endl;
    }
    else if (atoi(argv[4]) == 1){
        importance_sampling = true;
        cout << "Important sampling" << endl;
    }
    monte_carlo(MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step);
    return 0;
}