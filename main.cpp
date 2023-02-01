#include "include/vmc.hpp"
using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
    if (argc != 4){
        cout << "Usage: ./main number_of_particles number_of_dimensions order_of_MC_cycles " << endl;
        exit(1);
    }
    double step = 1.0;
    int N_particles = atoi(argv[1]);
    int N_dimensions = atoi(argv[2]);
    int MC_cycles = pow(10, atoi(argv[3]));
    monte_carlo(MC_cycles, step, N_particles, N_dimensions);
    return 0;
}