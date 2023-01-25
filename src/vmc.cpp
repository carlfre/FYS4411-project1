#include "../include/vmc.hpp"
using namespace arma;
using namespace std;


double local_energy(double position, double alpha){
    double energy = alpha + position*position*(0.5 - 2*alpha*alpha);
    return energy;
}

double probability_ratio(double old_position, double new_position, double alpha){
    double exponent = 2*alpha*(old_position*old_position - new_position*new_position);
    return min(exp(exponent), 1.0);
}

void monte_carlo(int MC_cycles, double step){
    //initialize the random number generator
    // Get a seed from the system clock
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();

    mt19937 generator;
    generator.seed(seed);
    uniform_real_distribution<double> rng_double(0,1);

    vec alpha_values = arma::linspace(0.1, 1.0, 10);
    vec energy_mean = vec(alpha_values.size());
    vec energy_variance = vec(alpha_values.size());
    for (int i = 0; i < alpha_values.size(); i++){
        double alpha = alpha_values[i];
        double energy = 0.0;
        double energy_squared = 0.0;
        double position = 0.0;
        for (int j = 0; j < MC_cycles; j++){
            double new_position = position + step*rng_double(generator);
            double ratio = probability_ratio(position, new_position, alpha);
            if (rng_double(generator) < ratio){
                position = new_position;
            }
            double new_energy = local_energy(position, alpha);
            energy += new_energy;
            energy_squared += new_energy*new_energy;
        }
        energy_mean[i] = energy/MC_cycles;
        energy_variance[i] = (energy_squared - energy*energy/MC_cycles)/(MC_cycles - 1);
        cout << "alpha = " << alpha_values[i] << " energy = " << energy_mean[i] << " variance = " << energy_variance[i] << endl;
    }

}

