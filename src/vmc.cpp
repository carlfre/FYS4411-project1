#include "../include/vmc.hpp"
using namespace arma;
using namespace std;


double local_energy(mat position, double alpha){
    double energy = 0;  
    for (int i = 0; i < position.n_rows; i++){ //loop over particles
        energy += arma::dot(position.row(i), position.row(i));
    }
    energy *= (0.5 - 2*alpha*alpha);
    energy += position.n_rows*position.n_cols*alpha;
    return energy;
}

double probability_ratio(mat old_position, mat new_position, double alpha){
    double exponent = 0.0;
    for (int i = 0; i < old_position.n_rows; i++){ //loop over particles
        exponent += (arma::dot(old_position.row(i), old_position.row(i)) - arma::dot(new_position.row(i), new_position.row(i)));
    }
    return min(exp(2*alpha*exponent), 1.0);
}

void monte_carlo(int MC_cycles, double step, int N_particles, int N_dimensions){
    //initialize the random number generator
    // Get a seed from the system clock
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();

    mt19937 generator;
    generator.seed(seed);
    arma_rng::set_seed_random();
    //arma::arma_rng::set_seed_random();
    uniform_real_distribution<double> rng_double(0,1);
    mat position = mat(N_particles, N_dimensions).fill(0.0);

    vec alpha_values = arma::linspace(0.1, 1.0, 10);
    for (int i = 0; i < alpha_values.size(); i++){
        double alpha = alpha_values[i];
        double energy = 0.0;
        double energy_squared = 0.0;
        for (int j = 0; j < MC_cycles; j++){
            mat random_walk = mat(N_particles, N_dimensions).randu() - mat(N_particles, N_dimensions).fill(0.5);
            mat new_position = position + step*random_walk;
            double ratio = probability_ratio(position, new_position, alpha);
            if (rng_double(generator) <= ratio){
                position = new_position;
            }
            double new_energy = local_energy(position, alpha);
            energy += new_energy;
            energy_squared += new_energy*new_energy;
        }
        energy = energy/MC_cycles;
        double energy2 = energy_squared/MC_cycles;
        double energy_variance = energy2 - energy*energy;
        cout << "alpha = " << alpha_values[i] << " energy = " << energy << " variance = " << energy_variance << endl;
    }
}

