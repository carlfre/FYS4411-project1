#include "../include/vmc.hpp"
using namespace arma;
using namespace std;


double local_energy(mat position, double alpha){
    double energy = (0.5 - 2*alpha*alpha)*accu(position%position) + position.n_rows*position.n_cols*alpha;
    /*  
    for (int i = 0; i < position.n_rows; i++){ //loop over particles
        energy += arma::dot(position.row(i), position.row(i));
    }
    
    energy *= (0.5 - 2*alpha*alpha);
    energy += position.n_rows*position.n_cols*alpha;
    */
    return energy;
}

double probability_ratio(mat old_position, mat new_position, double alpha){
    double exponent = accu(old_position%old_position) - accu(new_position%new_position);
    /*
    for (int i = 0; i < old_position.n_rows; i++){ //loop over particles
        exponent += (arma::dot(old_position.row(i), old_position.row(i)) - arma::dot(new_position.row(i), new_position.row(i)));
    }
    */
    return min(exp(2*alpha*exponent), 1.0);
}

// quantum force - untested
mat quantum_force(mat& position, int index, double alpha){

    // expression: 2* grad psi_T / psi_T
    // for non-interactive case, simplifies to -4*alpha*position

    // create matrix with same dimensions as position, but with all zeros
    mat force;
    force.zeros(position.n_rows, position.n_cols);

    // set force for particle index
    force.row(index) = -4*alpha*position.row(index);
    return force;
}

double greens_ratio(mat& qf_old, mat& qf_new, mat& position_old, mat& position_new, double step, double D, int index){
    // print dimensions of matrix
    double exponent = 0.0;
    for (int j=0; j<qf_old.n_cols; j++){
        double term = (0.5
                        * (qf_old(index, j) + qf_new(index, j))
                        * (D * step * 0.5 * (qf_old(index, j) - qf_new(index, j)) 
                            + position_old(index, j) - position_new(index, j)));
        exponent += term;
    }
    return exp(exponent);
}


void monte_carlo(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling){
    
    //prepare output file
    string filename = "output/N="+ to_string(N_particles) + "_d=" + to_string(N_dimensions) + "_energy.csv";
    ofstream ofile;
    ofile.open(filename);
    ofile << "N,d,alpha,energy,variance" << endl;

    //initialize the random number generator
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator;
    generator.seed(seed);
    
    arma::arma_rng::set_seed_random();
    uniform_real_distribution<double> rng_double(0,1);

    mat position = mat(N_particles, N_dimensions).fill(0.0); //initialize all particles at the origin
    vec stepping_vector = vec(N_dimensions).fill(0.5); //vector filled with 0.5 to get rand(-0.5, 0.5) = rand(0, 1) - 0.5

    vec alpha_values = arma::linspace(0.1, 1.0, 10);

    double D = 0.5; // diffusion constant for importance sampling
    for (int i = 0; i < alpha_values.size(); i++){
        //loop over alpha values
        double alpha = alpha_values[i];
        double energy = 0.0;
        double energy_squared = 0.0;
        cout << "alpha: " << alpha << endl;
        for (int j = 0; j < MC_cycles; j++){
            //loop over MC cycles
            for (int k = 0; k < N_particles; k++){
                //loop over particles
                vec random_walk = vec(N_dimensions).randu() - stepping_vector;
                mat new_position = position;
                new_position.row(k) += step*random_walk.t(); //move particle k with a random walk

                double greens = 1.0;
                if (importance_sampling){
                    mat qf = quantum_force(position, k, alpha);
                    qf.print();
                    position.print();
                    new_position += D * qf * step;

                    mat qf_new = quantum_force(new_position, k, alpha);
                    greens = greens_ratio(qf, qf_new, position, new_position, step, D, k);

                    cout << "greens: " << greens << endl;
                }

                double ratio = greens * probability_ratio(position, new_position, alpha);

                if (rng_double(generator) <= ratio){
                    //check whether or not to move particle k
                    position = new_position;
                }
            }
            double new_energy = local_energy(position, alpha);
            energy += new_energy;
            energy_squared += new_energy*new_energy;
        }
        energy = energy/MC_cycles;
        double energy2 = energy_squared/MC_cycles;
        double energy_variance = energy2 - energy*energy;
        ofile << N_particles << "," << N_dimensions << "," << alpha_values[i] << "," << energy << "," << energy_variance << endl;
    }
    ofile.close();
}

void monte_carlo(int MC_cycles, double step, int N_particles, int N_dimensions){
    monte_carlo(MC_cycles, step, N_particles, N_dimensions, false);
}