#include "../include/vmc.hpp"
using namespace arma;
using namespace std;


double local_energy(mat& position, double alpha){
    double energy = (0.5 - 2*alpha*alpha)*accu(position%position) + position.n_rows*position.n_cols*alpha;
    //accu(A%A) is elementwise multiplication and sum over all elements

    return energy;
}

double probability_ratio(mat& old_position, mat& new_position, double alpha){
    double exponent = accu(old_position%old_position) - accu(new_position%new_position);
    //accu(A%A) is elementwise multiplication and sum over all elements

    return exp(2*alpha*exponent);
}

// quantum force - untested
vec quantum_force(mat& position, int index, double alpha){

    // expression: 2* grad psi_T / psi_T
    // for non-interactive case, simplifies to -4*alpha*position
    // set force for particle index
    vec force = -4*alpha*position.col(index);
    return force;
}

double greens_ratio(vec& qf_old, vec& qf_new, mat& position_old, mat& position_new, double time_step, double D, int index){
    // print dimensions of matrix
    vec tmp = 0.5*(0.5*D*time_step*(qf_old - qf_new) + position_old.col(index) - position_new.col(index));
    double exponent = accu(tmp%(qf_old + qf_new)); //elementwise multiplication
    
    return exp(exponent);
}


void monte_carlo(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step){
    
    //prepare output file
    string filename;
    if (importance_sampling){
        filename = "output/N="+ to_string(N_particles) + "_d=" + to_string(N_dimensions) + "_IS_energy.csv";
    } 
    else{
        filename = "output/N="+ to_string(N_particles) + "_d=" + to_string(N_dimensions) + "_energy.csv";
    }
    ofstream ofile;
    ofile.open(filename);
    ofile << "MC,N,d,alpha,energy,variance" << endl;

    //initialize the random number generator
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator;
    generator.seed(seed);
    
    arma::arma_rng::set_seed_random();
    uniform_real_distribution<double> rng_double(0,1);

    mat position = mat(N_dimensions, N_particles).fill(0.0); //initialize all particles at the origin 
    mat new_position = mat(N_dimensions, N_particles).fill(0.0); //initialize all particles at the origin 
    vec stepping_vector = vec(N_dimensions).fill(0.5); //vector filled with 0.5 to get rand(-0.5, 0.5) = rand(0, 1) - 0.5

    vec alpha_values = arma::linspace(0.1, 1.0, 10);

    double D = 0.5; // diffusion constant for importance sampling
    for (int i = 0; i < alpha_values.size(); i++){
        //loop over alpha values
        double alpha = alpha_values[i];
        double energy = 0.0;
        double energy_squared = 0.0;
        double new_energy = 0.0;
        for (int j = 0; j < MC_cycles; j++){
            //loop over MC cycles
            for (int k = 0; k < N_particles; k++){
                //loop over particles
        
                //mat new_position = position;
                double greens = 1.0;
                if (importance_sampling){
                    //qf.print();
                    //position.print();
                    vec qf = quantum_force(position, k, alpha);
                    vec gaussian_random = vec(N_dimensions).randn();
                    new_position.col(k) = position.col(k) + D * qf * time_step + gaussian_random*sqrt(time_step);

                    vec qf_new = quantum_force(new_position, k, alpha);
                    greens = greens_ratio(qf, qf_new, position, new_position, time_step, D, k);

                    //cout << "greens: " << greens << endl;
                }
                else{
                    vec random_walk = vec(N_dimensions).randu() - stepping_vector;
                    new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
                }

                double ratio = greens * probability_ratio(position, new_position, alpha);

                if (rng_double(generator) <= ratio){
                    //check whether or not to move particle k
                    position.col(k) = new_position.col(k);
                    new_energy = local_energy(position, alpha); //update new energy
                }
                else{
                    new_position.col(k) = position.col(k); //move particle k back to original position
                }
            }
            //double new_energy = local_energy(position, alpha);
            energy += new_energy;
            energy_squared += new_energy*new_energy;
        }
        energy = energy/MC_cycles;
        double energy2 = energy_squared/MC_cycles;
        double energy_variance = energy2 - energy*energy;
        ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha_values[i] << "," << energy << "," << energy_variance << endl;
    }
    ofile.close();
}
