#include "../include/vmc.hpp"
using namespace arma;
using namespace std;


double local_energy(mat& position, double alpha){
    double energy = (0.5 - 2*alpha*alpha)*accu(position%position) + position.n_rows*position.n_cols*alpha;
    //accu(A%A) is elementwise multiplication and sum over all elements

    return energy;
}

double probability_ratio(mat& old_position, mat& new_position, double alpha, int index){
    double exponent = accu(old_position.col(index)%old_position.col(index)) - accu(new_position.col(index)%new_position.col(index));
    //accu(A%A) is elementwise multiplication and sum over all elements

    return exp(2*alpha*exponent);
}


double trial_wavefunction(mat& position, double alpha){
    double exponent = accu(position%position);
    //accu(A%A) is elementwise multiplication and sum over all elements
    return exp(-alpha*exponent);
}

double numerical_double_derivative(mat& position, double alpha, double h, int dimension_index, int particle_index){
    double wf_mid = trial_wavefunction(position, alpha);
    //move particle in dimension dimension_index h forwards
    position(dimension_index, particle_index) += h;
    double wf_plus = trial_wavefunction(position, alpha);
    //move particle in dimension dimension_index h backwards
    position(dimension_index, particle_index) -= 2*h;
    double wf_minus = trial_wavefunction(position, alpha);
    //move particle back to original position
    position(dimension_index, particle_index) += h;
    double double_derivative = (wf_plus - 2*wf_mid + wf_minus)/(h*h);
    return double_derivative;
}

double local_energy_numerical(mat& position, double alpha, double h){
    double energy = 0;
    for (int i = 0; i < position.n_cols; i++){
        for (int j = 0; j < position.n_rows; j++){
            //kinetic energy
            energy += -0.5*numerical_double_derivative(position, alpha, h, j, i);
            //potential energy
        }
    }
    energy /= trial_wavefunction(position, alpha);
    energy += 0.5*accu(position%position);
    return energy;
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

vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step){
    //default is non-numerical double derivative
    return monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, false, -1.0);
}

vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h){
    //initialize the random number generator
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator;
    generator.seed(seed);
    
    arma::arma_rng::set_seed_random();
    uniform_real_distribution<double> rng_double(0,1);

    mat position = mat(N_dimensions, N_particles).fill(0.0); //initialize all particles at the origin 
    mat new_position = mat(N_dimensions, N_particles).fill(0.0); //initialize all particles at the origin 
    vec stepping_vector = vec(N_dimensions).fill(0.5); //vector filled with 0.5 to get rand(-0.5, 0.5) = rand(0, 1) - 0.5

    double D = 0.5; // diffusion constant for importance sampling
    double energy = 0.0;
    double energy_squared = 0.0;
    double new_energy = 0.0;
    //values for the derivative of the local energy with respect to alpha
    double der_wf = 0.0;
    double der_wf_energy = 0.0;
    double new_der_wf = 0.0;
    int accepted_moves = 0;
    for (int j = 0; j < MC_cycles; j++){
        //loop over MC cycles
        for (int k = 0; k < N_particles; k++){
            //loop over particles
            double greens = 1.0;
            if (importance_sampling){
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

            double ratio = greens * probability_ratio(position, new_position, alpha, k);

            if (rng_double(generator) <= ratio){
                accepted_moves += 1;
                //check whether or not to move particle k
                position.col(k) = new_position.col(k);
                if (numerical_double_derivative){
                    new_energy = local_energy_numerical(position, alpha, ndd_h);
                }
                else{
                    new_energy = local_energy(position, alpha);
                }
                new_der_wf = dwf_dalpha(new_position); //update new derivative of wave function
            }
            else{
                new_position.col(k) = position.col(k); //move particle k back to original position
            }
        }
        //double new_energy = local_energy(position, alpha);
        energy += new_energy;
        energy_squared += new_energy*new_energy;
        der_wf += new_der_wf;
        der_wf_energy += new_der_wf*new_energy;
    }
    energy = energy/MC_cycles;
    double energy2 = energy_squared/MC_cycles;
    double energy_variance = energy2 - energy*energy;

    der_wf = der_wf/MC_cycles;
    der_wf_energy = der_wf_energy/MC_cycles;
    double alpha_derivative = 2*(der_wf_energy - der_wf*energy);
    vec result = {energy, energy_variance, alpha_derivative};
    //cout << "fraction of accepted moves: " << (accepted_moves + 0.0)/(N_particles*MC_cycles) << endl;
    return result;
}

double dwf_dalpha(arma::mat& position){
    return -accu(position%position);
}

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter){
    //prepare output file
    string filename;
    if (importance_sampling){
        filename = "output/N="+ to_string(N_particles) + "_d=" + to_string(N_dimensions) + "_GD_IS_energy.csv";
    } 
    else{
        filename = "output/N="+ to_string(N_particles) + "_d=" + to_string(N_dimensions) + "_GD_energy.csv";
    }
    ofstream ofile;
    ofile.open(filename);
    ofile << "MC,N,d,alpha,energy,variance" << endl;

    double alpha = 0.1; //starting value for alpha

    vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step);
    double old_energy = result[0];
    double energy_variance = result[1];

    ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << old_energy << "," << energy_variance << endl;

    alpha -= learning_rate*result[2];
    
    result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step);
    double new_energy = result[0];
    energy_variance = result[1];

    ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << new_energy << "," << energy_variance << endl;

    alpha -= learning_rate*result[2];

    int iter = 0;
    while (iter < max_iter and abs(old_energy - new_energy) > 0.0001){
        result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step);
        alpha -= learning_rate*result[2];
        old_energy = new_energy; //update old energy
        new_energy = result[0]; //update new energy
        energy_variance = result[1]; //update energy variance

        ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << new_energy << "," << energy_variance << endl;
        iter ++;
    }

    ofile.close();
    if (iter == max_iter){
        cout << "Did not converge after " << iter << " iterations. The final parameters are:" << endl;
        cout << "alpha: " << alpha << endl;
        cout << "energy: " << new_energy << endl;
        cout << "variance: " << energy_variance << endl;
    }
    else{
        cout << "Converged after " << iter << " iterations. The final parameters are:" << endl;
        cout << "alpha: " << alpha << endl;
        cout << "energy: " << new_energy << endl;
        cout << "variance: " << energy_variance << endl;
    }
}
