#include "../include/vmc.hpp"
#include "../include/interaction.hpp"

using namespace arma;
using namespace std;


double local_energy(mat& position, double alpha){
    //k is the index of the particle we are calculating the energy for
    double energy = (0.5 - 2*alpha*alpha)
                        *accu(position % position) 
                        + position.n_cols*position.n_rows*alpha;
    //accu(A%A) is elementwise multiplication and sum over all elements
    return energy;
}

double probability_ratio(mat& old_position, mat& new_position, double alpha, int particle_index){
    double exponent = dot(old_position.col(particle_index), old_position.col(particle_index)) 
                    - dot(new_position.col(particle_index), new_position.col(particle_index));

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
    //loop over all particles
    for (int i = 0; i < position.n_cols; i++){
        //loop over all dimensions
        for (int j = 0; j < position.n_rows; j++){
            energy += -0.5*numerical_double_derivative(position, alpha, h, j, i);
        }
    }
    energy /= trial_wavefunction(position, alpha);
    //potential energy
    energy += 0.5*accu(position % position);
    return energy;
}

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
    //option for non-numerical double derivative and no interactions
    return monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, false, -1.0, false, -1.0, -1.0, -1.0);
}

vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h){
    //option for numerical double derivative and no interactions
    return monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, true, ndd_h, false, -1.0, -1.0, -1.0);
}

vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool interactions, double gamma, double beta, double hard_core_radius){
    //option for interactions and no numerical double derivative
    return monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, false, -1.0, interactions, gamma, beta, hard_core_radius);
}

void sampling(mat& position, mat& new_position, double alpha, int k, double step, double time_step, double D, bool importance_sampling, int& accepted_moves, vec& stepping_vector, double random_number, bool interactions, mat& relative_position, mat& relative_position_new, double gamma, double beta, double hard_core_radius){
    double greens = 1.0;
    double ratio = 1.0;
    int N_dimensions = position.n_rows;
    int N_particles = position.n_cols;
    if (importance_sampling){
        vec qf = quantum_force(position, k, alpha);
        vec gaussian_random = vec(N_dimensions).randn();
        new_position.col(k) = position.col(k) + D * qf * time_step + gaussian_random*sqrt(time_step);

        vec qf_new = quantum_force(new_position, k, alpha);
        greens = greens_ratio(qf, qf_new, position, new_position, time_step, D, k);
        ratio = greens * probability_ratio(position, new_position, alpha, k);
        if (random_number <= ratio){
            accepted_moves += 1;
            //check whether or not to move particle k
            position.col(k) = new_position.col(k);
        }
        else{
            new_position.col(k) = position.col(k); //move particle k back to original position
        }
    }
    else{
        if (interactions){  
            //interaction case
            vec random_walk = vec(N_dimensions).randu() - stepping_vector;
            new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
            double rel_pos = 0.0;
            for (int i = 0; i < N_particles; i++ ){
                if (i != k){
                    rel_pos = norm(new_position.col(k) - new_position.col(i));
                    if (rel_pos < hard_core_radius){ //if particles are too close, reject move
                        new_position.col(k) = position.col(k); //move particle k back to original position
                        relative_position_new = relative_position;
                        return; //break out of function
                    }
                    else{
                        relative_position_new(max(i, k), min(i ,k)) = rel_pos;
                    }
                }
            }
            ratio = probability_ratio_naive(position, position, relative_position, relative_position_new, alpha, beta, hard_core_radius, k);
            if (random_number <= ratio){
                accepted_moves += 1;
                //check whether or not to move particle k
                position.col(k) = new_position.col(k);
                relative_position = relative_position_new;
            }
            else{
                new_position.col(k) = position.col(k); //move particle k back to original position
                relative_position_new = relative_position;
            }

        }
        else{
            vec random_walk = vec(N_dimensions).randu() - stepping_vector;
            new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
            ratio = probability_ratio(position, new_position, alpha, k);
            if (random_number <= ratio){
                accepted_moves += 1;
                //check whether or not to move particle k
                position.col(k) = new_position.col(k);
            }
            else{
                new_position.col(k) = position.col(k); //move particle k back to original position
            }
        }
    }
}
vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius){
    //initialize the random number generator
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator;
    generator.seed(seed);
    
    arma::arma_rng::set_seed_random();
    uniform_real_distribution<double> rng_double(0,1);
    mat position;
    mat relative_position = mat(N_particles, N_particles).fill(0.0); 
    if (interactions){
        position = 2*hard_core_radius*(mat(N_dimensions, N_particles).randu() - 0.5); //initialize all particles at random positions
        //position = mat(N_dimensions, N_particles).randu() - 0.5; 
        relative_position = trimatl(relative_position, 1); //strict lower triangular matrix to store relative positions
        for (int i = 0; i < relative_position.n_cols; i++){
            for (int j = i + 1; j < relative_position.n_rows; j++){
                relative_position(j, i) = norm(position.col(i) - position.col(j));
            }
        }
    }
    else{
         position = mat(N_dimensions, N_particles).fill(0.0); //initialize all particles at the origin 
    }
    //relative_position.print();
    
    mat new_position = position; //initialize all particles at the origin 
    mat relative_position_new = relative_position;
    vec stepping_vector = vec(N_dimensions).fill(0.5); //vector filled with 0.5 to get rand(-0.5, 0.5) = rand(0, 1) - 0.5
    double D = 0.5; // diffusion constant for importance sampling
    int accepted_moves = 0;
    // to equilibrate the system
    for (int j = 0; j < 1000000; j++){
        for (int k = 0; k < N_particles; k++){
            double random_number = rng_double(generator);
            sampling(position, new_position, alpha, k, step, time_step, D, importance_sampling, accepted_moves, stepping_vector, random_number, interactions, relative_position, relative_position_new, gamma, beta, hard_core_radius);
        }
    }
    //relative_position.print();
    cout << "Equilibration done" << endl;
    cout << "Accepted moves: " << accepted_moves << endl;
    
    
    double energy = 0.0; //accumulator for the energy
    double energy_squared = 0.0; //accumulator for the energy squared
    double new_energy = 0.0; //storing the local energy
    //values for the derivative of the local energy with respect to alpha
    double der_wf = 0.0; //acculmulator for the derivative of the wave function
    double new_der_wf = 0.0; //storing the current derivative of wf
    
    double der_wf_energy = 0.0;
    accepted_moves = 0;
    for (int j = 0; j < MC_cycles; j++){
        //loop over MC cycles
        for (int k = 0; k < N_particles; k++){
            //loop over particles and sample new positions
            double random_number = rng_double(generator);
            sampling(position, new_position, alpha, k, step, time_step, D, importance_sampling, accepted_moves, stepping_vector, random_number, interactions, relative_position, relative_position_new, gamma, beta, hard_core_radius);
            
        }
        if (numerical_double_derivative && !interactions){
            new_energy = local_energy_numerical(position, alpha, ndd_h); //update new energy
        }
        else if(!interactions && !numerical_double_derivative){
            new_energy = local_energy(position, alpha);
        }
        else if(interactions){
            //new_energy = local_energy(position, alpha);
            //new_energy = local_energy_interactions(position, relative_position, alpha, beta, gamma, hard_core_radius);
            new_energy = local_energy_naive(position, relative_position, alpha, beta, gamma, hard_core_radius);
        }
        new_der_wf = dwf_dalpha(new_position); //update new derivative of wave function 
        //update values for the energy and the derivative of the wave function
        energy += new_energy;
        energy_squared += new_energy*new_energy;
        der_wf += new_der_wf;
        der_wf_energy += new_der_wf*new_energy;
    }
    energy = energy/(MC_cycles); 
    double energy2 = energy_squared/(MC_cycles);
    double energy_variance = energy2 - energy*energy;

    der_wf = der_wf/(MC_cycles);
    der_wf_energy = der_wf_energy/(MC_cycles);
    double alpha_derivative = 2*(der_wf_energy - der_wf*energy);
    vec result = {energy, energy_variance, alpha_derivative};
    cout << "fraction of accepted moves: " << (accepted_moves + 0.0)/(N_particles*MC_cycles) << endl;
    return result;
}

double dwf_dalpha(arma::mat& position){
    return -accu(position % position);
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
