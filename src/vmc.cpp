#include "../include/vmc.hpp"
#include "../include/interaction.hpp"
#include "../include/no_interaction.hpp"
#include "../include/sampling.hpp"

using namespace arma;
using namespace std;


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

vec monte_carlo_parallelized(int N_walkers, double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius){
    int mc_per_walker = MC_cycles / N_walkers;

    cout << "N_walkers: " << N_walkers << endl;
    vector<vec> results;
    for (int i=0; i<N_walkers; i++){
        results.push_back(vec());
    }

    cout.setstate(ios_base::failbit); //turn off cout
    #pragma omp parallel for
    for (int i=0; i<N_walkers; i++){
        results[i] = monte_carlo(alpha, mc_per_walker, step, N_particles, N_dimensions, importance_sampling, time_step, numerical_double_derivative, ndd_h, interactions, gamma, beta, hard_core_radius);
    }
    cout.clear(); //turn on cout

    // Average results from all walkers
    vec sum_results({0, 0, 0});
    for (int i=0; i<N_walkers; i++){
        sum_results += results[i];
    }
    vec average = sum_results / N_walkers;
    return average;

}

vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius){
    //initialize the random number generator
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    //cout << seed << endl;
    mt19937 generator;
    generator.seed(seed);
    string filename_for_density = "output/density_" + to_string(interactions) + "_" + to_string(N_particles) + ".csv";
    ofstream ofile_density;
    ofile_density.open(filename_for_density);
    arma::arma_rng::set_seed_random();
    uniform_real_distribution<double> rng_double(0,1);
    mat position = mat(N_dimensions, N_particles);
    mat relative_position = mat(N_particles, N_particles).fill(0.0); 
    if (interactions){
        set_initial_state(position, relative_position, N_particles, N_dimensions, hard_core_radius);  
    }
    else{
         position = mat(N_dimensions, N_particles).fill(0.0); //initialize all particles at the origin 
         //not necessary to equilibrate the system as it starts out in the most probable state
    }
    mat new_position = position;
    mat relative_position_new = relative_position;
    vec stepping_vector = vec(N_dimensions).fill(0.5); //vector filled with 0.5 to get rand(-0.5, 0.5) = rand(0, 1) - 0.5
    double D = 0.5; // diffusion constant for importance sampling
    // to equilibrate the system
    int accepted_moves = 0; //number of accepted moves
    double fraction = 0.0;
    if (interactions){ // FIXME: case without importance sampling. (infinite loop currently)
        //not necessary to equilibrate the system for the non-interactive case as it starts out in the most probable state
        int N_equilibration = 1000;
        while(0.8 > fraction | fraction > 0.95){
            accepted_moves = 0;
            for (int j = 0; j < N_equilibration; j++){
                for (int k = 0; k < N_particles; k++){
                    double random_number = rng_double(generator);
                    sampling(position, new_position, alpha, k, step, time_step, D, importance_sampling, accepted_moves, stepping_vector, random_number, interactions, relative_position, relative_position_new, gamma, beta, hard_core_radius);
                   }
            }
            fraction = (accepted_moves + 0.0)/(N_particles*N_equilibration);
            if (importance_sampling & fraction < 0.8){
                //if the acceptance rate is too low, decrease the step size
                time_step *= 0.8;
            }
            else if (importance_sampling & fraction > 0.95){
                //if the acceptance rate is too high, increase the step size
                time_step *= 1.2;
            }
        }
        //relative_position.print();
        cout << "Equilibration done" << endl;
        cout << "Fraction of accepted moves: " << (accepted_moves + 0.0)/(N_particles*N_equilibration) << endl;    
        cout << "Final time step: " << time_step << endl;
    }
    double energy = 0.0; //accumulator for the energy
    double energy_squared = 0.0; //accumulator for the energy squared
    double new_energy = 0.0; //storing the local energy
    //values for the derivative of the local energy with respect to alpha
    double der_wf = 0.0; //acculmulator for the derivative of the wave function
    double new_der_wf = 0.0; //storing the current derivative of wf
    
    double der_wf_energy = 0.0;
    accepted_moves = 0; //reset the number of accepted moves after equilibration
    for (int j = 0; j < MC_cycles; j++){
        //loop over MC cycles
        for (int k = 0; k < N_particles; k++){
            //loop over particles and sample new positions
            double random_number = rng_double(generator);
            sampling(position, new_position, alpha, k, step, time_step, D, importance_sampling, accepted_moves, stepping_vector, random_number, interactions, relative_position, relative_position_new, gamma, beta, hard_core_radius);
            
        }
        ofile_density << position.t();
    
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
            //cout << alpha << "," << new_energy << endl;

        }
        new_der_wf = dwf_dalpha(new_position); //update new derivative of wave function 
        //update values for the energy and the derivative of the wave function
        energy += new_energy;
        energy_squared += new_energy*new_energy;
        der_wf += new_der_wf;
        der_wf_energy += new_der_wf*new_energy;
    }
    ofile_density.close();
    energy = energy/(MC_cycles); 
    double energy2 = energy_squared/(MC_cycles);
    double energy_variance = energy2 - energy*energy;

    der_wf = der_wf/(MC_cycles);
    der_wf_energy = der_wf_energy/(MC_cycles);
    double alpha_derivative = 2*(der_wf_energy - der_wf*energy);
    vec result = {energy, energy_variance, alpha_derivative};
    cout << "fraction of accepted moves: " << (accepted_moves + 0.0)/(N_particles*MC_cycles) << endl;
    // string filename = "./../output/density/particles_" + to_string(N_particles) + "_" + to_string(alpha) + ".bin";
    return result;
}



// vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool interactions, double gamma, double beta, double hard_core_radius){

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter, bool interactions, double gamma, double beta, double hard_core_radius){
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
    
    vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, interactions, gamma, beta, hard_core_radius);
    double energy = result[0];
    double energy_variance = result[1];
    vector<double> alpha_values;
    alpha_values.push_back(alpha);
    ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << energy << "," << energy_variance << endl;

    double momentum = 0.3;
    int iter = 1;
    alpha_values.push_back(alpha_values[0] - learning_rate*result[2]);
    while (iter < max_iter and abs(alpha_values[iter] - alpha_values[iter - 1]) > 0.00001){ 
        result = monte_carlo(alpha_values[iter], MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, interactions, gamma, beta, hard_core_radius);
        //old_alpha = new_alpha;
        alpha_values.push_back(alpha_values[iter] - learning_rate*result[2] + momentum*(alpha_values[iter] - alpha_values[iter-1]));
        energy = result[0]; //update new energy
        energy_variance = result[1]; //update energy variance

        ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha_values[iter] << "," << energy << "," << energy_variance << endl;
        iter ++;
    }

    ofile.close();
    if (iter == max_iter){
        cout << "Did not converge after " << iter << " iterations. The final parameters are:" << endl;
        cout << "alpha: " << alpha_values.back() << endl;
        cout << "energy: " << energy << endl;
        cout << "variance: " << energy_variance << endl;
    }
    else{
        cout << "Converged after " << iter << " iterations. The final parameters are:" << endl;
        cout << "alpha: " << alpha_values.back() << endl;
        cout << "energy: " << energy << endl;
        cout << "variance: " << energy_variance << endl;
    }
}

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter){
    // no interactions
    minimize_parameters(MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, learning_rate, max_iter, false, -2.0, -2.0, -2.0);
}