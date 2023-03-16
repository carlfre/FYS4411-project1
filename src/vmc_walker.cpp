#include "../include/vmc_walker.hpp"
#include "../include/interaction.hpp"
#include "../include/no_interaction.hpp"

using namespace arma;
using namespace std;

VMCWalker::VMCWalker(
            double alpha, 
            double beta,
            double gamma, 
            double step, // for no importance sampling
            double time_step, // for importance sampling
            double hard_core_radius, // for interactions
            double ndd_h, // h for finite difference double derivative
            int N_particles, 
            int N_dimensions, 
            bool importance_sampling, 
            bool numerical_double_derivative, 
            bool interactions
        )
{
    // set values of class variables
    this->alpha = alpha;
    this->beta = beta;
    this->gamma = gamma;
    this->step = step;
    this->time_step = time_step;
    this->D = 0.5;
    this->hard_core_radius = hard_core_radius;
    this->N_particles = N_particles;
    this->N_dimensions = N_dimensions;
    this->importance_sampling = importance_sampling;
    this->numerical_double_derivative = numerical_double_derivative;
    this->interactions = interactions;
    this->adjust_step_automatically=true;
    this->ndd_h = ndd_h;

    // initialize random number generator 
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    generator = mt19937_64(seed);
    rng_double = uniform_real_distribution<double>(0.0, 1.0);

    // initialize stepping vector (for no importance sampling)
    // this is a vector filled with 0.5 to get rand(-0.5, 0.5) = rand(0, 1) - 0.5
    this->stepping_vector = vec(N_dimensions).fill(0.5); 
}

void VMCWalker::set_step_adjustment(bool use_adjustment){
    adjust_step_automatically=use_adjustment;
}

vec VMCWalker::walk(int MC_cycles, string density_filename, string energy_filename)
{
    // if density_filename is specified, open a file for tracking positions
    ofstream ofile_density, ofile_energy;
    if (density_filename != ""){
        ofile_density.open("output/" + density_filename + ".csv"); 
    }
    // if energy_filename is specified, open a file for tracking energy
    if (energy_filename != ""){
        ofile_energy.open("output/" + energy_filename + ".csv"); 
    }

    accepted_moves = 0; 

    // Initialize particle positions, and burn in
    initialize();

    double energy = 0.0; //accumulator for the energy
    double energy_squared = 0.0; //accumulator for the energy squared
    double new_energy = 0.0; //storing the local energy
    //values for the derivative of the local energy with respect to alpha
    double der_wf = 0.0; //acculmulator for the derivative of the wave function
    double new_der_wf = 0.0; //storing the current derivative of wf
    double der_wf_energy = 0.0;

    accepted_moves = 0; //reset the number of accepted moves after equilibration

    // MCMC loop
    for (int j = 0; j < MC_cycles; j++){
        //loop over MC cycles
        for (int k = 0; k < N_particles; k++){
            sampling(k);            
        }
    
        if (numerical_double_derivative && !interactions){
            new_energy = local_energy_numerical(position, alpha, ndd_h); //update new energy
        }
        else if(!interactions && !numerical_double_derivative){
            new_energy = local_energy(position, alpha);
        }
        else if(interactions && !numerical_double_derivative){
            new_energy = local_energy_naive(position, relative_position, alpha, beta, gamma, hard_core_radius); // TODO: Switch to other implementation?
        }
        else{
            throw "Numerical double derivative + interactions is not implemented";
        }
        
        //update values for the energy and the derivative of the wave function
        new_der_wf = dwf_dalpha(new_position); //update new derivative of wave function
        energy += new_energy;
        energy_squared += new_energy*new_energy;

        der_wf += new_der_wf;
        der_wf_energy += new_der_wf*new_energy;

        // write all particle positions to file
        if (density_filename != ""){
            ofile_density << new_position.t();
        }
        // write energy to file
        if (energy_filename != ""){
            ofile_energy << new_energy << endl;
        }
    }

    ofile_density.close();
    ofile_energy.close();

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

vec VMCWalker::walk(int MC_cycles){
    return walk(MC_cycles, "", "");
}

void VMCWalker::minimize_parameters(int MC_cycles, double learning_rate, int max_iter){
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

    vec result = walk(MC_cycles);
    double energy = result[0];
    double energy_variance = result[1];
    vector<double> alpha_values;
    alpha_values.push_back(alpha);
    ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << energy << "," << energy_variance << endl;

    double momentum = 0.3;
    int iter = 1;
    alpha = alpha_values[0] - learning_rate*result[2];
    alpha_values.push_back(alpha);
    while (iter < max_iter and abs(alpha_values[iter] - alpha_values[iter - 1]) > 0.00001){ 
        result = walk(MC_cycles);
        alpha = alpha_values[iter] - learning_rate*result[2] + momentum*(alpha_values[iter] - alpha_values[iter-1]);
        alpha_values.push_back(alpha);
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