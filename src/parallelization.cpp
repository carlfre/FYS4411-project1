#include "../include/parallelization.hpp"

using namespace arma;
using namespace std;


void merge_files(string filename, int n_walkers){
    string line;
    if (filename != ""){
        // Merge to one file
        ofstream file;
        file.open("output/" + filename + ".csv");
        for (int i=0; i<n_walkers; i++){
            ifstream walker_file;
            walker_file.open("output/" + filename + "_" + to_string(i) + ".csv");
            while (getline(walker_file, line)){
                file << line << endl;
            }
            walker_file.close();
        }
        file.close();
        // Remove individual files
        for (int i=0; i<n_walkers; i++){
            string walker_filename = "output/" + filename + "_" + to_string(i) + ".csv";
            remove(walker_filename.c_str());
        }
    }
    
}

vec parallelized_mcmc(
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
    bool interactions,
    int MC_cycles,
    int n_walkers,
    string density_filename,
    string energy_filename
){
    vector<string> density_filenames = {};
    vector<string> energy_filenames = {};
    vector<VMCWalker> walkers = {};
    vector<vec> results = {};

    int cycles_per_walker = MC_cycles / n_walkers;

    for (int i = 0; i < n_walkers; i++){
        string density_filename_i;
        if (density_filename == ""){
            density_filename_i = "";
        }
        else{
            density_filename_i = density_filename + "_" + to_string(i);
        }

        string energy_filename_i;
        if (energy_filename == ""){
            energy_filename_i = "";
        }
        else{
            energy_filename_i = energy_filename + "_" + to_string(i);
        }

        VMCWalker walker(
            alpha, 
            beta,
            gamma, 
            step, // for no importance sampling
            time_step, // for importance sampling
            hard_core_radius, // for interactions
            ndd_h, // h for finite difference double derivative
            N_particles, 
            N_dimensions, 
            importance_sampling, 
            numerical_double_derivative, 
            interactions);
        walkers.push_back(walker);
        density_filenames.push_back(density_filename_i);
        energy_filenames.push_back(energy_filename_i);
        results.push_back(vec());
    }

    
    // cout.setstate(ios_base::failbit); // Turn off cout

    #pragma omp parallel for
    for (int i = 0; i < n_walkers; i++){
        results[i] = walkers[i].walk(cycles_per_walker, density_filenames[i], energy_filenames[i]);
    }

    // cout.clear(); // Turn on cout

    merge_files(density_filename, n_walkers);
    merge_files(energy_filename, n_walkers);

    // Average results
    vec average_results = zeros<vec>(results[0].n_elem);
    for (int i = 0; i < n_walkers; i++){
        average_results += results[i];
    }
    average_results /= n_walkers;
    return average_results;
}

vec parallelized_mcmc(
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
    bool interactions,
    int MC_cycles,
    int n_walkers)
{
    return parallelized_mcmc(
        alpha, 
        beta,
        gamma, 
        step, // for no importance sampling
        time_step, // for importance sampling
        hard_core_radius, // for interactions
        ndd_h, // h for finite difference double derivative
        N_particles, 
        N_dimensions, 
        importance_sampling, 
        numerical_double_derivative, 
        interactions,
        MC_cycles,
        n_walkers,
        "",
        ""
    );

}