// #include "include/vmc.hpp"
// #include "include/interaction.hpp"
#include "include/vmc_walker.hpp"
#include "include/parallelization.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <chrono>

using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage: taskname " << endl;
        exit(1);
    }
    string task = string(argv[1]);
    string filename = "configs/" + string(argv[1]) + ".txt";
    // A vector of vectors to store the the values in the input file
    vector<double> input_data;
    fstream myfile;
    myfile.open(filename, ios::in);
    if (myfile.is_open())
    { // This checks that the file was opened OK
        string line;
        double param;
        // Read file line by line
        while (getline(myfile, line))
        {
            // Skip lines with "#" at the first position
            if (line.at(0) == '#')
            {
                continue;
            }
            else
            {
                // Parse the string (line) and interpret it as one variable
                stringstream mysstream(line);
                mysstream >> param;
                input_data.push_back(param);
            }
        }
    }
    else
    {
        cout << "Unable to open the file " << filename << endl;
        exit(1);
    }
    myfile.close();

    // TODO: Fix block below
    if (task == "analytical")
    {
        double beta = input_data[0];
        double gamma = input_data[1];
        double step = input_data[2];
        double time_step = input_data[3];
        double hard_core_radius = input_data[4];
        double ndd_h = -1;
        int N_particles = input_data[5];
        int N_dimensions = input_data[6];
        bool importance_sampling = input_data[7];
        bool numerical_double_derivative = false;
        bool interactions = input_data[8];
        int MC_cycles = pow(10, input_data[9]);

        string filename;
        if (importance_sampling)
        {
            importance_sampling = true;
            filename = "output/N=" + to_string(N_particles) +
                       "_d=" + to_string(N_dimensions) + "_ana_IS.csv";
            cout << "Importance sampling" << endl;
        }
        else
        {
            importance_sampling = false;
            cout << "No importance sampling" << endl;
            filename = "output/N=" + to_string(N_particles) +
                       "_d=" + to_string(N_dimensions) + "_ana.csv";
        }

        ofstream ofile;
        ofile.open(filename);
        cout << filename << endl;
        ofile << "MC,N,d,alpha,energy,variance" << endl;
        vec alpha_values = linspace(0.1, 1.0, 10);
        for (double alpha : alpha_values)
        {
            VMCWalker walker(
                alpha,
                beta,
                gamma,
                step,
                time_step,
                hard_core_radius,
                ndd_h,
                N_particles,
                N_dimensions,
                importance_sampling,
                numerical_double_derivative,
                interactions);

            vec result = walker.walk(MC_cycles);
            ofile << MC_cycles << "," << N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
        }
    }
    else if (task == "numerical") 
    {
        double beta = -1;
        double gamma = -1;
        double step = input_data[0];
        double time_step = input_data[1];
        double hard_core_radius = -1;
        double ndd_h = input_data[2]; // h for finite difference double derivative
        int N_particles = input_data[3];
        int N_dimensions = input_data[4];
        bool importance_sampling = input_data[5];
        bool numerical_double_derivative = true;
        bool interactions = false;
        int MC_cycles = pow(10, (int)input_data[6]);

        string filename;
        if (importance_sampling)
        {
            importance_sampling = true;
            filename = "output/numerical_N=" + to_string(N_particles) +
                       "_d=" + to_string(N_dimensions) + "num_IS.csv";
            cout << "Importance sampling" << endl;
        }
        else
        {
            importance_sampling = false;
            cout << "No importance sampling" << endl;
            filename = "output/N=" + to_string(N_particles) +
                       "_d=" + to_string(N_dimensions) + "_num.csv";
        }

        ofstream ofile;
        ofile.open(filename);
        cout << filename << endl;
        ofile << "MC,N,d,alpha,energy,variance" << endl;
        vec alpha_values = linspace(0.1, 1.0, 10);
        for (double alpha : alpha_values)
        {
            VMCWalker walker(
                alpha,
                beta,
                gamma,
                step,
                time_step,
                hard_core_radius,
                ndd_h,
                N_particles,
                N_dimensions,
                importance_sampling,
                numerical_double_derivative,
                interactions);

            vec result = walker.walk(MC_cycles);
            // vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, true, ndd_h);
            ofile << MC_cycles << "," << N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
        }
    }

    else if (task == "density")
    {
        double alpha = input_data[0];
        double beta = input_data[1];
        double gamma = input_data[2];
        double step = input_data[3];
        double time_step = input_data[4];
        double hard_core_radius = input_data[5];
        double ndd_h = -1; // not in use
        int N_particles = input_data[6];
        int N_dimensions = input_data[7];
        bool importance_sampling = input_data[8];
        bool numerical_double_derivative = false; // not in use
        bool interactions = input_data[9];

        // int MC_cycles = input_data[10];
        int MC_cycles = pow(2, input_data[10]);
        int n_walkers = input_data[11];
        string density_filename = "density_N=" + to_string(N_particles); // + "_r=" + to_string(hard_core_radius).substr(0, 4)

        parallelized_mcmc(
            alpha,
            beta,
            gamma,
            step,
            time_step,
            hard_core_radius,
            ndd_h,
            N_particles,
            N_dimensions,
            importance_sampling,
            numerical_double_derivative,
            interactions,
            MC_cycles,
            n_walkers,
            density_filename,
            "");
    }
    else if (task == "statistics")
    {
        double alpha = input_data[0];
        double beta = input_data[1];
        double gamma = input_data[2];
        double step = input_data[3];
        double time_step = input_data[4];
        double hard_core_radius = input_data[5];
        double ndd_h = -1; // not in use
        int N_particles = input_data[6];
        int N_dimensions = input_data[7];
        bool importance_sampling = input_data[8];
        bool numerical_double_derivative = false; // not in use
        bool interactions = input_data[9];
        int MC_cycles = pow(2, input_data[10]); // NOTICE: base 2 here, to work better with blocking!
        int n_walkers = input_data[11];
        string energy_filename = "energy_statistics_N=" + to_string(N_particles) + "_d=" + to_string(N_dimensions);
        if (interactions){
            energy_filename += "_interactions";
        }
        else{
            energy_filename += "_no_interactions";
        }
        cout << "Filename: " << energy_filename << endl;
        cout << "MC cycles: " << MC_cycles << endl;

        parallelized_mcmc(
            alpha,
            beta,
            gamma,
            step,
            time_step,
            hard_core_radius,
            ndd_h,
            N_particles,
            N_dimensions,
            importance_sampling,
            numerical_double_derivative,
            interactions,
            MC_cycles,
            n_walkers,
            "",
            energy_filename);
    }
    else if (task == "gradient")
    {
        double alpha = input_data[0];
        double beta = input_data[1];
        double gamma = input_data[2];
        double step = input_data[3];
        double time_step = input_data[4];
        double hard_core_radius = input_data[5];
        double ndd_h = -1; // not in use
        int N_particles = input_data[6];
        int N_dimensions = input_data[7];
        bool importance_sampling = input_data[8];
        bool numerical_double_derivative = false;
        bool interactions = input_data[9];
        int MC_cycles = pow(10, input_data[10]);
        double learning_rate = input_data[11];
        int max_iterations = input_data[12];

        if (importance_sampling)
        {
            cout << "Importance sampling" << endl;
        }
        else
        {
            cout << "No importance sampling" << endl;
        }

        VMCWalker walker(
            alpha,
            beta,
            gamma,
            step,             // for no importance sampling
            time_step,        // for importance sampling
            hard_core_radius, // for interactions
            ndd_h,            // h for finite difference double derivative
            N_particles,
            N_dimensions,
            importance_sampling,
            numerical_double_derivative,
            interactions);

        walker.minimize_parameters(MC_cycles, learning_rate, max_iterations);
    }

    else if (task == "interactions")
    {
        double beta = input_data[0];
        double gamma = input_data[1];
        double step = input_data[2];
        double time_step = input_data[3];
        double hard_core_radius = input_data[4];
        double ndd_h = -1;
        int N_particles = input_data[5];
        int N_dimensions = input_data[6];
        bool importance_sampling = input_data[7];
        bool numerical_double_derivative = false;
        bool interactions = true;
        int MC_cycles = pow(10, input_data[8]);
        int n_walkers = input_data[9];

        if (importance_sampling)
        {
            filename = "output/N=" + to_string(N_particles) +
                       "_d=" + to_string(N_dimensions) + "_int_IS.csv";
            cout << "Importance sampling" << endl;
        }
        else
        {
            cout << "No importance sampling" << endl;
            filename = "output/N=" + to_string(N_particles) +
                       "_d=" + to_string(N_dimensions) + "_int.csv";
        }
        ofstream ofile;
        ofile.open(filename);
        ofile << "MC,N,d,alpha,energy,variance" << endl;

        vec alpha_values = linspace(0.45, 0.55, 10);
        // vec step_sizes = step * linspace(0.4, 1.0, 10);
        // vec time_steps = {0.7, 0.7, 0.65, 0.6, 0.5, 0.2, 0.1, 0.07, 0.04, 0.01};
        // // vec time_steps = time_step*linspace(0.5, 1.5, 10);
        // double alpha = 0.5; // TODO: what task should this solve? loop over alpha values? step_sizes?
        
        // result.print();
        int N_alpha = alpha_values.size();
        for (int i = 0; i < N_alpha; i++){
            double alpha = alpha_values[i];
            vec result = parallelized_mcmc(
            alpha,
            beta,
            gamma,
            step,
            time_step,
            hard_core_radius,
            ndd_h,
            N_particles,
            N_dimensions,
            importance_sampling,
            numerical_double_derivative,
            interactions,
            MC_cycles,
            n_walkers);
            // vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, interactions, gamma, beta, hard_core_radius);
            ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
        }
    }

      else if(task == "interactions_gradient"){ 

        double alpha_0 = input_data[0];
        double beta = input_data[1];
        double gamma = input_data[2];
        double step = input_data[3];
        double time_step = input_data[4];
        double hard_core_radius = input_data[5];
        double ndd_h = -1;
        int N_particles = input_data[6];
        int N_dimensions = input_data[7];
        bool importance_sampling = input_data[8];
        bool numerical_double_derivative = false;
        bool interactions = true;
        int MC_cycles = pow(10, input_data[9]);
        double learning_rate = input_data[10];
        int max_iterations = input_data[11];


        int N_averages = 32;
        int n_cores = 16;

        vector<string> filenames = {};
        vector<vec> results = {};
        for (int i=0; i<N_averages; i++){
            string filename = "output/N=" + to_string(N_particles) + "_int_grad_" + to_string(i) + ".csv";
            filenames.push_back(filename);
            results.push_back(vec(2));
            
        }

        for (int j=0; j<N_averages; j+=n_cores){
            #pragma omp parallel for
            for(int i=0; i<n_cores; i++){
                string filename = filenames[i+j];
                VMCWalker walker(
                    alpha_0,
                    beta,
                    gamma,
                    step,             // for no importance sampling
                    time_step,        // for importance sampling
                    hard_core_radius, // for interactions
                    ndd_h,            // h for finite difference double derivative
                    N_particles,
                    N_dimensions,
                    importance_sampling,
                    numerical_double_derivative,
                    interactions);
                vec res = walker.minimize_parameters(MC_cycles, learning_rate, max_iterations, filename);
                results[i+j] = res;
            }
        }
        // #pragma omp parallel for
        // for(int i=0; i<N_averages; i++){
        //     string filename = filenames[i];
        //     VMCWalker walker(
        //         alpha_0,
        //         beta,
        //         gamma,
        //         step,             // for no importance sampling
        //         time_step,        // for importance sampling
        //         hard_core_radius, // for interactions
        //         ndd_h,            // h for finite difference double derivative
        //         N_particles,
        //         N_dimensions,
        //         importance_sampling,
        //         numerical_double_derivative,
        //         interactions);
        //     vec res = walker.minimize_parameters(MC_cycles, learning_rate, max_iterations, filename);
        //     results[i] = res;
        // }

        // write to file
        ofstream ofile;
        ofile.open("output/interactions_gradient_N=" + to_string(N_particles) + ".csv");
        ofile << "alpha,iterations" << endl;
        for (int i=0; i<N_averages; i++){
            ofile << results[i][0] << "," << results[i][1] << endl;
        }
      }
    else if (task == "parallel"){ // Worth noting: each walker has individual burn-in period.
        double alpha = input_data[0];
        double beta = -1;
        double gamma = -1;
        double step = input_data[1];
        double time_step = input_data[2];
        double hard_core_radius = -1;
        double ndd_h = -1;
        int N_particles = input_data[3];
        int N_dimensions = input_data[4];
        bool importance_sampling = input_data[5];
        bool numerical_double_derivative = false;
        bool interactions = 0;
        int MC_cycles = pow(10, input_data[6]);
        int n_walkers = input_data[7];
        int n_averages = input_data[8];

        double t1 = 0; // parallel time
        double t2 = 0; // serial time

        for (int i = 0; i < n_averages; i++)
        {

            // time parallelized monte carlo
            auto start1 = chrono::high_resolution_clock::now();
            vec r1 = parallelized_mcmc(
                alpha,
                beta,
                gamma,
                step,
                time_step,
                hard_core_radius,
                ndd_h,
                N_particles,
                N_dimensions,
                importance_sampling,
                numerical_double_derivative,
                interactions,
                MC_cycles,
                n_walkers);
            auto finish1 = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed1 = finish1 - start1;

            // time serial monte carlo
            auto start2 = chrono::high_resolution_clock::now();
            VMCWalker walker(
                alpha,
                beta,
                gamma,
                step,             // for no importance sampling
                time_step,        // for importance sampling
                hard_core_radius, // for interactions
                ndd_h,            // h for finite difference double derivative
                N_particles,
                N_dimensions,
                importance_sampling,
                numerical_double_derivative,
                interactions);
            vec r2 = walker.walk(MC_cycles);
            auto finish2 = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed2 = finish2 - start2;

            // tally up time
            t1 += elapsed1.count();
            t2 += elapsed2.count();

            // Verify that the results are sensible
            if (i == 0)
            {
                cout << "for i=0:" << endl;
                cout << "r1 = " << r1 << endl;
                cout << "r2 = " << r2 << endl;
            }
        }
        cout << "Time parallelized: " << t1 / (double)n_averages << " s" << endl;
        cout << "Time serial: " << t2 / (double)n_averages << " s" << endl;
    }
    else
    {
        cout << "Unknown task" << endl;
    }
    return 0;
}