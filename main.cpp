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
    // if (task == "analytical"){
    //     double beta = input_data[0];
    //     double gamma = input_data[1];
    //     double step = input_data[2];
    //     double time_step = input_data[3];
    //     double hard_core_radius = input_data[4];
    //     double ndd_h = -1;
    //     int N_particles = input_data[5];
    //     int N_dimensions = input_data[6];
    //     bool importance_sampling = input_data[7];
    //     bool numerical_double_derivative = false;
    //     bool interactions = input_data[8];
    //     int MC_cycles = pow(10, input_data[9]);

    //     string filename;
    //     if (importance_sampling){
    //         importance_sampling = true;
    //         filename = "output/N=" + to_string(N_particles) +
    //                    "_d=" + to_string(N_dimensions) + "_ana_IS.csv";
    //         cout << "Importance sampling" << endl;
    //     }
    //     else{
    //         importance_sampling = false;
    //         cout << "No importance sampling" << endl;
    //         filename = "output/N=" + to_string(N_particles) +
    //                    "_d=" + to_string(N_dimensions) + "_ana.csv";
    //     }

    //     ofstream ofile;
    //     ofile.open(filename);
    //     ofile << "MC,N,d,alpha,energy,variance" << endl;
    //     vec alpha_values = linspace(0.1, 1.0, 10);
    //     for (double alpha : alpha_values)
    //     {
    //         VMCWalker walker(
    //             alpha,
    //             beta,
    //             gamma,
    //             step,
    //             time_step,
    //             hard_core_radius,
    //             ndd_h,
    //             N_particles,
    //             N_dimensions,
    //             importance_sampling,
    //             numerical_double_derivative,
    //             interactions);
            
    //         vec result = walker.walk(MC_cycles);
    //         ofile << MC_cycles << "," << N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
    //     }
    // }
    if (task == "numerical") // TODO: fix variance calculation, all zeros currently.
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
        bool numerical_double_derivative = input_data[6];
        bool interactions = false;
        int MC_cycles = pow(10, input_data[7]);

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

            vec result = walker.walk(MC_cycles, "", "");
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
        int MC_cycles = pow(10, input_data[10]);
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

        // int MC_cycles = input_data[10];
        int MC_cycles = pow(10, input_data[10]);
        int n_walkers = input_data[11];
        string energy_filename = "energy_statistics_N=" + to_string(N_particles);

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

    //   else if(task == "interactions"){
    //     if (input_data.size() != 6){
    //         cout << "Wrong number of parameters in config file" << endl;
    //         exit(1);
    //     }
    //     int MC_cycles = pow(10, input_data[0]);
    //     int N_particles = input_data[1];
    //     int N_dimensions = input_data[2];
    //     bool importance_sampling = (bool) input_data[3];
    //     double time_step = input_data[4];
    //     double step = input_data[5];
    //     bool interactions = true;
    //     double gamma = 2.82843;
    //     double beta = 2.82843;
    //     double hard_core_radius = 0.0043;
    //     if (input_data[3] == 0){
    //           importance_sampling = false;
    //           cout << "No importance sampling" << endl;
    //           filename = "output/N=" + to_string(N_particles) +
    //                             "_d=" + to_string(N_dimensions) + "_int.csv";
    //       }
    //     else{
    //           importance_sampling = true;
    //           filename = "output/N=" + to_string(N_particles) +
    //                             "_d=" + to_string(N_dimensions) + "_int_IS.csv";
    //           cout << "Importance sampling" << endl;
    //     }
    //     ofstream ofile;
    //     ofile.open(filename);
    //     ofile << "MC,N,d,alpha,energy,variance" << endl;
    //     vec alpha_values = linspace(0.1, 1.0, 10);
    //     vec step_sizes = step*linspace(0.4, 1.0, 10);
    //     vec time_steps = {0.7, 0.7, 0.65, 0.6, 0.5, 0.2, 0.1, 0.07, 0.04, 0.01};
    //     //vec time_steps = time_step*linspace(0.5, 1.5, 10);
    //     vec result = monte_carlo(0.5, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, interactions, gamma, beta, hard_core_radius);
    //     // int N_alpha = alpha_values.size();
    //     // for (int i = 0; i < N_alpha; i++){
    //     //     double alpha = alpha_values[i];
    //     //     double step = step_sizes[N_alpha - i - 1];
    //     //     double time_step = time_steps[i];
    //     //     vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, interactions, gamma, beta, hard_core_radius);
    //     //     ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
    //     // }
    //   }

    //   else if(task == "interactions_gradient"){
    //     if (input_data.size() != 7){
    //         cout << "Wrong number of parameters in config file" << endl;
    //         exit(1);
    //     }
    //     int MC_cycles = pow(10, input_data[0]);
    //     int N_particles = input_data[1];
    //     int N_dimensions = input_data[2];
    //     bool importance_sampling = (bool) input_data[3];
    //     double time_step = input_data[4];
    //     double learning_rate = pow(10, input_data[5]);
    //     int max_iter = input_data[6];
    //     bool interactions = true;
    //     double gamma = 2.82843;
    //     double beta = 2.82843;
    //     double hard_core_radius = 0.0043;
    //     minimize_parameters(MC_cycles, -2.0, N_particles, N_dimensions, importance_sampling, time_step, learning_rate, max_iter, interactions, gamma, beta, hard_core_radius);
    //   }
    //   else if(task == "parallel"){  // TODO: set up config file for this.
    //                                 // This task should be used to run timing tests.
    //                                 // Worth noting: each walker has individual burn-in period.
    //         int MC_cycles = 1'000'000;
    //         int N_particles = 10;
    //         int N_dimensions = 3;
    //         double time_step = 0.5;
    //         double step = 1.5;
    //         bool importance_sampling = false;
    //         bool interactions = false;
    //         bool numerical_double_derivative = false;
    //         double gamma = 2.82843;
    //         double beta = 2.82843;
    //         double ndd_h = 0.0001;
    //         double hard_core_radius = 0.0043;
    //         int N_cycles = 100;
    //         double alpha = 0.7;
    //         int N_walkers = 10;

    //         // Run timing tests

    //         // time parallelized monte carlo
    //         auto start1 = chrono::high_resolution_clock::now();
    //         vec v1 = monte_carlo_parallelized(N_walkers, alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, numerical_double_derivative, ndd_h, interactions, gamma, beta, hard_core_radius);
    //         auto finish1 = chrono::high_resolution_clock::now();
    //         chrono::duration<double> elapsed1 = finish1 - start1;

    //         // time serial monte carlo
    //         auto start2 = chrono::high_resolution_clock::now();
    //         vec v2 = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, numerical_double_derivative, ndd_h, interactions, gamma, beta, hard_core_radius);
    //         auto finish2 = chrono::high_resolution_clock::now();
    //         chrono::duration<double> elapsed2 = finish2 - start2;

    //         cout << "Time parallelized: " << elapsed1.count() << " s" << endl;
    //         cout << "Time serial: " << elapsed2.count() << " s" << endl;

    //         // Verify that the results are the same
    //         cout << "v1 = " << v1 << endl;
    //         cout << "v2 = " << v2 << endl;
    //     }
    //   else{
    //       cout << "Unknown task" << endl;
    //   }
    //   return 0;
    // }

    // double alpha, int MC_cycles, int N_walkers, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius)
    // double alpha, int MC_cycles, int N_walkers, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius)
    // vec v = monte_carlo_parallelized(alpha, MC_cycles, N_walkers, step, N_particles, N_dimensions, importance_sampling, time_step, numerical_double_derivative, ndd_h, interactions, gamma, beta, hard_core_radius);
}