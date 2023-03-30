#pragma once 
#include <iostream>
#include <iomanip> 
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

/* Monte Carlo walker class. 
After initializing object x with wanted parameters, 
run monte carlo simulation with x.walk(MC_cycles). 
Filenames for density problem can also be passed - then position is written to file. 
Similarly, filename for energy estimation/statistics problemcan also be passed - then energy is written to file. 
*/
class VMCWalker{
    public:
        arma::mat position;
        arma::mat new_position;
        arma::mat relative_position;
        arma::mat relative_position_new;
        arma::vec qf;
        arma::vec qf_new;
        // arma::vec gaussian_random;
        arma::vec stepping_vector;
        int N_particles;
        int N_dimensions;
        double alpha;
        double beta;
        double gamma;
        double hard_core_radius;
        double time_step;
        double step;
        double D;
        double random_number;
        double ndd_h;
        int accepted_moves;
        bool interactions;
        bool importance_sampling;
        bool numerical_double_derivative;
        bool adjust_step_automatically;
        // random number generation
        // mersenne twister
        std::mt19937_64 generator;
        // uniform distribution
        std::uniform_real_distribution<double> rng_double;
        // normal distribution
        std::normal_distribution<double> rng_normal;


        // vmc_walker.cpp

        /* Set parameters for Monte Carlo simulation.*/
        VMCWalker(
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

        );

        /* Can be used to disable step/time step size adjustment*/
        void set_step_adjustment(bool use_adjustment);

        /* Run Monte Carlo simulation.*/
        arma::vec walk(int MC_cycles, std::string density_filename, std::string energy_filename);
        arma::vec walk(int MC_cycles);

        /* Gradient descent with momentum to find alpha-value that minimizes energy*/
        arma::vec minimize_parameters(int MC_cycles, double learning_rate, int max_iter, std::string filename);
        arma::vec minimize_parameters(int MC_cycles, double learning_rate, int max_iter);


        // initialization.cpp

        /* Initialize position and relative position */
        void set_initial_state_interaction();
        void set_initial_state_no_interaction();
        void set_initial_state();

        /* Automatic adjustment of timestep to achieve good acceptance ratio */
        void adjust_timestep_importance_sampling();
        void adjust_step_brute_force_sampling();

        /* Run 10,000 MC-cycles of equilibration to reach higher-probability state */
        void burn_in();

        /* Run initialization the above initialization functions */
        void initialize();

        // sampling.cpp

        /* Run 1 step of metropolis-algorithm for particle k. 
        Ie. we draw a new position, and check if we accept it. 
        4 separate implementations, for cases with/without interactions
        and with/without importance sampling */
        void importance_sampling_with_interactions(int k);
        void importance_sampling_without_interactions(int k);
        void brute_force_sampling_with_interactions(int k);
        void brute_force_sampling_without_interactions(int k);

        /* Calls one of four functions above, according to parameters of the run. */
        void sampling(int k);
};