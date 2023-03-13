// #include "../include/initialization.hpp"
#include "../include/vmc_walker.hpp"

using namespace arma;
using namespace std;

void VMCWalker::set_initial_state_interaction(){
    int spread = 50;
    position = spread*hard_core_radius*(mat(N_dimensions, N_particles).randu() - 0.5); //initialize all particles at random positions
    relative_position = trimatl(relative_position, 1); //strict lower triangular matrix to store relative positions
    for (int i = 0; i < N_particles; i++){
        vec r_i = position.col(i);
        for (int j = i + 1; j < N_particles; j++){
            double rel_pos = norm(r_i - position.col(j));
            while (rel_pos < hard_core_radius){
                //if two particles are too close, get new random positions
                position.col(j) = spread*hard_core_radius*(vec(N_dimensions).randu() - 0.5);
                rel_pos = norm(r_i - position.col(j));    
            }
            relative_position(j, i) = rel_pos;
        }
    }
}

void VMCWalker::set_initial_state_no_interaction(){
    position = mat(N_dimensions, N_particles).fill(0.0);
}

void VMCWalker::set_initial_state(){
    if (interactions){
        set_initial_state_interaction();
    }
    else{
        set_initial_state_no_interaction();        
    }
}

void VMCWalker::burn_in_interaction()
{
    double fraction = 0.0;
    int N_equilibration = 1000;
    if (importance_sampling){
        while(0.8 > fraction | fraction > 0.95){
            accepted_moves = 0;
            for (int j = 0; j < N_equilibration; j++){
                for (int k = 0; k < N_particles; k++){
                    sampling(k);
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
    }
    else{
        for (int j = 0; j < N_equilibration; j++){
            for (int k = 0; k < N_particles; k++){
                sampling(k);
            }
        }
        fraction = (accepted_moves + 0.0)/(N_particles*N_equilibration);
    }
    cout << "Equilibration done" << endl;
    cout << "Fraction of accepted moves: " << (accepted_moves + 0.0)/(N_particles*N_equilibration) << endl;    
    cout << "Final time step: " << time_step << endl;
}

void VMCWalker::initialize(){
    set_initial_state();
    new_position = position;
    relative_position_new = relative_position;

    // With interactions, the system needs to equilibriate before the simulation starts
    if (interactions){
        burn_in_interaction();
    }
}
