// #include "../include/initialization.hpp"
#include "../include/vmc_walker.hpp"

using namespace arma;
using namespace std;

void VMCWalker::set_initial_state_interaction(){
    position = (mat(N_dimensions, N_particles).randn()); //initialize all particles at random positions
    
    relative_position = mat(N_particles, N_particles).fill(0.0);
    relative_position = trimatl(relative_position, 1); //strict lower triangular matrix to store relative positions
    for (int i = 0; i < N_particles; i++){
        vec r_i = position.col(i);
        for (int j = i + 1; j < N_particles; j++){
            double rel_pos = norm(r_i - position.col(j));
            while (rel_pos < hard_core_radius){
                //if two particles are too close, get new random positions
                position.col(j) = vec(N_dimensions).randn();
                rel_pos = norm(r_i - position.col(j));    
            }
            relative_position(j, i) = rel_pos;
        }
    }
}

void VMCWalker::set_initial_state_no_interaction(){
    position = mat(N_dimensions, N_particles).fill(0.0);
    relative_position = mat(N_particles, N_particles).fill(0.0); // not in use - for compatibility only
}

void VMCWalker::set_initial_state(){
    if (interactions && hard_core_radius > 0){ 
        set_initial_state_interaction();
    }
    else{
        set_initial_state_no_interaction();        
    }
}

void VMCWalker::adjust_timestep_importance_sampling()
{
    double fraction = 0.0;
    int N_equilibration = 1'000;
    while(0.80 > fraction | fraction > 0.90){
        accepted_moves = 0;
        for (int j = 0; j < N_equilibration; j++){
            for (int k = 0; k < N_particles; k++){
                sampling(k);
                }
        }
        fraction = (accepted_moves + 0.0)/(N_particles*N_equilibration);
        if ( fraction < 0.80){
            //if the acceptance rate is too low, decrease the step size
            time_step *= 0.8;
        }
        else if (fraction > 0.90){
            //if the acceptance rate is too high, increase the step size
            time_step *= 1.2;
        }
    }

}

void VMCWalker::adjust_step_brute_force_sampling(){
    double fraction = 0.0;
    int N_equilibration = 1000;
    while(0.4 > fraction | fraction > 0.6){
        accepted_moves = 0;
        for (int j = 0; j < N_equilibration; j++){
            for (int k = 0; k < N_particles; k++){
                sampling(k);
                }
        }
        fraction = (accepted_moves + 0.0)/(N_particles*N_equilibration);
        if (fraction < 0.4){
            //if the acceptance rate is too low, decrease the step size
            step *= 0.8;
        }
        else if (fraction > 0.6){
            //if the acceptance rate is too high, increase the step size
            step *= 1.2;
        }
    }
}

void VMCWalker::burn_in(){
    int N_equilibration = 10'000;
    accepted_moves = 0;
    for (int j = 0; j < N_equilibration; j++){
        for (int k = 0; k < N_particles; k++){
            sampling(k);
        }
    }
    double fraction = (accepted_moves + 0.0)/(N_particles*N_equilibration);
}

void VMCWalker::initialize(){
    set_initial_state();
    new_position = position;
    relative_position_new = relative_position;
    if (adjust_step_automatically){
        if (importance_sampling){
            adjust_timestep_importance_sampling();
        }
        else{
            adjust_step_brute_force_sampling();;
        }
    }
    burn_in();
    if (adjust_step_automatically){
        if (importance_sampling){
            adjust_timestep_importance_sampling();
        }
        else{
            adjust_step_brute_force_sampling();;
        }
    }
}
