#include "../include/interaction.hpp"
using namespace arma;
using namespace std;


double local_energy_interaction(mat& position, mat& relative_position, double alpha, double beta, double gamma, double hard_core_radius){
    double kinetic_energy = 0;
    double potential_energy = 0;
    assert(relative_position.is_trimatl());
    
    for (int k = 0; k < position.n_cols; k++){
        vec r_k = position.col(k);
        kinetic_energy += -2*alpha*(beta + 2);
        kinetic_energy += 4*alpha*alpha*(r_k(0)*r_k(0) + r_k(1)*r_k(1) + beta*beta*r_k(2)*r_k(2));
        if (hard_core_radius > 0){

            vec kinetic_energy_contrib = vec(3, fill::zeros);
            for (int l = 0; l < position.n_cols; l++){
                if (l != k){
                    vec r_l = position.col(l);
                    double r_kl = max(relative_position.col(k)(l), relative_position.col(l)(k));
                    assert(r_kl != 0);
                    assert(r_kl > hard_core_radius);
                    double H_kl = hard_core_radius/(r_kl*r_kl*(r_kl-hard_core_radius));
                    kinetic_energy_contrib += H_kl*(r_k - r_l);
                }
            }
            vec tmp = {r_k(0), r_k(1), beta*r_k(2)};
            kinetic_energy += - 4*alpha*dot(tmp, kinetic_energy_contrib);

            double kinetic_energy_contrib2 = 0.0;
            for (int i = 0; i < position.n_cols; i++){
                if (i != k){
                    vec r_i = position.col(i);
                    double r_ki = max(relative_position.col(k)(i), relative_position.col(i)(k));
                    assert(r_ki != 0);
                    assert(r_ki > hard_core_radius);
                    double H_ki = hard_core_radius/(r_ki*r_ki*(r_ki-hard_core_radius));
                    for (int j = 0; j < position.n_cols; j++){
                        if (j != k){
                            vec r_j = position.col(j);
                            double r_kj = max(relative_position.col(k)(j), relative_position.col(j)(k));
                            assert(r_kj != 0);
                            assert(r_kj > hard_core_radius);
                            double H_kj = hard_core_radius/(r_kj*r_kj*(r_kj-hard_core_radius));
                            kinetic_energy_contrib2 += dot(r_k - r_i, r_k - r_j) * H_kj * H_ki;
                        }
                    } 
                    // kinetic_energy_contrib2 *= H_ki;               
                }
            }
            kinetic_energy += kinetic_energy_contrib2;
            
            for (int l = 0; l < position.n_cols; l++){
                if (l != k){
                    double r_kl = max(relative_position.col(k)(l), relative_position.col(l)(k));
                    vec r_l = position.col(l);
                    assert(r_kl != 0);
                    assert(r_kl > hard_core_radius);
                    double H_kl = hard_core_radius/(r_kl*(r_kl-hard_core_radius));
                    kinetic_energy -= H_kl*H_kl;
                }
            }
            potential_energy += r_k(0)*r_k(0) + r_k(1)*r_k(1) + gamma*gamma*r_k(2)*r_k(2);
            
        }
    }
    kinetic_energy *= -0.5;
    potential_energy *= 0.5;

    return kinetic_energy + potential_energy;
}

double probability_ratio_interaction(mat& old_position, mat& new_position, mat& relative_position, mat& relative_position_new, double alpha, double beta, double hard_core_radius, int particle_index){
    double ratio = 1.0;
    vec r_k_old = old_position.col(particle_index);
    vec r_k_new = new_position.col(particle_index);

    double exponent = (r_k_new(0)*r_k_new(0) + r_k_new(1)*r_k_new(1) + beta*r_k_new(2)*r_k_new(2)) 
                - (r_k_old(0)*r_k_old(0) + r_k_old(1)*r_k_old(1) + beta*r_k_old(2)*r_k_old(2));
    ratio *= exp(-2*alpha*exponent);

    if (hard_core_radius > 0){
        for (int j = 0; j < old_position.n_cols; j++){
            if (j != particle_index){
                double r_old = max(relative_position.col(particle_index)(j), relative_position.col(j)(particle_index));
                double r_new = max(relative_position_new.col(particle_index)(j), relative_position_new.col(j)(particle_index));
                double H = r_old*(r_new - hard_core_radius)/(r_new*(r_old - hard_core_radius));
                ratio *= H*H;
            }
        }
    }
    //cout << "ratio: " << ratio << endl;
    return ratio;
}

vec quantum_force_interaction(mat& position, mat& relative_position, double alpha, double beta, double hard_core_radius, int particle_index){
    vec r_k = position.col(particle_index);
    vec qf = {r_k[0], r_k[1], beta*r_k[2]};
    qf *= -4*alpha;
    if (hard_core_radius > 0){
        for (int l = 0; l < position.n_cols; l++){
            if (l != particle_index){
                vec r_l = position.col(l);
                double r_kl = max(relative_position.col(particle_index)(l), relative_position.col(l)(particle_index));
                // assert(r_kl != 0);
                // assert(r_kl > hard_core_radius);
                double H_kl = hard_core_radius/(r_kl*r_kl*(r_kl-hard_core_radius)); 
                qf += 2*H_kl*(r_k - r_l);
            }
        }
    }
    // cout << "qf: " << qf << endl;
    return qf;
}
