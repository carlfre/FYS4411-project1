#include "include/vmc.hpp"
#include "include/interaction.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
    if (argc != 3){
        cout << "Usage: task config_file_name " << endl;
        exit(1);
    }
    string task = string(argv[1]);
    string filename = "configs/" + string(argv[2]);
    // A vector of vectors to store the the values in the input file
    vector<double> input_data;
    fstream myfile;
    myfile.open(filename, ios::in);
    if (myfile.is_open()){  // This checks that the file was opened OK
        string line;
        double param;
        // Read file line by line
        while (getline(myfile, line)){
        // Skip lines with "#" at the first position
        if (line.at(0) == '#'){
            continue;
        }
        else{
            // Parse the string (line) and interpret it as one variable
            stringstream mysstream(line);
            mysstream >> param;
            input_data.push_back(param);
      }
    }
  }
  else{
    cout << "Unable to open the file " << filename << endl;
    exit(1);
  }
  // Close the input file
  myfile.close();
  if (task == "numerical"){
      double step = 1.0;
      if (input_data.size() != 6){
          cout << "Wrong number of parameters in config file" << endl;
          exit(1);
      }
      int N_particles = input_data[0];
      int N_dimensions = input_data[1];
      int MC_cycles = pow(10, input_data[2]);
      double ndd_h = input_data[3];
      bool importance_sampling;
      double time_step = input_data[5];
      string filename;
      if (input_data[4] == 0){
          importance_sampling = false;
          cout << "No importance sampling" << endl;
          filename = "output/N=" + to_string(N_particles) +
                            "_d=" + to_string(N_dimensions) + "_num.csv";   
      }
      else{
          importance_sampling = true;
          filename = "output/numerical_N=" + to_string(N_particles) +
                            "_d=" + to_string(N_dimensions) + "num_IS.csv";
          cout << "Importance sampling" << endl;
      }
      ofstream ofile;
      ofile.open(filename);
      ofile << "MC,N,d,alpha,energy,variance" << endl;
      vec alpha_values = linspace(0.1, 1.0, 10);
      for (double alpha : alpha_values){
          vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, true, ndd_h);
          ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
       }
 
  }
  else if(task == "analytical"){
      double step = 1.5;
      if (input_data.size() != 5){
          cout << "Wrong number of parameters in config file" << endl;
          exit(1);
      }
      int MC_cycles = pow(10, input_data[0]);
      int N_particles = input_data[1];
      int N_dimensions = input_data[2];
      double time_step = input_data[4];
      bool importance_sampling = (bool) input_data[3];
      string filename;
      if (input_data[3] == 0){
          importance_sampling = false;
          cout << "No importance sampling" << endl;
          filename = "output/N=" + to_string(N_particles) +
                            "_d=" + to_string(N_dimensions) + "_ana.csv";   
      }
      else{
          importance_sampling = true;
          filename = "output/N=" + to_string(N_particles) +
                            "_d=" + to_string(N_dimensions) + "_ana_IS.csv";
          cout << "Importance sampling" << endl;
      }
      ofstream ofile;
      ofile.open(filename);
      ofile << "MC,N,d,alpha,energy,variance" << endl;
      vec alpha_values = linspace(0.1, 1.0, 10);
      for (double alpha : alpha_values){
          vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step);
          ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
       }
  }
  else if(task == "gradient"){
      if (input_data.size() != 7){
          cout << "Wrong number of parameters in config file" << endl;
          exit(1);
      }
      double step = 1.0;
      int MC_cycles = pow(10, input_data[0]);
      int N_particles = input_data[1];
      int N_dimensions = input_data[2];
      double time_step = input_data[4];
      double learning_rate = pow(10, input_data[5]);
      int max_iter = input_data[6];
      bool importance_sampling = (bool) input_data[3];
      if (input_data[4] == 0){
          importance_sampling = false;
          cout << "No importance sampling" << endl;
      }
      else if (input_data[4] == 1){
          importance_sampling = true;
          cout << "Importance sampling" << endl;
      }
      minimize_parameters(MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, learning_rate, max_iter);
  }
  else if(task == "interactions"){
    if (input_data.size() != 6){
        cout << "Wrong number of parameters in config file" << endl;
        exit(1);
    }
    int MC_cycles = pow(10, input_data[0]);
    int N_particles = input_data[1];
    int N_dimensions = input_data[2];
    bool importance_sampling = (bool) input_data[3];
    double time_step = input_data[4];
    double step = input_data[5];
    bool interactions = true;
    double gamma = 2.82843;
    double beta = 2.82843;
    double hard_core_radius = 0.0043;
    if (input_data[3] == 0){
          importance_sampling = false;
          cout << "No importance sampling" << endl;
          filename = "output/N=" + to_string(N_particles) +
                            "_d=" + to_string(N_dimensions) + "_int.csv";   
      }
    else{
          importance_sampling = true;
          filename = "output/N=" + to_string(N_particles) +
                            "_d=" + to_string(N_dimensions) + "_int_IS.csv";
          cout << "Importance sampling" << endl;
    }
    ofstream ofile;
    ofile.open(filename);
    ofile << "MC,N,d,alpha,energy,variance" << endl;
    vec alpha_values = linspace(0.1, 1.0, 10);
    vec step_sizes = step*linspace(0.4, 1.0, 10);
    vec time_steps = time_step*linspace(0.4, 1.0, 10);
    int N_alpha = alpha_values.size(); 
    for (int i = 0; i < N_alpha; i++){
        double alpha = alpha_values[i];
        double step = step_sizes[N_alpha - i - 1];
        double time_step = time_steps[N_alpha - i - 1];
        vec result = monte_carlo(alpha, MC_cycles, step, N_particles, N_dimensions, importance_sampling, time_step, interactions, gamma, beta, hard_core_radius);
        ofile << MC_cycles << "," <<  N_particles << "," << N_dimensions << "," << alpha << "," << result[0] << "," << result[1] << endl;
    }
  }
  else{
      cout << "Unknown task" << endl;
  }
  return 0;
}

