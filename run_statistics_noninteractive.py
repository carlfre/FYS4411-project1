import numpy as np
import subprocess

alphas = 1/2 + np.array([[3.604e-4, -1.150e-5, -1.140e-5, 2.150e-5],
          [1.588e-4, -7.900e-6, 9.200e-6, 2.300e-6],
          [-9.900e-5, -1.210e-5, 1.180e-5, 5.620e-5],])

N_values = [1, 10, 100, 500]
d_values = [1, 2, 3]

for j in range(len(N_values)):
    for i in range(len(d_values)):
        N = N_values[j]
        d = d_values[i]
        alpha = alphas[i][j]

        msg = f"""# alpha
{alpha}
# beta 
2.82843
# gamma
2.82843
# step
1.0
# time_step
0.1
# hard_core_radius
0.0043
# N_particles
{N}
# N_dimensions
{d}
# importance_sampling
1
# interactions
0
# log_2 (MC_cycles) - NOTICE: base 2 here, to work better with blocking!
23
#12
# N_walkers
16"""

        # write to file
        with open("configs/statistics.txt", "w") as f:
            f.write(msg)
        
        # run program
        subprocess.run(["./main.exe", "statistics"])