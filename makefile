build:
	g++ main.cpp src/vmc.cpp src/interaction.cpp src/no_interaction.cpp src/sampling.cpp -I include/ -larmadillo -O2 -o main.exe -fopenmp
