build:
	g++ main.cpp src/vmc.cpp src/interaction.cpp -I include/ -larmadillo -O2 -o main.exe -fopenmp
