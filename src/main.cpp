//
// reaction diffusion simulation of urethane grafted acrylate 
// polymer resins with custom additive particles for 3D printing.
// Created by Brian Howell on 02/24/23
// MSOL UC Berkeley
// bhowell@berkeley.edu

#include <iostream>
#include "Voxel.h"

int main() {
    /*
        coarse simulation:
            - node = 21
            - dt   = 1e-3
            - VoxelSystem1.Simulate(10.0, TFINAL, DT, 1);
            - outputs 300 .vtk concentration files

        fine simulation
            - node = 51
            - dt   = 5e-5
            - VoxelSystem1.Simulate(10.0, TFINAL, DT, 1);
            - outputs 600 .vtk concentration files
            
    */

    const int NODE = 21;
    const double DT = 1e-3;
    // const int NODE = 51;
    // const double DT = 5e-5;
    double TFINAL = 30.;
    

    auto start = std::chrono::high_resolution_clock::now();
    Voxel VoxelSystem1(NODE);
    VoxelSystem1.ComputeParticles(0.00084 / 10, 0.70);
    VoxelSystem1.Density2File();

    VoxelSystem1.Simulate(10.0, TFINAL, DT, 2);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << " --- Simulation time: " << duration / 60 << "min ---" << std::endl;

    std::cout << "\nHello, World!" << std::endl;

    return 0;
}
