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
    float TFINAL = 30.;
    
    // parameter sweeps
    const int num_sweeps          = 5;
    float I_UV[num_sweeps]        = {2, 6, 10, 50, 100};
    float temp_amb[num_sweeps]    = {298.15, 303.15, 315.15, 325.15, 335.15};
    double dt_sweep[num_sweeps]   = {1e-5, 5e-5, 1e-4, 5e-4, 1e-3};
    double node_sweep[num_sweeps] = {51, 41, 31, 21, 11};


    int   save_voxel = 0;  // 0: off | 1: on
    int   method     = 2;  // 0: forward euler | 1: backward euler | 2: trap
    int   sim_count  = 1;

    for (int i = 0; i < num_sweeps; i++){

        std::cout << " --- --------- Simulation: " << i << "----------- ---" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        // SINGLE SIMULATION
        // Voxel VoxelSystem1(I_UV[2], TFINAL, 5e-4, NODE, sim_count);

        // UV SWEEP
        // Voxel VoxelSystem1(I_UV[i], TFINAL, DT, NODE, sim_count);
        
        // DT SWEEP
        // Voxel VoxelSystem1(I_UV[2], TFINAL, dt_sweep[i], NODE, sim_count);
        
        // TEMP SWEEP
        Voxel VoxelSystem1(I_UV[2], TFINAL, DT, NODE, sim_count, temp_amb[i]);

        VoxelSystem1.ComputeParticles(0.00084 / 10, 0.70);
        VoxelSystem1.Density2File();

        VoxelSystem1.Simulate(method, save_voxel);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
        std::cout << " --- Simulation time: " << duration / 60 << "min ---" << std::endl;
        std::cout << " --- ----------------------------- ---" << std::endl;

        sim_count++; 
    }

    std::cout << "\nHello, World!" << std::endl;

    return 0;
}
