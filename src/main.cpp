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
    const int num_sweeps          = 7;
    float I_UV[num_sweeps]        = {2, 6, 10, 50, 100};
    float temp_amb[num_sweeps]    = {298.15, 303.15, 315.15, 325.15, 335.15};
    //0.0032666666666666673, 0.0016666666666666666, 0.0008166666666666668, 0.0005671296296296297, 0.00036296296296296294, 0.0002041666666666667, 0.00013066666666666668
    // double dt_sweep[num_sweeps]   = {3e-3, 1e-3,   7.5e-4, 5.5e-4, 3.5e-4, 2e-4, 1e-5};
    // dt:  [0.0014949152542372884, 0.0007627118644067797, 0.0003737288135593221, 0.00025953389830508483, 0.00016610169491525423, 9.343220338983053e-05, 5.979661016949154e-05]
    double dt_sweep[num_sweeps]   = {1.5e-3, 7.5e-4,   3.5e-4, 2.5e-4, 1.5e-4, 9.25e-5, 5.9e-5};
    double node_sweep[num_sweeps] = {11,     15,     21,   25,     31,     41,   51};


    int   save_voxel = 0;  // 0: off | 1: on
    int   method     = 2;  // 0: forward euler | 1: backward euler | 2: trap
    int   sim_count  = 1;

    for (int i = 0; i < num_sweeps; i++){

        std::cout << " --- --------- Simulation: " << i << "----------- ---" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        // SINGLE SIMULATION
        Voxel VoxelSystem1(I_UV[2], TFINAL, dt_sweep[i], node_sweep[i], i, temp_amb[1]);

        // UV SWEEP
        // Voxel VoxelSystem1(I_UV[i], TFINAL, dt_sweep[2], node_sweep[2], i, temp_amb[1]);
        
        // CONVERGENCE SWEEP
        // Voxel VoxelSystem1(I_UV[2], TFINAL, dt_sweep[sim_count], node_sweep[sim_count], sim_count, temp_amb[1]);

        // TEMP SWEEP
        // Voxel VoxelSystem1(I_UV[2], TFINAL, dt_sweep[2], node_sweep[2], i, temp_amb[i]);

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
