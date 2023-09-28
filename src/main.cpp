// Description: main file for ugap_simulation
//              reaction diffusion simulation of 
//              urethane grafted acrylate polymer (UPAP)
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
    const bool mthread = false; 
    const int NODE = 21;
    const double DT = 1e-3;

    std::string file_path = "/Users/brianhowell/Desktop/Berkeley/MSOL/ugap_simulation/output/"; 

    // const int NODE = 51;
    // const double DT = 5e-5;
    float TFINAL = 30.;
    
    // parameter sweeps
    const int num_sweeps          = 7;
    float I_UV[num_sweeps]        = {2, 6, 10, 50, 100};
    float temp_amb[num_sweeps]    = {298.15, 303.15, 315.15, 325.15, 335.15};

    double dt_sweep[num_sweeps]   = {1.5e-3, 7.5e-4,   3.5e-4, 2.5e-4, 1.5e-4, 9.25e-5, 5.9e-5};
    double node_sweep[num_sweeps] = {11,     15,     21,   25,     31,     41,   51};


    int   save_voxel = 1;  // 0: off | 1: on
    int   method     = 2;  // 0: forward euler | 1: backward euler | 2: trap
    int   sim_count  = 1;

    // SINGLE SIMULATION
    auto start = std::chrono::high_resolution_clock::now();
    Voxel VoxelSystem1(TFINAL, 
                       dt_sweep[1],
                       node_sweep[1],
                       method,
                       340.484,
                       76.1647,
                       3.0,
                       file_path,
                       mthread);
    VoxelSystem1.ComputeParticles(5.51541e-5, 0.596627);
    VoxelSystem1.Density2File();
    VoxelSystem1.Simulate(method, save_voxel);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;

    std::cout << " --- Simulation time: " << duration / 60 << "min ---" << std::endl;
    std::cout << " --- ----------------------------- ---" << std::endl;

    return 0;
}
