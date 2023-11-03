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
    // std::string file_path = "/Users/brianhowell/Desktop/Berkeley/MSOL/ugap_simulation/output/";
    std::string file_path = "/home/brian/Documents/brian/ugap_simulation/output/";

    // SINGLE SIMULATION
    bopt default_bopt;
    default_bopt.temp = 303.15;
    default_bopt.rp   = 0.00084 / 10;
    default_bopt.vp   = 0.7;
    default_bopt.uvi  = 10.;
    default_bopt.uvt  = 30.;
    
    sim default_sim;
    default_sim.time_stepping = 0;
    default_sim.update_time_stepping_values();

    const bool mthread = false; 
    int   save_voxel = 1;  // 0: off | 1: on

    // auto start = std::chrono::high_resolution_clock::now();
    // std::cout << "===== RUNNING DEFAULT PARAMETERS =====" << std::endl;
    // Voxel VoxelSystem1(default_sim.tfinal,  // tot sim time
    //                    default_sim.dt,      // time step
    //                    default_sim.node,    // num nodes
    //                    default_sim.method,  // sim id
    //                    default_bopt.temp,   // amb temp
    //                    default_bopt.uvi,    // uv intensity
    //                    default_bopt.uvt,    // uv exposure time
    //                    file_path,
    //                    mthread);
    // VoxelSystem1.computeParticles(default_bopt.rp, default_bopt.vp);
    // VoxelSystem1.density2File();
    // VoxelSystem1.simulate(default_sim.method, save_voxel);
    // double default_objective = VoxelSystem1.getObjective();
    // std::cout << "default_objective: " << default_objective << std::endl;
    // auto stop = std::chrono::high_resolution_clock::now();
    // auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;

    // std::cout << " --- Simulation time: " << duration / 60 << "min ---" << std::endl;
    // std::cout << " --- ----------------------------- ---" << std::endl;

    std::cout << "\n===== RUNNING OPT PARAMETERS =====" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    bopt opt_bopt;
    opt_bopt.temp = 348.281; 
    opt_bopt.rp   = 5.74817e-05;
    opt_bopt.vp   = 0.798616;
    opt_bopt.uvi  = 60.3975;
    opt_bopt.uvt  = 1.11566;

    Voxel VoxelSystem2(default_sim.tfinal,  // tot sim time
                       default_sim.dt,      // time step
                       default_sim.node,    // num nodes
                       default_sim.method,  // sim id
                       default_bopt.temp,   // amb temp
                       default_bopt.uvi,    // uv intensity
                       default_bopt.uvt,    // uv exposure time
                       file_path,
                       mthread);
    VoxelSystem2.computeParticles(default_bopt.rp, default_bopt.vp);
    VoxelSystem2.density2File();
    VoxelSystem2.simulate(default_sim.method, save_voxel);

    double opt_objective = VoxelSystem2.getObjective();
    std::cout << "opt_objective: " << opt_objective << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;

    std::cout << " --- Simulation time: " << duration / 60 << "min ---" << std::endl;
    std::cout << " --- ----------------------------- ---" << std::endl;

    // std::cout << "\n==== ANALYSIS ====" << std::endl;
    // std::cout << "pcnt improvement: " << 100 * (default_objective - opt_objective) / default_objective << std::endl;

    return 0;
}
