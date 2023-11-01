// Copyright 2023 Brian Howell
// MIT License
// Project: BayesOpt

#ifndef SRC_COMMON_H_
#define SRC_COMMON_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <algorithm>


// input/output BOpt data structure
typedef struct bopt {
    float temp;    // temperature
    float rp;      // particle radius
    float vp;      // volume fraction
    float uvi;     // initial UV intensity
    float uvt;     // uv exposure time
    float obj;     // objective
} bopt;

// sim settings
typedef struct sim {
    int method;         // simulation method
    int save_voxel;     // save voxel data
    int save_density;   // save density data
    bool bootstrap;     // bootstrap
    int time_stepping;  // representing with dt/node pair

    static constexpr double DT[7]   = {1.5e-3,
                                       7.5e-4,
                                       3.5e-4,
                                       2.5e-4,
                                       1.5e-4,
                                       9.25e-5,
                                       5.9e-5};
    static constexpr int    NODE[7] = {11,
                                       15,
                                       21,
                                       25,
                                       31,
                                       41,
                                       51};

    float  tfinal;   // final time
    double dt;              // time step
    int    node;            // number of nodes

    // default constructor
    sim() {
        method        = 2;      // forward euler | 1: backward euler | 2: trap
        save_voxel    = 0;      // save voxel data
        save_density  = 0;      // save density data
        bootstrap     = false;  // bootstrap
        time_stepping = 1;      // representing with dt/node pair
        dt            = DT[time_stepping];
        node          = NODE[time_stepping];
        tfinal        = 30.;
    }

    // overload constructor
    sim(int method,
        int save_voxel,
        int save_density,
        int bootstrap,
        int time_stepping)
        : method(method),
          save_voxel(save_voxel),
          save_density(save_density),
          bootstrap(bootstrap),
          time_stepping(time_stepping) {
        // set dt and node spacing
        dt   = DT[time_stepping];
        node = NODE[time_stepping];
        tfinal = 30.;
    }

    void update_time_stepping_values() {
        dt = DT[time_stepping];
        node = NODE[time_stepping];
    }
} sim;

#endif  // SRC_COMMON_H_
