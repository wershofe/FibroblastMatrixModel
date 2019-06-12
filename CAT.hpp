//
//  chemoattr_pde.hpp
//  
//
//  Created by Xiao Fu on 16/08/2017.
//
//CAT

#ifndef CAT_hpp
#define CAT_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
//#include "gnuplot-iostream/gnuplot-iostream.h"  // git clone https://github.com/dstahlke/gnuplot-iostream.git; add linker when complie: -lboost_iostreams -lboost_system -lboost_filesystem

typedef struct {int i,j;} iPair;

const int nblock = 20;            /* number of blocks in each dimension */
const double CHEMOATTR_MAX = 1;     /* max chemoattractant level */
const double tt = 1;             /* total time */
const double dt = 0.001;            /* time step */
const int    N_DT_SAVE_FIG = 1000; /* # of time steps to save figures */

const double k_dif = 0.1;           /* diffusion coefficient */
const double k_dec = .0025;         /* decay/degredation rate constant: 0.025, 0.0025 */
const double k_sec = 10.;           /* secretion rate constant */
const double k_upt = 0;           /* uptake rate constant */


using namespace std;

// declare function to simulate one step of chemoattractant diffusion-reaction
void singleStep(double chemoattr[][nblock], /* 2D array of chemoattractant level */
                const double xdim,          /* x dimension */
                const double ydim,          /* y dimension */
                vector<iPair> source,       /* pixel location of chemoattractant sources (cancer cells) */
                vector<iPair> sink          /* pixel location of chemoattractant sinks (CAF) */
                );

// declare secretion method
void do_secretion(double chemoattr[][nblock], vector<iPair> source);

// declare decay method
void do_decay(double chemoattr[][nblock]);

// declare uptake method
void do_uptake(double chemoattr[][nblock], vector<iPair> sink);

// declare simple finite difference method
void do_finite_diff_simple(double chemoattr_tmp[][nblock], double chemoattr[][nblock], double Kx, double Ky);

#endif /* chemoattr_pde_hpp */
