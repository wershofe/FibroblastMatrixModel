//
//  chemoattr_pde.cpp
//  
//
//  Created by Xiao Fu on 16/08/2017.
//
//

#include "CAT.hpp"

void singleStep(double chemoattr[][nblock], /* 2D array of chemoattractant level */
                const double xdim,          /* x dimension */
                const double ydim,          /* y dimension */
                vector<iPair> source,       /* pixel location of chemoattractant sources (cancer cells) */
                vector<iPair> sink          /* pixel location of chemoattractant sinks (CAF) */
)
{
    double dx = xdim/nblock, dy = ydim/nblock;
    double Kx = .5 * k_dif * dt / dx, Ky = .5 * k_dif * dt / dy;
    double chemoattr_tmp[nblock][nblock];
    
    // secretion
    do_secretion(chemoattr, source);
    // uptake by CAF
    //do_uptake(chemoattr, sink);
    // decay/degredation
    do_decay(chemoattr);
    
    // diffusion
    memcpy(chemoattr_tmp[0], chemoattr[0], nblock*nblock*sizeof chemoattr[0][0]);   // copy data to temp array
    do_finite_diff_simple(chemoattr_tmp, chemoattr, Kx, Ky);
    
}

void do_secretion(double chemoattr[][nblock], vector<iPair> source)
{
    for (vector<iPair>::iterator it = source.begin(); it!=source.end(); it++)
    {
        iPair loc = *it;
        if (loc.i < nblock && loc.j < nblock)
        {
            chemoattr[loc.i][loc.j] += k_sec*dt* (CHEMOATTR_MAX-chemoattr[loc.i][loc.j]);
        }
    }
}

void do_decay(double chemoattr[][nblock])
{
    for (int j = 0; j < nblock; j++)
    {
        for (int i = 0; i < nblock; i++)
        {
            chemoattr[i][j] -= k_dec * dt * chemoattr[i][j];
        }
    }
    
}

void do_uptake(double chemoattr[][nblock], vector<iPair> sink)
{
    for (vector<iPair>::iterator it = sink.begin(); it!=sink.end(); it++)
    {
        iPair loc = *it;
        if (loc.i < nblock && loc.j < nblock)
        {
            chemoattr[loc.i][loc.j] -= k_upt*dt* chemoattr[loc.i][loc.j];
        }
    }
}

void do_finite_diff_simple(double chemoattr_tmp[][nblock], double chemoattr[][nblock], double Kx, double Ky)
{
    for (int j = 0; j < nblock; j++)
    {
        for (int i = 0; i < nblock; i++)
        {
            chemoattr[i][j] = (1 - 2*Kx - 2*Ky) * chemoattr_tmp[i][j] + \
                              Kx * chemoattr_tmp[ ((i-1)+nblock)%nblock ][ j ] + \
                              Kx * chemoattr_tmp[ ((i+1)+nblock)%nblock ][ j ] + \
                              Ky * chemoattr_tmp[ i ][ ((j-1)+nblock)%nblock ] + \
                              Ky * chemoattr_tmp[ i ][ ((j+1)+nblock)%nblock ];
        }
    }
}
