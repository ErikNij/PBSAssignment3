#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr;
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        dr.x = p_vectors->v[i].x * p_parameters->dt;
        dr.y = p_vectors->v[i].y * p_parameters->dt;
        dr.z = p_vectors->v[i].z * p_parameters->dt;
        p_vectors->dr[i] = dr;
        p_vectors->r[i].x += dr.x; //updating positions
        p_vectors->r[i].y += dr.y;
        p_vectors->r[i].z += dr.z;
        p_nbrlist->dr[i].x += dr.x; // update displacements with respect to neighbor list creation are updated
        p_nbrlist->dr[i].y += dr.y;
        p_nbrlist->dr[i].z += dr.z;
        p_nbrlist->dr[i].sq = (p_nbrlist->dr[i].x) * (p_nbrlist->dr[i].x) + (p_nbrlist->dr[i].y) * (p_nbrlist->dr[i].y) + (p_nbrlist->dr[i].z) * (p_nbrlist->dr[i].z);
    }
}

double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
/* Update velocities for half a time step using forces. This function returns the kinetic energy.
   v(t+dt/2) = v(t) + 0.5*F(t)/m*dt, or
   v(t+dt) = v(t+dt/2)+0.5*F(t+dt)/m*dt
*/
{
    double Ekin = 0.0;
    const double factor = 0.5 / p_parameters->mass * p_parameters->dt;
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x += factor * p_vectors->f[i].x;
        p_vectors->v[i].y += factor * p_vectors->f[i].y;
        p_vectors->v[i].z += factor * p_vectors->f[i].z;
        Ekin += p_vectors->v[i].x * p_vectors->v[i].x + p_vectors->v[i].y * p_vectors->v[i].y + p_vectors->v[i].z * p_vectors->v[i].z;
    }
    Ekin = 0.5 * Ekin * p_parameters->mass;
    return Ekin;
}

void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* Apply boundary conditions. In this of periodic BCs case particles are put back in the box */
{
    struct Vec3D invL;

    invL.x = 1.0 / p_parameters->L.x;
    invL.y = 1.0 / p_parameters->L.y;
    invL.z = 1.0 / p_parameters->L.z;

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->r[i].x -= p_parameters->L.x * floor(p_vectors->r[i].x * invL.x);
        p_vectors->r[i].y -= p_parameters->L.y * floor(p_vectors->r[i].y * invL.y);
        p_vectors->r[i].z -= p_parameters->L.z * floor(p_vectors->r[i].z * invL.z);
    }
}

void Radialcalculation(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double xdist, ydist, zdist;
    struct Vec3D rij;
    double disttot;
    size_t nbin;
    for(size_t i=0; i<p_parameters->num_part; i++)
    {
        for(size_t j=0; j<p_parameters->num_part; j++)
        {
            rij.x = p_vectors->r[i].x-p_vectors->r[j].x;
            rij.y = p_vectors->r[i].y-p_vectors->r[j].y;
            rij.z = p_vectors->r[i].z-p_vectors->r[j].z;

            rij.x = rij.x - p_parameters->L.x *floor((rij.x/p_parameters->L.x) +0.5);
            rij.y = rij.y - p_parameters->L.y *floor((rij.y/p_parameters->L.y) +0.5);
            rij.z = rij.z - p_parameters->L.z *floor((rij.z/p_parameters->L.z) +0.5);


            xdist = (p_vectors->r[i].x-p_vectors->r[j].x)*(p_vectors->r[i].x-p_vectors->r[j].x);
            ydist = (p_vectors->r[i].y-p_vectors->r[j].y)*(p_vectors->r[i].y-p_vectors->r[j].y);
            zdist = (p_vectors->r[i].z-p_vectors->r[j].z)*(p_vectors->r[i].z-p_vectors->r[j].z);

            disttot = sqrt(xdist+ydist+zdist);
            nbin = floor(disttot);

            if(nbin>0 && nbin<100) 
            {
                p_parameters->rdf[nbin] +=1;
            }
        }
    }

}

