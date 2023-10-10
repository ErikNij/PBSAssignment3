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

void thermostat(struct Parameters *p_parameters, struct Vectors *p_vectors, double Ekin)
/* Change velocities by thermostatting */
{
}