#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;
    for (size_t i = 0; i < num_part; i++)
        // initialize the forces to zero
        f[i] = (struct Vec3D){0.0, 0.0, 0.0}; /*initialize forces to zero*/

    double Epot = calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_angle(p_parameters, p_vectors);
    Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);
    return Epot;
}

double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0;
    struct Bond * bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi={0};
    for(size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x*floor(rij.x/L.x+0.5); //apply minimum image convenction for bonded particles
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y*floor(rij.y/L.y+0.5); 
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z*floor(rij.z/L.z+0.5);

        /*
            Here provide the force calculation
        */

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;
    }
    return Epot;
}

double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0;
    struct Angle * angles = p_vectors->angles;
    size_t num_angles = p_vectors->num_angles;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj;
    struct Vec3D fi={0}, fk={0};
    for(size_t q = 0; q < num_angles; ++q)
    {
        size_t i = angles[q].i;
        size_t j = angles[q].j;
        size_t k = angles[q].k;

        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x*floor(rij.x/L.x+0.5); //apply minimum image convenction for bonded particles
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y*floor(rij.y/L.y+0.5); 
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z*floor(rij.z/L.z+0.5);

        rkj.x = r[k].x - r[j].x;
        rkj.x = rkj.x - L.x*floor(rkj.x/L.x+0.5); //apply minimum image convenction for bonded particles
        rkj.y = r[k].y - r[j].y;
        rkj.y = rkj.y - L.y*floor(rkj.y/L.y+0.5); 
        rkj.z = r[k].z - r[j].z;
        rkj.z = rkj.z - L.z*floor(rkj.z/L.z+0.5);

        /*
            Here provide the force calculation
        */

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= (fi.x + fk.x);
        f[j].y -= (fi.y + fk.y);
        f[j].z -= (fi.z + fk.z);
        f[k].x += fk.x;
        f[k].y += fk.y;
        f[k].z += fk.z;
    }
    return Epot;
}

double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0;
    struct Dihedral * dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    for(size_t q = 0; q < num_dihedrals; ++q)
    {
        size_t i = dihedrals[q].i;
        size_t j = dihedrals[q].j;
        size_t k = dihedrals[q].k;
        size_t l = dihedrals[q].l;

    }
    return Epot;
}


double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
/* Compute non-bonded forces on particles using the pairs in a neighbor list.
This function returns the total potential energy of the system. */
{
    struct Vec3D df;
    double r_cutsq, sigmasq, sr2, sr6, sr12, fr, prefctr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;
    sigmasq = p_parameters->sigma * p_parameters->sigma;
    double epsilon = p_parameters->epsilon;

    double Epot = 0.0, Epot_cutoff;
    sr2 = sigmasq / r_cutsq;
    sr6 = sr2 * sr2 * sr2;
    sr12 = sr6 * sr6;
    Epot_cutoff = sr12 - sr6;

    for (size_t k = 0; k < num_nbrs; k++)
    {
        // for each pair in the neighbor list compute the pair forces
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        if (rij.sq < r_cutsq)
        // Compute forces if the distance is smaller than the cutoff distance
        {
            // pair forces are given by the LJ interaction
            sr2 = sigmasq / rij.sq;
            sr6 = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;
            Epot += 4.0 * epsilon * (sr12 - sr6 - Epot_cutoff);
            fr = 24.0 * epsilon * (2.0 * sr12 - sr6) / rij.sq; //force divided by distance
            df.x = fr * rij.x;
            df.y = fr * rij.y;
            df.z = fr * rij.z;
            f[i].x += df.x;
            f[i].y += df.y;
            f[i].z += df.z;
            f[j].x -= df.x;
            f[j].y -= df.y;
            f[j].z -= df.z;
        }
    }

    return Epot;
}
