#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

void initialise_types(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    for (size_t i = 0; i < p_parameters->num_part; i++)
        p_vectors->type[i] = 0; // Specify particle type
}

void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    size_t num_bonds = 0;
    struct Bond *bonds = (struct Bond *)malloc(num_bonds * sizeof(struct Bond));

    /* 
        Provide code to specify bonds as pairs of particles: bond[k].i and bond[k].j
    */

    p_vectors->num_bonds = num_bonds;
    p_vectors->bonds = bonds;
}

void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    initialise_bond_connectivity(p_parameters, p_vectors);
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    size_t num_part = p_parameters->num_part;
    size_t *cnt = (size_t *)calloc(num_part + 1, sizeof(size_t));
    size_t cnt_tot = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        ++cnt[bonds[i].i + 1];
        ++cnt[bonds[i].j + 1];
    }
    size_t *head12 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head12[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head12[i] = cnt[i] + head12[i - 1];
        cnt[i] = head12[i];
    }
    size_t *pairs12 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_bonds; ++i)
    {
        pairs12[cnt[bonds[i].i]++] = bonds[i].j;
        pairs12[cnt[bonds[i].j]++] = bonds[i].i;
    }
    size_t num_angles = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        size_t num_pairs = head12[i] - head12[i - 1];
        num_angles += (num_pairs * (num_pairs - 1)) / 2;
    }
    struct Angle *angles = (struct Angle *)malloc(num_angles * sizeof(struct Angle));
    size_t m = 0;
    for (size_t i = 0; i < num_part; ++i)
        for (size_t j = head12[i]; j <= head12[i + 1]; ++j)
            for (size_t k = j + 1; k < head12[i + 1]; ++k)
            {
                angles[m].i = pairs12[j];
                angles[m].j = i;
                angles[m].k = pairs12[k];
                ++m;
            }
    size_t num_dihdr = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        num_dihdr += (head12[bonds[i].i + 1] - head12[bonds[i].i] - 1) * (head12[bonds[i].j + 1] - head12[bonds[i].j] - 1);
    }
    struct Dihedral *dihedrals = (struct Dihedral *) malloc(num_dihdr * sizeof(struct Dihedral));
    size_t k = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        struct Dihedral dihdr;
        dihdr.j = bonds[i].i;
        dihdr.k = bonds[i].j;
        for (size_t j = head12[dihdr.j]; j < head12[dihdr.j + 1]; ++j)
        {
            if (pairs12[j] == dihdr.k)
                continue;
            dihdr.i = pairs12[j];
            for (size_t j = head12[dihdr.k]; j < head12[dihdr.k + 1]; ++j)
            {
                if (pairs12[j] != dihdr.j)
                {
                    dihdr.l = pairs12[j];
                    dihedrals[k++] = dihdr;
                }
            }
        }
    }
    if (p_parameters->exclude_13_nb)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;
        for (size_t i = 0; i < num_angles; ++i)
        {
            ++cnt[angles[i].i + 1];
            ++cnt[angles[i].k + 1];
        }
        size_t *head13 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
        head13[0] = 0;
        for (size_t i = 1; i <= num_part; ++i)
        {
            head13[i] = cnt[i] + head13[i - 1];
            cnt[i] = head13[i];
        }
        size_t *pairs13 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
        for (size_t i = 0; i < num_angles; ++i)
        {
            pairs13[cnt[angles[i].i]++] = angles[i].k;
            pairs13[cnt[angles[i].k]++] = angles[i].i;
        }
        p_nbrlist->head13 = head13;
        p_nbrlist->pairs13 = pairs13;
    }

    for (size_t i = 0; i <= num_part; ++i)
        cnt[i] = 0;
    for (size_t i = 0; i < num_dihdr; ++i)
    {
        ++cnt[dihedrals[i].i + 1];
        ++cnt[dihedrals[i].l + 1];
    }
    size_t *head14 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head14[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head14[i] = cnt[i] + head14[i - 1];
        cnt[i] = head14[i];
    }
    size_t *pairs14 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_dihdr; ++i)
    {
        pairs14[cnt[dihedrals[i].i]++] = dihedrals[i].l;
        pairs14[cnt[dihedrals[i].l]++] = dihedrals[i].i;
    }
    p_nbrlist->head14 = head14;
    p_nbrlist->pairs14 = pairs14;
    p_vectors->num_angles = num_angles;
    p_vectors->angles = angles;
    p_vectors->num_dihedrals = num_dihdr;
    p_vectors->dihedrals = dihedrals;
    if (p_parameters->exclude_12_nb)
    {
        p_nbrlist->head12 = head12;
        p_nbrlist->pairs12 = pairs12;
    }
    else
    {
        free(head12);
        free(pairs12);
    }
    free(cnt);
}

void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, size_t *p_step, double *p_time)
/* Initialize all the variables of an MD simulation. */
{
    initialise_structure(p_parameters, p_vectors, p_nbrlist);
    srand(13); // Positive integer as seed for random number generator
    initialise_positions(p_parameters, p_vectors);
    initialise_velocities(p_parameters, p_vectors);
    *p_step = 0;   // Initialize the step as zero
    *p_time = 0.0; // Initialize the time to zero
    return;
}

void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors)
/*  Initialize positions of all particles.
    Particles are initialized on a square lattice.
*/
{
    struct Vec3D dr;
    struct Index3D n;
    double dl;
    int ipart;

    dl = pow(p_parameters->L.x * p_parameters->L.y * p_parameters->L.z / ((double)p_parameters->num_part), 1.0 / 3.0);
    n.i = (int)ceil(p_parameters->L.x / dl);
    n.j = (int)ceil(p_parameters->L.y / dl);
    n.k = (int)ceil(p_parameters->L.z / dl);
    dr.x = p_parameters->L.x / (double)n.i;
    dr.y = p_parameters->L.y / (double)n.j;
    dr.z = p_parameters->L.z / (double)n.k;
    ipart = 0;
    for (size_t i = 0; i < n.i; ++i)
        for (size_t j = 0; j < n.j; ++j)
            for (size_t k = 0; k < n.k; ++k, ++ipart)
            {
                if (ipart >= p_parameters->num_part)
                    break;
                p_vectors->r[ipart].x = (i + 0.5) * dr.x;
                p_vectors->r[ipart].y = (j + 0.5) * dr.y;
                p_vectors->r[ipart].z = (k + 0.5) * dr.z;
                //      p_vectors->r[ipart].x = p_parameters->L.x*generate_uniform_random();
                //      p_vectors->r[ipart].y = p_parameters->L.y*generate_uniform_random();
                //      p_vectors->r[ipart].z = p_parameters->L.z*generate_uniform_random();
            }
}

void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
/*  Initialize the velocities of all particles.
    Initial velocities are sampled according to the Maxwell-Boltzmann distribution
*/
{
    double sqrtktm = sqrt(p_parameters->kT / p_parameters->mass);
    struct Vec3D sumv = (struct Vec3D){
        0.0};

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x = sqrtktm * gauss();
        p_vectors->v[i].y = sqrtktm * gauss();
        p_vectors->v[i].z = sqrtktm * gauss();
        sumv.x += p_vectors->v[i].x;
        sumv.y += p_vectors->v[i].y;
        sumv.z += p_vectors->v[i].z;
    }

    sumv.x /= ((double)(p_parameters->num_part)); /* remove average velocity */
    sumv.y /= ((double)(p_parameters->num_part)); /* so total momentum is zero */
    sumv.z /= ((double)(p_parameters->num_part));
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x -= sumv.x;
        p_vectors->v[i].y -= sumv.y;
        p_vectors->v[i].z -= sumv.z;
    }
}
