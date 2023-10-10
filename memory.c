#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "nbrlist.h"

void alloc_vectors(struct Vectors *p_vectors, size_t sz)
/* Allocate the arrays in 'vectors' needed to store information of all particles */
{
    p_vectors->size = sz;
    p_vectors->type = (int *)malloc(sz * sizeof(size_t));
    p_vectors->r = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->dr = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->v = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->f = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->bonds = NULL;
    p_vectors->num_bonds = 0;
    p_vectors->angles = NULL;
    p_vectors->num_angles = 0;
    p_vectors->dihedrals = NULL;
    p_vectors->num_dihedrals = 0;
}

void free_vectors(struct Vectors *p_vectors)
/* Free the arrays in 'vectors' */
{
    free(p_vectors->type);
    p_vectors->type = NULL;
    free(p_vectors->r);
    p_vectors->r = NULL;
    free(p_vectors->dr);
    p_vectors->dr = NULL;
    free(p_vectors->v);
    p_vectors->v = NULL;
    free(p_vectors->f);
    p_vectors->f = NULL;
    free(p_vectors->bonds);
    p_vectors->bonds = NULL;
    free(p_vectors->angles);
    p_vectors->angles = NULL;
    free(p_vectors->dihedrals);
    p_vectors->dihedrals = NULL;
    p_vectors->size = 0;
    p_vectors->num_bonds = 0;
    p_vectors->num_angles = 0;
    p_vectors->num_dihedrals = 0;
}

void alloc_memory(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Allocate all variables needed in the MD simulation */
{    
    alloc_vectors(p_vectors, p_parameters->num_part);
    alloc_nbrlist(p_parameters, p_nbrlist);
}

void free_memory(struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Free the memory allocated by alloc_memory */
{
    free_vectors(p_vectors);
    free_nbrlist(p_nbrlist);
}
