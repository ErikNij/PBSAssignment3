#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

void record_trajectories_pdb(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
/*  Write the particle positions to a pdf file
    The filename (without extension) is given by p_parameters->filename_pdb.
    If reset = 1 the data is written to the file deleting data it possibly contained.
    If reset = 0 the data is appended. */
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_pdb, ".pdb");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "MODEL\n");
  fprintf(fp_traj, "REMARK TIME = %f\n", time);
  fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs*p_parameters->L.x, rs*p_parameters->L.y, rs*p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "HETATM %5u  C 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", (unsigned int)i % 100000, rs*p_vectors->r[i].x, rs*p_vectors->r[i].y, rs*p_vectors->r[i].z);
  }
  fprintf(fp_traj, "ENDMDL\n");

  fclose(fp_traj);
}

void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
/*  Write the particle positions to a xyz file
    The filename (without extension) is given by p_parameters->filename_xyz.
    If reset = 1 the data is written to the file deleting data it possibly contained.
    If reset = 0 the data is appended. */
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_xyz, ".dat");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "%lu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", time);
  struct Vec3D *r = p_vectors->r;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f\n", rs*r[i].x, rs*r[i].y, rs*r[i].z);
  }

  fclose(fp_traj);
}

void record_radial(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fprad;
  char filename[1024];

  snprintf(filename, 1024, "%s%s", p_parameters->filename_radial, ".dat");
  if (reset == 1)
  {
    fprad = fopen(filename, "w");
  }
  else
  {
    fprad = fopen(filename, "a");
  }

double Volume = p_parameters->L.x * p_parameters->L.y * p_parameters->L.z;
double frac = 1.333; // 4/3 in numbers
double Dens = p_parameters->num_part/Volume;


for(size_t numbin=0; numbin<100; numbin++)
{
  double numbin1 = numbin+1;
  double rad = (numbin1*numbin1*numbin1)-(numbin*numbin*numbin);
  p_parameters->rdf[numbin] = p_parameters->rdf[numbin]/(Dens * frac * PI * rad);
}


  fprintf(fprad, "MODEL\n");
  fprintf(fprad, "REMARK TIME = %f\n", time);
  fprintf(fprad, "CRYST1%9.3f\n", p_parameters->rdf, "P 1", "1");
  for (size_t i = 0; i < 100; i++)
  {
    fprintf(fprad, "HETATM %5u  C 1 UNK     1    %7.4f\n", (unsigned int)i % 100000, p_parameters->rdf[i]);
  }
  fprintf(fprad, "ENDMDL\n");

  fclose(fprad);
}


void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* save arrays in vectors to binary file */
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* load arrays in vectors from binary file */
{
  FILE* p_file = fopen( p_parameters->restart_in_filename, "rb" );
  size_t num_part;
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part*sizeof(struct Vec3D);
  alloc_vectors(p_vectors,num_part);
  p_parameters->num_part = num_part;
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}


