#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
  fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs * p_parameters->L.x, rs * p_parameters->L.y, rs * p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
  char typeLetter;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    if (i < p_parameters->num_part * p_parameters->moleFrac)
    {
      typeLetter = 'C';
    }
    else
    {
      typeLetter = 'H';
    }
    //printf("%c",typeLetter);
    fprintf(fp_traj, "HETATM %5u  %c 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", (unsigned int)i % 100000,typeLetter, rs * p_vectors->r[i].x, rs * p_vectors->r[i].y, rs * p_vectors->r[i].z);
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

  snprintf(filename, 1024, "%s%s", p_parameters->filename_xyz, ".xyz");
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
    fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f\n", rs * r[i].x, rs * r[i].y, rs * r[i].z);
  }

  fclose(fp_traj);
}

void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* save arrays in vectors to binary file */
{
  FILE *p_file = fopen(p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part * sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v, sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* load arrays in vectors from binary file */
{
  FILE *p_file = fopen(p_parameters->restart_in_filename, "rb");
  size_t num_part;
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part * sizeof(struct Vec3D);
  alloc_vectors(p_vectors, num_part);
  p_parameters->num_part = num_part;
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

void density_profile(struct Parameters *p_parameters, struct Vectors *p_vectors, int step, FILE **fp, double **abcxyz)
{
  double *ax = abcxyz[0];
  double *ay = abcxyz[1];
  double *az = abcxyz[2];
  double *bx = abcxyz[3];
  double *by = abcxyz[4];
  double *bz = abcxyz[5];

  int a = 0;
  int b = 0;

  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    double test = p_parameters->num_part * p_parameters->moleFrac;
    if (i < p_parameters->num_part * p_parameters->moleFrac)
    {
      ax[(int)floor(p_vectors->r[i].x * p_parameters->resolutionDensity)]++;
      ay[(int)floor(p_vectors->r[i].y * p_parameters->resolutionDensity)]++;
      az[(int)floor(p_vectors->r[i].z * p_parameters->resolutionDensity)]++;
      a++;
    }
    else
    {
      bx[(int)floor(p_vectors->r[i].x * p_parameters->resolutionDensity)]++;
      by[(int)floor(p_vectors->r[i].y * p_parameters->resolutionDensity)]++;
      bz[(int)floor(p_vectors->r[i].z * p_parameters->resolutionDensity)]++;
      b++;
    }
  }
  
  if (step == p_parameters->num_dt_steps)
  {
    for (int i = 0; i < ceil(p_parameters->L.x * p_parameters->resolutionDensity); i++)
    {
      fprintf(fp[0], "%f,", ax[i]);
      fprintf(fp[3], "%f,", bx[i]);
    }
    for (int i = 0; i < ceil(p_parameters->L.y * p_parameters->resolutionDensity); i++)
    {
      fprintf(fp[1], "%f,", ay[i]);
      fprintf(fp[4], "%f,", by[i]);
    }
    for (int i = 0; i < ceil(p_parameters->L.z * p_parameters->resolutionDensity); i++)
    {
      fprintf(fp[2], "%f,", az[i]);
      fprintf(fp[5], "%f,", bz[i]);
    }
    for (int i = 0; i < 6; i++)
    {
      fprintf(fp[i], "%d,%d \n", a, b);
    }
  }
}

