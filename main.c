/******************************************************************************/
/*                                                    			              */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 4.0, 18/9/2018    			              */
/******************************************************************************/

/*
 * For the 2023 PBS assigment the code needs to be extended
 *
 * -Implement a Berendsen thermostat in dynamics.c
 * -Implement bonds in initialise_bonds in file initialise.c
 * -Initialize vectors.type such that particles get the proper type 
 * -Implement the bonded and non-bonded force in force.c. (Make the forces type dependent)
 * -Change the particle position initialisation such that it takes into account bond lengths and angles
 * -Implement the needed on-the-fly data analysis
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"

/**
 * @brief main The main of the MD code. After initialization,
 * a velocity-Verlet scheme is executed for a specified number of time steps.
 *
 * @return int 0 if successful
 */
int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    size_t step;
    double Ekin, Epot, time;

    set_parameters(&parameters);
    alloc_memory(&parameters, &vectors, &nbrlist);
    if (parameters.load_restart == 1)
    {
        load_restart(&parameters, &vectors);
        initialise_structure(&parameters, &vectors, &nbrlist);
        step = 0;
        time = 0.0;
    }
    else
        initialise(&parameters, &vectors, &nbrlist, &step, &time);
    build_nbrlist(&parameters, &vectors, &nbrlist);
    Epot = calculate_forces(&parameters, &nbrlist, &vectors);
    record_trajectories_pdb(1, &parameters, &vectors, time);

    while (step < parameters.num_dt_steps) // start of the velocity-Verlet loop
    {
        step++;
        time += parameters.dt;
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);
        thermostat(&parameters, &vectors, Ekin);
        update_positions(&parameters, &nbrlist, &vectors);
        boundary_conditions(&parameters, &vectors);
        update_nbrlist(&parameters, &vectors, &nbrlist);
        Epot = calculate_forces(&parameters, &nbrlist, &vectors);
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);

        printf("Step %zu, Time %f, Epot %f, Ekin %f, Etot %f\n", step, time, Epot, Ekin, Epot + Ekin);
        if (step % parameters.num_dt_pdb == 0) record_trajectories_pdb(0, &parameters, &vectors, time);
        if (step % parameters.num_dt_restart == 0) save_restart(&parameters,&vectors); 
    }
    save_restart(&parameters, &vectors);
    free_memory(&vectors, &nbrlist);

    return 0;
}
