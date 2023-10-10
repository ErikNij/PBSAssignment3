#ifndef INITIALISE_H_
#define INITIALISE_H_


/**
 * @brief Initialise_bond_connectivity defines what particles are connected by bonds. This function is called by initialise_structure
 * 
 * @param p_parameters parameters used for initialization
 * @param[out] p_vectors member bonds is used to store the bond information
 * @see initialise_structure
 */
void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialise_structure creates bonds, angles and dihedrals used for force computations. Computes 1-2, 1-3 and 1-4 connectivity used in the neighbor list for non-bonded force computation
 * 
 * @param p_parameters parameters used for initialization
 * @param[out] p_vectors members bonds, angles and dihedrals are used to store connectivity information
 * @param[out] p_neighbors members head12, pairs12, head13, pairs13, head14, pairs14 store 1-2, 1-3 and 1-4 pairs
 * @see initialise_bonds, initialise, is_connected_12, is_connected_13, is_connected_14
 */
void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Initialise also calls initialize_bond_connectivity, initialise_bonds, initialise_structure, initialise_positions and initialise_velocities
 * 
 * @param p_parameters parameters used for initialization
 * @param p_vectors used members: r, v
 * @param p_nbrlist used members: head12, pair12, head13, pair13, head14, pair14 
 * @see initialise_bonds, initialise_structure initialise_positions, initialise_velocities
 */
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, size_t *p_step, double *p_time);

/**
 * @brief Initialises positions on a cubic lattice
 * 
 * @param p_parameters used members: L
 * @param p_vectors used members: r
 */
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises velocities according to Maxwell-Boltzmann distribution
 * 
 * @param p_parameters used members: kT, mass
 * @param p_vectors used members: r
 */
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* INITIALISE_H_ */
