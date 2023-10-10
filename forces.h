#ifndef FORCES_H_
#define FORCES_H_

#include "structs.h"

/**
 * @brief Calculate total forces on paricles
 * @param p_parameters used members: num_part, L
 * @param p_nbrlist used members: num_nbrs, nbr
 * @param p_vectors used members: r, f
 * @return double potential energy
 */
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Calculate non-bonded forces on paricles
 * @param p_parameters used members: num_part
 * @param p_nbrlist used members: num_nbrs, nbr
 * @param[out] p_vectors used members: f
 * @return double potential energy
 */
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Calculate bond-stretch forces on paricles
 * @param p_parameters used members: num_part, L
 * @param p_vectors used members: r, f
 * @return double potential energy
 */
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculate angle-bend forces on paricles
 * @param p_parameters used members: num_part, L
 * @param p_vectors used members: r, f
 * @return double potential energy
 */
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculate dihedral-torsion forces on paricles
 * @param p_parameters used members: num_part, L
 * @param p_vectors used members: r, f
 * @return double potential energy
 */
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* FORCES_H_ */