/*
 *      BEM_Poisson_Boltzmann.h
 *	Author: Robin Richardson, University of Leeds
 *	Email: pyrar@leeds.ac.uk
 */

#ifndef BEM_POISSON_BOLTZMANN_H_INCLUDED
#define BEM_POISSON_BOLTZMANN_H_INCLUDED

#include <math.h>

#include "FFEA_return_codes.h"
#include "NearestNeighbourLinkedListCube.h"
#include "SparseMatrixUnknownPattern.h"

#include "GaussianQuadrature_tri.h"
#include "GaussianQuadrature_1d.h"

class BEM_Poisson_Boltzmann: GaussianQuadrature_tri, GaussianQuadrature_1d
{
	public:
		BEM_Poisson_Boltzmann();
		~BEM_Poisson_Boltzmann();
		int init(NearestNeighbourLinkedListCube *lookup);
		/*
		 * Sets the inverse debye screening length for the system
		 */
		void set_kappa(scalar kappa);
		void build_BEM_matrices();
		void perform_integrals_for_lookup_cell_self(LinkedListNode<Face> *l_i, vector3 gqp[4]);
		void perform_integrals_for_lookup_cell_relative(LinkedListNode<Face> *l_i, vector3 gqp[4], int dx, int dy, int dz);
		void print_matrices();
		SparseMatrixUnknownPattern * get_C(); 
		SparseMatrixUnknownPattern * get_D();

	private:

		/* Nearest neighbour lookup data structure containing all faces in the system */
		NearestNeighbourLinkedListCube *lookup;

		/* Number of faces in system */
		int num_faces;

		/* BEM matrices */
		SparseMatrixUnknownPattern *mat_C, *mat_D;

		/* The inverse Debye-screening length, kappa */
		scalar kappa;

		/* Returns the value of the fundamental solution u multiplied by 4*pi */
		scalar u_4pi(scalar r);

		/* Returns the radial component of grad of u multiplied by 4*pi (all other components are zero) */
		scalar grad_u_4pi(scalar r, scalar r2);

/*
		scalar screened_R_theta(scalar r_perp_mag, scalar half_theta_max, scalar theta_bar, scalar xi);
*/

		void gauss_quadrature_4_point(vector3 gqp[4], vector3 *p, scalar *int_u, scalar *int_du, Face *f);

		scalar self_term(vector3 *n0, vector3 *n1, vector3 *n2, int precision);

		scalar f_1d(scalar r);

		scalar f_3d(vector3 *p, vector3 *q);
};

#endif
