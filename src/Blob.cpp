//
//  This file is part of the FFEA simulation package
//
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file.
//
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
//
//  To help us fund FFEA development, we humbly ask that you cite
//  the research papers on the package.
//

#include "Blob.h"

Blob::Blob() {

    /* Initialise everything to zero */
    blob_index = 0;
    conformation_index = 0;
    previous_conformation_index = 0;
    state_index = 0;
    previous_state_index = 0;
    num_nodes = 0;
    num_elements = 0;
    num_surface_faces = 0;
    num_beads = 0;
    num_l_ctf = 0;
    num_r_ctf = 0;
    num_sltotal_ctf = 0;
    num_slsets_ctf = 0;
    num_binding_sites = 0;
    num_surface_nodes = 0;
    num_interior_nodes = 0;
    num_surface_elements = 0;
    num_interior_elements = 0;
    ctf_l_nodes = NULL;
    ctf_r_nodes = NULL;
    ctf_sl_faces = NULL;
    ctf_slsurf_ndx = NULL;
    ctf_sl_faces = NULL;
    ctf_sl_surfsize = NULL;
    ctf_l_forces = NULL;
    ctf_r_forces = NULL;
    ctf_r_axis = NULL;
    ctf_r_type = NULL;
    ctf_sl_forces = NULL;
    mass = 0;
    blob_state = FFEA_BLOB_IS_STATIC;
    node = NULL;
    node_position = NULL;
    elem = NULL;
    surface = NULL;
    binding_site = NULL;
    solver = NULL;
    linear_solver = 0;
    mass_in_blob = false;
    ssint_on_blob = false;
    springs_on_blob = false;
    beads_on_blob = false;
    force = NULL;
    rng = NULL;
    poisson_solver = NULL;
    phi_Omega = NULL;
    phi_Gamma = NULL;
    q = NULL;
    nodal_q = NULL;
    poisson_rhs = NULL;
    num_pinned_nodes = 0;
    pinned_nodes_list = NULL;
    bsite_pinned_nodes_list.clear();

    num_contributing_faces = NULL;

    pbc_count[0]= 0;
    pbc_count[1]= 0;
    pbc_count[2]= 0;

    toBePrinted_nodes = NULL;

    poisson_surface_matrix = NULL;
    poisson_interior_matrix = NULL;
}

Blob::~Blob() {
    /* Release the node, element and surface arrays */
    delete[] node;
    node = NULL;
    delete[] node_position;
    node_position = NULL;
    delete[] elem;
    elem = NULL;
    delete[] surface;
    surface = NULL;

    delete[] binding_site;
    binding_site = NULL;

    /* Release the force vector */
    delete[] force;
    force = NULL;

    /* Release the Solver */
    delete solver;
    solver = NULL;
    linear_solver = 0;
    mass_in_blob = false;
    ssint_on_blob = false;
    springs_on_blob = false;
    beads_on_blob = false;

    /* Release the Poisson Solver */
    delete poisson_solver;
    poisson_solver = NULL;
    delete[] phi_Omega;
    phi_Omega = NULL;
    delete[] phi_Gamma;
    phi_Gamma = NULL;
    delete[] q;
    q = NULL;
    delete[] nodal_q;
    nodal_q = NULL;
    delete[] poisson_rhs;
    poisson_rhs = NULL;
    delete[] pinned_nodes_list;
    pinned_nodes_list = NULL;

    /* delete precomp stuff */
    if (num_beads > 0) {
        num_beads = 0;
        delete[] bead_position;
        bead_position = NULL;
        delete[] bead_type;
        bead_type = NULL;
    }

    /* delete ctforces stuff */
    //  first linear:
    if (num_l_ctf > 0) {
        num_l_ctf = 0;
        delete[] ctf_l_nodes;
        ctf_l_nodes = NULL;
        delete[] ctf_l_forces;
        ctf_l_forces = NULL;
    }
    // then rotational:
    if (num_r_ctf > 0) {
        num_r_ctf = 0;
        delete[] ctf_r_nodes;
        ctf_r_nodes = NULL;
        delete[] ctf_r_forces;
        ctf_r_forces = NULL;
        delete[] ctf_r_axis;
        ctf_r_axis = NULL;
        delete[] ctf_r_type;
        ctf_r_type = NULL;
    }
    // and then linear onto surfaces:
    if (num_sltotal_ctf > 0) {
        num_sltotal_ctf = 0;
        num_slsets_ctf = 0;
        delete[] ctf_slsurf_ndx;
        ctf_slsurf_ndx = NULL;
        delete[] ctf_sl_faces;
        ctf_sl_faces = NULL;
        delete[] ctf_sl_surfsize;
        ctf_sl_surfsize = NULL;
        delete[] ctf_sl_forces;
        ctf_sl_forces = NULL;
    }

    /* Set relevant data to zero */
    conformation_index = 0;
    previous_conformation_index = 0;
    state_index = 0;
    previous_state_index = 0;
    num_nodes = 0;
    num_elements = 0;
    num_surface_elements = 0;
    num_interior_elements = 0;
    num_surface_faces = 0;
    num_binding_sites = 0;
    num_surface_nodes = 0;
    num_interior_nodes = 0;
    blob_state = FFEA_BLOB_IS_STATIC;
    mass = 0;
    rng = NULL;
    num_pinned_nodes = 0;
    bsite_pinned_nodes_list.clear();

    delete[] toBePrinted_nodes;
    toBePrinted_nodes = NULL;

    delete[] num_contributing_faces;
    num_contributing_faces = NULL;

    delete poisson_surface_matrix;
    delete poisson_interior_matrix;
    poisson_surface_matrix = NULL;
    poisson_interior_matrix = NULL;
}


int Blob::config(const int blob_index, const int conformation_index, string node_filename,
             string topology_filename, string surface_filename, string material_params_filename,
             string stokes_filename, string ssint_filename, string pin_filename,
             string binding_filename, string beads_filename, scalar scale, scalar calc_compress,
             scalar compress, int linear_solver, int blob_state, SimulationParams *params,
             PreComp_params *pc_params, SSINT_matrix *ssint_matrix,
             BindingSite_matrix *binding_matrix, RngStream rng[]){

    // Which blob and conformation am i?
    this->blob_index = blob_index;
    this->conformation_index = conformation_index;

    // Input files:
    this->s_node_filename = node_filename;
    this->s_topology_filename = topology_filename;
    this->s_surface_filename = surface_filename;
    this->s_material_params_filename = material_params_filename;
    this->s_stokes_filename = stokes_filename;
    this->s_ssint_filename = ssint_filename;
    this->s_beads_filename = beads_filename;
    this->s_binding_filename = binding_filename;
    this->s_pin_filename = pin_filename;

    // scaling coordinates:
    this->scale = scale;

    // compressing:
    this->calc_compress = calc_compress;
    this->compress = compress;

    // precomputed configuration:
    this->pc_params = pc_params;

    // BindingSite_matrix:
    this->binding_matrix = binding_matrix;

    // lennard-jones interaction matrix:
    this->ssint_matrix = ssint_matrix;

    // Store the pointer to the simulation parameters class
    this->params = params;
    this->pc_params = pc_params;

    // Store the pointer to the random number generator array
    this->rng = rng;

    // Get the blob state
    this->blob_state = blob_state;

    // Need to know solver type
    this->linear_solver = linear_solver;


    return FFEA_OK;
}

int Blob::init(){

    //Load the node, topology, surface, materials and stokes parameter files.
    if (load_nodes(s_node_filename.c_str(), scale) == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Error when loading Blob nodes.\n");
    }

    if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        if (load_topology(s_topology_filename.c_str()) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading Blob topology.\n");
        }
        if (load_surface(s_surface_filename.c_str(), params) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading Blob surface.\n");
        }
    } else {
        if (load_surface_no_topology(s_surface_filename.c_str(), params) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading Blob surface (ignoring topology).\n");
        }
    }

    if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        if (load_material_params(s_material_params_filename.c_str()) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading Blob element material params.\n");
        }
        if (load_stokes_params(s_stokes_filename.c_str(), scale) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading Blob stokes parameter file.\n");
        }
    }

    if (params->calc_ssint == 1 || params->calc_steric == 1) {
        if (load_ssint(s_ssint_filename.c_str(), ssint_matrix->get_num_types(), params->ssint_type) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading ssint parameter file.\n")
        }
    }

    if (params->calc_preComp == 1) {
        if (strcmp(s_beads_filename.c_str(), "") == 0) {
            FFEA_CAUTION_MESSG("No beads file was assigned to Blob %d\n", blob_index);
        } else if (load_beads(s_beads_filename.c_str(), pc_params, scale) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading beads file.\n")
        }
    }

    if (params->calc_ctforces == 1) {
        if (load_ctforces(params->ctforces_fname) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading constant forces\n");
        }
    }

    // Kinetic binding sites are still a structural property
    if ( load_binding_sites() == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Error when loading binding sites file.\n")
    }

    // This can be defaulted by leaving it out of the script
    if (blob_state == FFEA_BLOB_IS_DYNAMIC && s_pin_filename != "") {
        if (load_pinned_nodes(s_pin_filename.c_str()) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when loading Blob pinned nodes file.\n");
        }
    } else {
        num_pinned_nodes = 0;
        pinned_nodes_list = NULL;
    }

    // Linearise all the elements
    for (int i = 0; i < num_elements; i++) {
        elem[i].linearise_element();
    }

    // Get the rest jacobian, rest volume etc. of this Blob and store it for later use
    calc_rest_state_info();

    // Calculate the connectivity of the mesh, giving each node a list of pointers to elements
    // it is a member of. The pointers will be to the exact memory location in that element's struct
    // in which can be found the contribution to the force on that node.
    if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        printf("\t\tCalculating node-element connectivity...");
        calculate_node_element_connectivity();
        printf("\t\tdone\n");

        // Run a check on parameters that are dependent upon the solver type
        if((params->calc_stokes == 0 && params->calc_noise == 0) && (params->calc_springs == 1 || params->calc_ctforces == 1)) {
            if(linear_solver == FFEA_NOMASS_CG_SOLVER) {
                FFEA_CAUTION_MESSG("For Blob %d, Conformation %d:\n\tUsing Springs / Constant forces in conjuction with the CG_nomass solver without calc_stokes or calc_noise may not converge, specially if there are very few ctforces or springs.\n\tIf you encounter this problem, consider setting either calc_stokes or calc_noise (you may try a low temperature) re-run the system :)\n", blob_index, conformation_index)
            }
        }

        // Create the chosen linear equation Solver for this Blob
	// Only initialise for dynamic blobs. Static don not need to be mechanically solved
	if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
		if (linear_solver == FFEA_DIRECT_SOLVER) {
		    solver = new(std::nothrow) SparseSubstitutionSolver();
		    mass_in_blob = true;
		} else if (linear_solver == FFEA_ITERATIVE_SOLVER) {
		    solver = new(std::nothrow) ConjugateGradientSolver();
		    mass_in_blob = true;
		} else if (linear_solver == FFEA_MASSLUMPED_SOLVER) {
		    solver = new(std::nothrow) MassLumpedSolver();
		    mass_in_blob = true;
		} else if (linear_solver == FFEA_NOMASS_CG_SOLVER) {
		    solver = new(std::nothrow) NoMassCGSolver();
		} else {
		    FFEA_ERROR_MESSG("Error in Blob initialisation: linear_solver=%d is not a valid solver choice\n", linear_solver);
		}
		if (solver == NULL) FFEA_ERROR_MESSG("No solver to work with\n");
	}

        // Initialise the Solver (whatever it may be)
        printf("\t\tBuilding solver:\n");
        if (solver->init(num_nodes, num_elements, node, elem, params, num_pinned_nodes, pinned_nodes_list, bsite_pinned_nodes_list) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error initialising solver.\n")
        }
    }


    // Allocate the force vector array for the whole Blob
    force = new(std::nothrow) vector3[num_nodes];
    if (force == NULL) FFEA_ERROR_MESSG("Could not store the force vector array\n");
    for (int i = 0; i < num_nodes; i++) {
        force[i].x = 0;
        force[i].y = 0;
        force[i].z = 0;
    }

    // Calculate how many faces each surface node is a part of
    num_contributing_faces = new(std::nothrow) int[num_surface_nodes];

    if (num_contributing_faces == NULL) FFEA_ERROR_MESSG("Failed to allocate num_contributing_faces\n");
    for (int i = 0; i < num_surface_nodes; i++) {
        num_contributing_faces[i] = 0;
    }
    for (int i = 0; i < num_surface_faces; i++) {
        for (int j = 0; j < 3; j++) {
            num_contributing_faces[surface[i].n[j]->index]++;
        }
    }

    // Store stokes drag on nodes, for use in viscosity matrix
    if (params->calc_stokes == 1) {
        for (int i = 0; i < num_nodes; ++i) {
            node[i].stokes_drag = 6.0 * ffea_const::pi * params->stokes_visc * node[i].stokes_radius;
        }
    }

    // Generate the Mass matrix for this Blob
    if (params->calc_es == 1) {
        build_mass_matrix();
    }

    // Create and initialise the poisson solver
    if (num_interior_nodes > 0) {

        // Calculate the Sparsity Pattern for the Poisson matrix (interior and exterior)
        printf("\t\tCalculating sparsity pattern for Poisson matrix and RHS 'knowns' matrix\n");
        SparsityPattern sparsity_pattern_knowns, sparsity_pattern_unknowns;
        sparsity_pattern_knowns.init(num_interior_nodes);
        sparsity_pattern_unknowns.init(num_interior_nodes);

        scalar *mem_loc;
        int ni_index, nj_index;
        for (int el = 0; el < num_elements; el++) {
            for (int ni = 0; ni < 10; ni++) {
                for (int nj = 0; nj < 10; nj++) {

                    ni_index = elem[el].n[ni]->index;
                    nj_index = elem[el].n[nj]->index;

                    mem_loc = elem[el].get_K_alpha_element_mem_loc(ni, nj);

                    /* We don't care about rows of the matrix before row num_surface_nodes */
                    if (ni_index >= num_surface_nodes) {

                        /* if the column is in the surface node area, add these contributions to the 'known' matrix */
                        if (nj_index < num_surface_nodes) {
                            sparsity_pattern_knowns.register_contribution(
                                ni_index - num_surface_nodes,
                                nj_index,
                                mem_loc);
                        }
                        /* otherwise add the contribution to the 'unknowns' matrix */
                        else {
                            sparsity_pattern_unknowns.register_contribution(
                                ni_index - num_surface_nodes,
                                nj_index - num_surface_nodes,
                                mem_loc);
                        }
                    }
                }
            }
        }

        // Use the sparsity patterns to create fixed pattern sparse matrices for use with the poisson solver
        poisson_surface_matrix = sparsity_pattern_knowns.create_sparse_matrix();
        poisson_interior_matrix = sparsity_pattern_unknowns.create_sparse_matrix();

        // Create a conjugate gradient solver for use with the 'unknowns' (interior) poisson matrix
        printf("\t\tCreating and initialising Poisson Solver...");
        poisson_solver = new CG_solver();
        if (poisson_solver->init(num_interior_nodes, params->epsilon2, params->max_iterations_cg) != FFEA_OK) {
            FFEA_ERROR_MESSG("Failed to initialise poisson_solver\n");
        }

        // Create the vector containing all values of the potential at each interior node
        phi_Omega = new(std::nothrow) scalar[num_interior_nodes];

        if (phi_Omega == NULL) {
            FFEA_ERROR_MESSG("Could not allocate memory (for phi_Omega array).\n");
        }

        for (int i = 0; i < num_interior_nodes; i++) {
            phi_Omega[i] = 0;
        }

        // Create the vector containing all values of the potential on each surface node
        phi_Gamma = new(std::nothrow) scalar[num_surface_nodes];

        if (phi_Gamma == NULL) {
            FFEA_ERROR_MESSG("Could not allocate memory (for phi_Gamma array).\n");
        }

        for (int i = 0; i < num_surface_nodes; i++) {
            phi_Gamma[i] = 0;
        }


        /*
                        // Create the vector containing the smoothed out charge distribution across the elements
                        q = new scalar[num_nodes];

                        if(q == NULL) {
                                FFEA_error_text();
                                printf("Could not allocate memory (for q array).\n");
                                return FFEA_ERROR;
                        }

                        for(int i = 0; i < num_nodes; i++) {
                                q[i] = 0;
                        }

                        // Create the vector containing all charges on each node of the Blob
                        nodal_q = new scalar[num_nodes];

                        if(nodal_q == NULL) {
                                FFEA_error_text();
                                printf("Could not allocate memory (for nodal_q array).\n");
                                return FFEA_ERROR;
                        }

                        for(int i = 0; i < num_nodes; i++) {
                                nodal_q[i] = 0;
                        }

                        // Set the charges on the Blob
                        for(int i = 0; i < num_nodes; i++) {
        //			nodal_q[i] = (node[i].pos.z - (2.9e-8/2.0)) * 1e35; //1e26;
                                nodal_q[i] = 1e26;
        //			printf("%e %e\n", node[i].pos.z, nodal_q[i]);
                        }

                        // Apply mass matrix to obtain smooth extrapolated charge density across elements
                        M->apply(nodal_q, q);

                        for(int i = 0; i < num_nodes; i++) {
                                node[i].rho = q[i];
                        }
         */

        // Create the vector containing the charge distribution across the elements
        q = new(std::nothrow) scalar[num_nodes];

        if (q == NULL) {
            FFEA_ERROR_MESSG("Could not allocate memory (for q array).\n");
        }

        for (int i = 0; i < num_nodes; i++) {
            q[i] = 0;
        }

        // const scalar charge_density = 6.0e25;
        const scalar charge_density = 6.0e25 * mesoDimensions::volume / mesoDimensions::charge;
        for (int n = 0; n < num_elements; n++) {
            for (int i = 0; i < 10; i++) {
                q[elem[n].n[i]->index] += charge_density * elem[n].vol_0;
            }
        }

        for (int i = 0; i < num_nodes; i++) {
            node[i].rho = q[i];
        }

        //		for(int i = 0; i < num_nodes; i++) {
        //			printf("q[%d] = %e\n", i, q[i]);
        //		}

        // Create the Right Hand Side vector for the Poisson solver
        poisson_rhs = new(std::nothrow) scalar[num_interior_nodes];

        if (poisson_rhs == NULL) {
            FFEA_ERROR_MESSG("Could not allocate memory (for poisson_rhs array).\n");
        }

        for (int i = 0; i < num_interior_nodes; i++) {
            poisson_rhs[i] = 0;
        }

        //		build_poisson_matrices_and_setup_for_solve();
        //		poisson_solver->solve(poisson_interior_matrix, phi_Omega, poisson_rhs);

        printf("\t\tdone.\n");
    }


    //If calc_compress option set in input script (default off), compress by factor specified in script
    if (calc_compress ==1) {
        compress_blob(compress);
    }

    if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
        toBePrinted_nodes = new scalar[10*num_nodes];
    } else {
        if (params->calc_es == 0) {
            toBePrinted_nodes = new scalar[3*num_nodes];
        } else {
            toBePrinted_nodes = new scalar[4*num_nodes];
        }
    }


    // Return FFEA_OK to indicate "success"
    return FFEA_OK;
}

int Blob::check_inversion() {

	int n;
	matrix3 J;

	vector<int> invEls;
	invEls.clear();

        for (n = 0; n < num_elements; n++) {

            // calculate jacobian for this element
            elem[n].calculate_jacobian(J);

            // get the 12 derivatives of the shape functions (by inverting the jacobian)
            // and also get the element volume. The function returns an error in the
            // case of an element inverting itself (determinant changing sign since last step)
            if (elem[n].calc_shape_function_derivatives_and_volume(J) == FFEA_ERROR) {
	        invEls.push_back(n);
            }
	}

	// Are we screwed?
	if(invEls.size() != 0) {
		printf("\n");
		FFEA_error_text();
		printf("%zu inverted elements: ", invEls.size());
		for(n = 0; n < invEls.size(); ++n) {
	                printf("%d ", invEls.at(n));
		}
		printf("\n");
		return FFEA_ERROR;
	}

	return FFEA_OK;
}

int Blob::update_internal_forces() {
    if (blob_state != FFEA_BLOB_IS_DYNAMIC) {
        return FFEA_OK;
    }

    if(num_pinned_nodes == num_nodes) {
        return FFEA_OK;
    }

    if (params->calc_ctforces) apply_ctforces();



    /* some "work" variables */
    matrix3 J; // Holds the Jacobian calculated for the *current* element being processed
    matrix3 stress; // Holds the current stress tensor (elastic stress, with thermal fluctuations)
    vector12 du; // Holds the force change for the current element
    int n; // Iterates through all the elements
    int tid; // Holds the current thread id (in parallel regions)
    int num_inversions = 0; // Counts the number of elements that have inverted (if > 0 then simulation has failed)

    // Element loop
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel default(none) private(J, stress, du, tid, n) reduction(+:num_inversions)
    {
#endif
#ifdef USE_OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif

#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp for schedule(guided)
#endif
        for (n = 0; n < num_elements; n++) {

            // calculate jacobian for this element
            elem[n].calculate_jacobian(J);

            // get the 12 derivatives of the shape functions (by inverting the jacobian)
            // and also get the element volume. The function returns an error in the
            // case of an element inverting itself (determinant changing sign since last step)
            if (elem[n].calc_shape_function_derivatives_and_volume(J) == FFEA_ERROR) {
                FFEA_error_text();
                printf("Element %d has inverted during update\n", n);
           	//elem[n].print_structural_details();
	        num_inversions++;
            }

            // create viscosity matrix
            elem[n].create_viscosity_matrix();

            // Now build the stress tensor from the shear elastic, bulk elastic and fluctuating stress contributions
            mat3_set_zero(stress);
            elem[n].add_shear_elastic_stress(J, stress);
            elem[n].add_bulk_elastic_stress(stress);

            if (params->calc_noise == 1) {
                elem[n].add_fluctuating_stress(params, rng, stress, tid);
            }

            elem[n].internal_stress_mag = sqrt(mat3_double_contraction_symmetric(stress));

            // Calculate internal forces of current element (or don't, depending on solver)
            if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
                elem[n].get_element_velocity_vector(du);
                mat12_apply(elem[n].viscosity_matrix, du);
            } else {
                vec12_set_zero(du);
            }

            elem[n].apply_stress_tensor(stress, du);

            // Store the contributions to the force on each of this element's nodes (Store them on
            // the element - they will be aggregated on the actual nodes outside of this parallel region)
            elem[n].add_element_force_vector(du);

            if (params->calc_es == 1) {
                elem[n].calculate_electrostatic_forces();
            }
        }

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    }
#endif

    // Check if any elements have inverted
    if (num_inversions != 0) {
        if (num_inversions == 1) {
            FFEA_ERROR_MESSG("1 element has inverted since the last step. Aborting simulation.\n");
        } else {
            FFEA_ERROR_MESSG("%d elements have inverted since the last step. Aborting simulation.\n", num_inversions);
        }
    }

    return FFEA_OK;
}

int Blob::update_positions() {

    // Aggregate forces on nodes from all elements (if not static)
    if (get_motion_state() != FFEA_BLOB_IS_DYNAMIC) {
	return FFEA_OK;
    }

    if (aggregate_forces_and_solve() == FFEA_ERROR) {
        FFEA_ERROR_MESSG("There was a problem in 'aggregate_forces_and_solve' function.\n");
    }

    // Update node velocities and positions
    euler_integrate();

    // todo: remove this
    //for (int i=0; i<this->get_num_elements(); i++){
    //    for (int j=0; j<9; j++){
    //        printf("%f", this->get_element(i)->n[j]->pos.data[0]);
    //    }
    //}

    // Linearise the 2nd order elements
    for (int n = 0; n < num_elements; n++) {
        elem[n].linearise_element();
    }

    return FFEA_OK;
}

int Blob::reset_solver() {

    // Delete and rebuild (to make sure everything is overwritten)

    if (solver->init(num_nodes, num_elements, node, elem, params, num_pinned_nodes, pinned_nodes_list, bsite_pinned_nodes_list) == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Error reinitialising solver.\n")
    } else {
        return FFEA_OK;
    }
    params->calc_noise = 0;
}

void Blob::translate_linear(vector3 *vec) {

    // Get a mapping from all node indices to just linear node indices
    int num_linear_nodes = get_num_linear_nodes();
    int map[num_linear_nodes];
    int i, j;
    j = 0;
    for(i = 0; i < num_nodes; ++i) {
        if(node[i].am_I_linear()) {
            map[j] = i;
            j++;
        }
    }

    // Translate linear nodes
    for(i = 0; i < num_linear_nodes; ++i) {
        node[map[i]].pos.x += vec[i].x;
        node[map[i]].pos.y += vec[i].y;
        node[map[i]].pos.z += vec[i].z;
    }

    // Sort secondary nodes
    linearise_elements();
}

// Rotate about x axis, then y axis, then z axis
void Blob::rotate(float xang, float yang, float zang, bool beads) {

    scalar r[3][3];

    // Convert to radians
    xang *= ffea_const::pi / 180.0;
    yang *= ffea_const::pi / 180.0;
    zang *= ffea_const::pi / 180.0;

    // Get rotation
    r[0][0] = cos(yang) * cos(zang);
    r[0][1] = sin(xang) * sin(yang) * cos(zang) - cos(xang) * sin(zang);
    r[0][2] = cos(xang) * sin(yang) * cos(zang) + sin(xang) * sin(zang);
    r[1][0] = cos(yang) * sin(zang);
    r[1][1] = sin(xang) * sin(yang) * sin(zang) + cos(xang) * cos(zang);
    r[1][2] = cos(xang) * sin(yang) * sin(zang) - sin(xang) * cos(zang);
    r[2][0] = -1 * sin(yang);
    r[2][1] = sin(xang) * cos(yang);
    r[2][2] = cos(xang) * cos(yang);

    rotate(r[0][0], r[0][1], r[0][2],
           r[1][0], r[1][1], r[1][2],
           r[2][0], r[2][1], r[2][2], beads);

}

void Blob::rotate(float r11, float r12, float r13, float r21, float r22, float r23, float r31, float r32, float r33, bool beads) {
    int i;
    vector3 com;
    scalar x, y, z;

    get_centroid(&com);

    // Move all nodes to the origin:
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) private(i) shared(com)
#endif
    for (i = 0; i < num_nodes; i++) {
        node[i].pos.x -= com.x;
        node[i].pos.y -= com.y;
        node[i].pos.z -= com.z;
    }

    // Do the actual rotation and bring the nodes back to its initial position:
    for (i = 0; i < num_nodes; i++) {
        x = node[i].pos.x;
        y = node[i].pos.y;
        z = node[i].pos.z;

        node[i].pos.x = x * r11 + y * r12 + z * r13 + com.x;
        node[i].pos.y = x * r21 + y * r22 + z * r23 + com.y;
        node[i].pos.z = x * r31 + y * r32 + z * r33 + com.z;


    }


    if (beads) {
        if (num_beads > 0) {
            // Move all beads to the origin:
            for (i = 0; i < num_beads; i++) {
                bead_position[3*i] -= com.x;
                bead_position[3*i+1] -= com.y;
                bead_position[3*i+2] -= com.z;
            }

            // Do the actual rotation and bring the beads back to its initial position:
            for (i = 0; i < num_beads; i++) {
                x = bead_position[3*i];
                y = bead_position[3*i+1];
                z = bead_position[3*i+2];

                bead_position[3*i  ] = x * r11 + y * r12 + z * r13 + com.x;
                bead_position[3*i+1] = x * r21 + y * r22 + z * r23 + com.y;
                bead_position[3*i+2] = x * r31 + y * r32 + z * r33 + com.z;

            }
        }
    }
}


vector3 Blob::position(scalar x, scalar y, scalar z) {
    int i;
    scalar centroid_x = 0.0, centroid_y = 0.0, centroid_z = 0.0;
    vector3 v;
    // scalar dx, dy, dz;

    // Calculate centroid of entire Blob mesh
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(shared) private(i) reduction(+:centroid_x,centroid_y,centroid_z)
#endif
    for (i = 0; i < num_nodes; i++) {
        centroid_x += node[i].pos.x;
        centroid_y += node[i].pos.y;
        centroid_z += node[i].pos.z;
    }

    centroid_x *= (1.0 / num_nodes);
    centroid_y *= (1.0 / num_nodes);
    centroid_z *= (1.0 / num_nodes);

    // Calculate displacement vector required to move centroid to requested position
    v.x = x - centroid_x;
    v.y = y - centroid_y;
    v.z = z - centroid_z;

    // Move all nodes in mesh by displacement vector
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) private(i) shared(v)
#endif
    for (i = 0; i < num_nodes; i++) {
        node[i].pos.x += v.x;
        node[i].pos.y += v.y;
        node[i].pos.z += v.z;
    }

    return v;

}

void Blob::position_beads(scalar x, scalar y, scalar z) {

    for (int i = 0; i < num_beads; i ++) {
        bead_position[3*i] += x;
        bead_position[3*i+1] += y;
        bead_position[3*i+2] += z;
    }

}

void Blob::move(scalar dx, scalar dy, scalar dz) {
    for (int i = 0; i < num_nodes; i++) {
        node[i].pos.x += dx;
        node[i].pos.y += dy;
        node[i].pos.z += dz;
    }
    for(int i = 0; i < get_num_faces(); ++i) {
        if(get_motion_state() != FFEA_BLOB_IS_DYNAMIC && surface[i].n[3] != NULL) {
            surface[i].n[3]->pos.x += dx;
            surface[i].n[3]->pos.y += dy;
            surface[i].n[3]->pos.z += dz;
            //fprintf(stderr, "Surface %d = %f %f %f\n", i, dx, dy, dz);
        }
    }
}

void Blob::get_CoM(vector3 *com) {
    com->x = 0;
    com->y = 0;
    com->z = 0;


    for (int n = 0; n < num_elements; n++) {
        com->x += (elem[n].mass * elem[n].n[0]->pos.x + elem[n].n[1]->pos.x + elem[n].n[2]->pos.x + elem[n].n[3]->pos.x) / 4;
        com->y += (elem[n].mass * elem[n].n[0]->pos.y + elem[n].n[1]->pos.y + elem[n].n[2]->pos.y + elem[n].n[3]->pos.y) / 4;
        com->z += (elem[n].mass * elem[n].n[0]->pos.z + elem[n].n[1]->pos.z + elem[n].n[2]->pos.z + elem[n].n[3]->pos.z) / 4;
    }
    if (num_elements == 0) {
        return;
    } else {
        com->x /= num_elements;
        com->y /= num_elements;
        com->z /= num_elements;
    }
}

void Blob::get_centroid(vector3 *com) {
    com->x = 0;
    com->y = 0;
    com->z = 0;
    for (int n = 0; n < num_nodes; n++) {
        com->x += node[n].pos.x;
        com->y += node[n].pos.y;
        com->z += node[n].pos.z;
    }
    com->x /= num_nodes;
    com->y /= num_nodes;
    com->z /= num_nodes;
}

void Blob::calc_and_store_centroid(vector3 &com) {

    get_centroid(&com);

    CoG[0] = com[0];
    CoG[1] = com[1];
    CoG[2] = com[2];

}

// This one returns an array rather than arsing about with pointers
vector3 Blob::calc_centroid() {

    vector3 com;
    com.x = 0.0;
    com.y = 0.0;
    com.z = 0.0;
    for (int n = 0; n < num_nodes; n++) {
        com.x += node[n].pos.x;
        com.y += node[n].pos.y;
        com.z += node[n].pos.z;
    }
    com.x /= num_nodes;
    com.y /= num_nodes;
    com.z /= num_nodes;
    return com;
}

vector3 ** Blob::get_actual_node_positions() {
    return node_position;
}

void Blob::copy_node_positions(vector3 *nodes) {

    for(int i = 0; i < num_nodes; ++i) {
        nodes[i].x = node[i].pos.x;
        nodes[i].y = node[i].pos.y;
        nodes[i].z = node[i].pos.z;
    }
}

void Blob::set_node_positions(vector3 *node_pos) {

    for (int i = 0; i < num_nodes; i++) {
        node[i].pos.x = node_pos[i].x;
        node[i].pos.y = node_pos[i].y;
        node[i].pos.z = node_pos[i].z;
    }
}

void Blob::set_pos_0() {

    CoG_0.x = 0.0;
    CoG_0.y = 0.0;
    CoG_0.z = 0.0;
    for (int i = 0; i < num_nodes; i++) {
        node[i].pos_0.x = node[i].pos.x;
        node[i].pos_0.y = node[i].pos.y;
        node[i].pos_0.z = node[i].pos.z;
        CoG_0.x += node[i].pos_0.x;
        CoG_0.y += node[i].pos_0.y;
        CoG_0.z += node[i].pos_0.z;
    }
    CoG_0.x /= num_nodes;
    CoG_0.y /= num_nodes;
    CoG_0.z /= num_nodes;

}

void Blob::kinetically_set_faces(bool state) {
    for(int i = 0; i < num_surface_faces; ++i) {
        surface[i].set_kinetic_state(state);
    }
}

void Blob::linearise_elements() {
    for(int i = 0; i < num_elements; ++i) {
        elem[i].linearise_element();
    }
}

void Blob::linearise_force() {

    int nIdx[NUM_NODES_QUADRATIC_TET];
    for (int i=0; i< num_elements; i++) {
        for (int j=0; j<NUM_NODES_QUADRATIC_TET; j++) {
            nIdx[j] = elem[i].n[j]->index;
        }
        for (int j=0; j<3; j++){
           force[nIdx[0]][j] += 0.5 * ( force[nIdx[4]][j] + force[nIdx[5]][j] + force[nIdx[6]][j]);
           force[nIdx[1]][j] += 0.5 * ( force[nIdx[4]][j] + force[nIdx[7]][j] + force[nIdx[8]][j]);
           force[nIdx[2]][j] += 0.5 * ( force[nIdx[5]][j] + force[nIdx[7]][j] + force[nIdx[9]][j]);
           force[nIdx[3]][j] += 0.5 * ( force[nIdx[6]][j] + force[nIdx[8]][j] + force[nIdx[9]][j]);
        }
        for (int j=4; j<NUM_NODES_QUADRATIC_TET; j++) {
            force[nIdx[j]][0]   = 0.;
            force[nIdx[j]][1]   = 0.;
            force[nIdx[j]][2]   = 0.;
        }
    }
}

void Blob::compress_blob(scalar compress) {


    vector3 cog,cogaft;
    this->get_centroid(&cog);
    //loop moves nodes in by
    for (int i = 0; i < num_nodes; i++) {
        node[i].pos.x *= compress;
        node[i].pos.y *= compress;
        node[i].pos.z *= compress;
    }
    this->get_centroid(&cogaft);

    if (!(fabs(cogaft.x - cog.x*compress) < 0.000001&&fabs(cogaft.y - cog.y*compress) < 0.000001&&fabs(cogaft.z - cog.z*compress) < 0.000001)) {
        printf("FRIENDLY WARNING: Centre of Geometry has moved during compression. Compression feature implemented for spherical objects so other shapes may experience issues.\n");
    }
}

int Blob::create_viewer_node_file(const char *node_filename, scalar scale) {

    FILE *out = NULL;
    char new_node_filename[] = "VIEWERNODE_";
    int i;

    // Name new node file
    strcat(new_node_filename, node_filename);


    //open the new node file
    if ((out = fopen(new_node_filename, "w")) == NULL) {
        FFEA_FILE_ERROR_MESSG(new_node_filename);
    }
    printf("\t\tWriting to viewer nodes file: %s\n", new_node_filename);

    fprintf(out, "ffea viewer node file\n");
    fprintf(out, "num_nodes %d\n", num_nodes);
    fprintf(out, "num_surface_nodes %d\n", num_surface_nodes);
    fprintf(out, "num_interior_nodes %d\n", num_interior_nodes);

    // Write all the nodes to file
    fprintf(out, "surface nodes:\n");
    for (i = 0; i < num_surface_nodes; i++) {
        fprintf(out, "%le %le %le\n", node[i].pos.x / scale, node[i].pos.y / scale, node[i].pos.z / scale);
    }

    fprintf(out, "interior nodes:\n");
    for (i = num_surface_nodes; i < num_nodes; i++) {
        fprintf(out, "%le %le %le\n", node[i].pos.x / scale, node[i].pos.y / scale, node[i].pos.z / scale);
    }

    fclose(out);
    printf("\t\t\tWrote %d nodes from %s\n", i, new_node_filename);

    return FFEA_OK;
}

void Blob::write_nodes_to_file(FILE *trajectory_out) {
    // If this is a static blob, then don't bother printing out all the node positions (since there will be no change)
    if (blob_state == FFEA_BLOB_IS_STATIC) {
        fprintf(trajectory_out, "STATIC\n");
        return;
    } else if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        fprintf(trajectory_out, "DYNAMIC\n");
    } else if (blob_state == FFEA_BLOB_IS_FROZEN) {
        fprintf(trajectory_out, "FROZEN\n");
    }

    if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
        for (int i = 0; i < num_nodes; i++) {
            fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                    node[i].pos.x*mesoDimensions::length, node[i].pos.y*mesoDimensions::length, node[i].pos.z*mesoDimensions::length,
                    node[i].vel.x*mesoDimensions::velocity, node[i].vel.y*mesoDimensions::velocity, node[i].vel.z*mesoDimensions::velocity,
                    node[i].phi,
                    force[i].x*mesoDimensions::force, force[i].y*mesoDimensions::force, force[i].z*mesoDimensions::force);
        }
    } else {
        if (params->calc_es == 0) {
            for (int i = 0; i < num_nodes; i++) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                        node[i].pos.x*mesoDimensions::length, node[i].pos.y*mesoDimensions::length, node[i].pos.z*mesoDimensions::length,
                        0., 0., 0., 0., 0., 0., 0.);
            }
        } else {
            for (int i = 0; i < num_nodes; i++) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                        node[i].pos.x*mesoDimensions::length, node[i].pos.y*mesoDimensions::length, node[i].pos.z*mesoDimensions::length,
                        0., 0., 0.,
                        node[i].phi, 0., 0., 0.);
            }
        }
    }
}

void Blob::write_pre_print_to_file(FILE *trajectory_out) {
    // If this is a static blob, then don't bother printing out all the node positions (since there will be no change)
    if (blob_state == FFEA_BLOB_IS_STATIC) {
        fprintf(trajectory_out, "STATIC\n");
        return;
    } else if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        fprintf(trajectory_out, "DYNAMIC\n");
    } else if (blob_state == FFEA_BLOB_IS_FROZEN) {
        fprintf(trajectory_out, "FROZEN\n");
    }

    if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
        for (int i = 0; i < num_nodes; i++) {
            fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                    toBePrinted_nodes[10*i], toBePrinted_nodes[10*i+1], toBePrinted_nodes[10*i+2],
                    toBePrinted_nodes[10*i+3], toBePrinted_nodes[10*i+4], toBePrinted_nodes[10*i+5],
                    toBePrinted_nodes[10*i+6], toBePrinted_nodes[10*i+7], toBePrinted_nodes[10*i+8],
                    toBePrinted_nodes[10*i+9]);
        }
    } else {
        if (params->calc_es == 0) {
            for (int i = 0; i < num_nodes; i++) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                        toBePrinted_nodes[3*i], toBePrinted_nodes[3*i+1], toBePrinted_nodes[3*i+2],
                        0., 0., 0., 0., 0., 0., 0.);
            }
        } else {
            for (int i = 0; i < num_nodes; i++) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                        toBePrinted_nodes[3*i], toBePrinted_nodes[3*i+1], toBePrinted_nodes[3*i+2],
                        0., 0., 0.,
                        toBePrinted_nodes[3*i+3], 0., 0., 0.);
            }
        }
    }
}

void Blob::pre_print() {
    // If this is a static blob, then don't bother printing out all the node positions (since there will be no change)
    if (blob_state == FFEA_BLOB_IS_STATIC) {
        return;
    }

    toBePrinted_conf[0] = previous_conformation_index;
    toBePrinted_conf[1] = conformation_index;
    toBePrinted_state[0] = previous_state_index;
    toBePrinted_state[1] = state_index;

    if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none)
#endif
        for (int i = 0; i < num_nodes; i++) {
            toBePrinted_nodes[10*i   ] = node[i].pos.x*mesoDimensions::length;
            toBePrinted_nodes[10*i +1] = node[i].pos.y*mesoDimensions::length;
            toBePrinted_nodes[10*i +2] = node[i].pos.z*mesoDimensions::length;
            toBePrinted_nodes[10*i +3] = node[i].vel.x*mesoDimensions::velocity;
            toBePrinted_nodes[10*i +4] = node[i].vel.y*mesoDimensions::velocity;
            toBePrinted_nodes[10*i +5] = node[i].vel.z*mesoDimensions::velocity;
            toBePrinted_nodes[10*i +6] = node[i].phi;
            toBePrinted_nodes[10*i +7] = force[i].x*mesoDimensions::force;
            toBePrinted_nodes[10*i +8] = force[i].y*mesoDimensions::force;
            toBePrinted_nodes[10*i +9] = force[i].z*mesoDimensions::force;
        }
    } else {
        if (params->calc_es == 0) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            #pragma omp parallel for default(none)
#endif
            for (int i = 0; i < num_nodes; i++) {
                toBePrinted_nodes[3*i   ] = node[i].pos.x*mesoDimensions::length;
                toBePrinted_nodes[3*i +1] = node[i].pos.y*mesoDimensions::length;
                toBePrinted_nodes[3*i +2] = node[i].pos.z*mesoDimensions::length;
            }
        } else {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            #pragma omp parallel for default(none)
#endif
            for (int i = 0; i < num_nodes; i++) {
                toBePrinted_nodes[4*i   ] = node[i].pos.x*mesoDimensions::length;
                toBePrinted_nodes[4*i +1] = node[i].pos.y*mesoDimensions::length;
                toBePrinted_nodes[4*i +2] = node[i].pos.z*mesoDimensions::length;
                toBePrinted_nodes[4*i +3] = node[i].phi;
            }
        }
    }
}

int Blob::read_nodes_from_file(FILE *trajectory_out) {
    char state_str[20];
    char *result = NULL;

    // If blob is static, don't read any nodes. Simply read the word "STATIC"
    if (blob_state == FFEA_BLOB_IS_STATIC) {
        result = fgets(state_str, 20, trajectory_out);
        if (result == NULL) {
            FFEA_ERROR_MESSG("Problem when reading 'STATIC' (expected) line in trajectory file\n")
        }
        if (strcmp(state_str, "STATIC\n") != 0) {
            FFEA_ERROR_MESSG("When restarting from trajectory file, expected to read 'STATIC', but instead found '%s...'\n", state_str)
        }
        return FFEA_OK;
    } else {
        result = fgets(state_str, 20, trajectory_out);
        if (result == NULL) {
            FFEA_ERROR_MESSG("Problem when reading state line in trajectory file\n")
        }
    }

    for (int i = 0; i < num_nodes; i++) {
        if (fscanf(trajectory_out, "%le %le %le %le %le %le %le %le %le %le\n", &node[i].pos.x, &node[i].pos.y, &node[i].pos.z, &node[i].vel.x, &node[i].vel.y, &node[i].vel.z, &node[i].phi, &force[i].x, &force[i].y, &force[i].z) != 10) {
            FFEA_ERROR_MESSG("(When restarting) Error reading from trajectory file, for node %d\n", i)
        } else {
            node[i].pos.x /= mesoDimensions::length;
            node[i].pos.y /= mesoDimensions::length;
            node[i].pos.z /= mesoDimensions::length;
            node[i].vel.x /= mesoDimensions::velocity;
            node[i].vel.y /= mesoDimensions::velocity;
            node[i].vel.z /= mesoDimensions::velocity;
            force[i].x /= mesoDimensions::force;
            force[i].y /= mesoDimensions::force;
            force[i].z /= mesoDimensions::force;
        }

    }
    return FFEA_OK;
}

int Blob::calculate_deformation() {

    int num_inversions = 0;
    matrix3 J;
    for (int n = 0; n < num_elements; n++) {

        // calculate jacobian for this element
        elem[n].calculate_jacobian(J);

        // get the 12 derivatives of the shape functions (by inverting the jacobian)
        // and also get the element volume. The function returns an error in the
        // case of an element inverting itself (determinant changing sign since last step)
        if (elem[n].calc_shape_function_derivatives_and_volume(J) == FFEA_ERROR) {
            FFEA_error_text();
            printf("Element %d has inverted during deformation calculation\n", n);
            num_inversions++;
        }

        // And F_ij
        elem[n].calc_deformation(J);
    }

    if(num_inversions > 0) {
        return FFEA_ERROR;
    } else {
        return FFEA_OK;
    }

}

scalar Blob::calc_volume() {

    scalar volume = 0.0;
    for(int i = 0; i < num_elements; ++i) {
        volume += elem[i].calc_volume();
    }
    return volume;
}

scalar Blob::calculate_strain_energy() {

    int n;
    scalar C, detF, strain_energy = 0.0;
    calculate_deformation();
    for(n = 0; n < num_elements; ++n) {
        C = elem[n].E - elem[n].G * 2.0 / 3.0;
        detF = elem[n].vol / elem[n].vol_0;
        strain_energy += elem[n].vol_0 * (elem[n].G * (mat3_double_contraction(elem[n].F_ij) - 3)
                                          + 0.5 * C * (detF * detF - 1)
                                          - ((2 * elem[n].G) + C) * log(detF)
                                         );
    }
    return 0.5 * strain_energy;
}

void Blob::make_measurements() {

    int n, i, j;
    scalar kenergy, senergy;
    vector3 total_ssint_xz_force;
    scalar total_ssint_xz_energy = 0.0, total_ssint_xz_area = 0.0;
    scalar temp1, temp2, temp3;
    scalar cogx = 0.0, cogy = 0.0, cogz = 0.0, lx = 0.0, ly = 0.0, lz = 0.0, brmsd = 0.0;
    scalar r[4][3];
    vector12 vec;

    kenergy = 0.0;
    senergy = 0.0;

    // OpenMP can't reduce members of classes :(
    //vector3_set_zero(L);
    //vector3_set_zero(CoG);
    vector3_set_zero(CoM);

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) reduction(+:kenergy, senergy) private(n, vec, temp1, temp2, temp3)
#endif
    for (n = 0; n < num_elements; n++) {

        if (linear_solver != FFEA_NOMASS_CG_SOLVER) {

            /*
             * Kinetic energy contribution:
             */

            // Read the u vector for this element
            elem[n].get_element_velocity_vector(vec);

            // Apply the mass matrix
            elem[n].apply_element_mass_matrix(vec);

            // Dot u with M.u to get the contribution to the kinetic energy
            kenergy += elem[n].n[0]->vel.x * vec[0] +
                       elem[n].n[1]->vel.x * vec[1] +
                       elem[n].n[2]->vel.x * vec[2] +
                       elem[n].n[3]->vel.x * vec[3] +
                       elem[n].n[0]->vel.y * vec[4] +
                       elem[n].n[1]->vel.y * vec[5] +
                       elem[n].n[2]->vel.y * vec[6] +
                       elem[n].n[3]->vel.y * vec[7] +
                       elem[n].n[0]->vel.z * vec[8] +
                       elem[n].n[1]->vel.z * vec[9] +
                       elem[n].n[2]->vel.z * vec[10] +
                       elem[n].n[3]->vel.z * vec[11];

            /*
             * Centre of Mass contribution:
             */
            //com_x += elem[n].mass * (elem[n].n[0]->pos.x + elem[n].n[1]->pos.x + elem[n].n[2]->pos.x + elem[n].n[3]->pos.x) / 4;
            //com_y += elem[n].mass * (elem[n].n[0]->pos.y + elem[n].n[1]->pos.y + elem[n].n[2]->pos.y + elem[n].n[3]->pos.y) / 4;
            //com_z += elem[n].mass * (elem[n].n[0]->pos.z + elem[n].n[1]->pos.z + elem[n].n[2]->pos.z + elem[n].n[3]->pos.z) / 4;
        }

        /*
         * Strain energy contribution:
         */

        scalar C = elem[n].E - (2.0 / 3.0) * elem[n].G;
        temp1 = elem[n].vol / elem[n].vol_0;
        senergy += elem[n].vol_0 * (elem[n].G * (mat3_double_contraction(elem[n].F_ij) - 3) + 0.5 * C * (temp1*temp1 - 1) - (C + 2 * elem[n].G) * log(temp1));
    }

    // And don't forget to multiply by a half
    kineticenergy = kenergy * 0.5;
    strainenergy = senergy * 0.5;
    get_CoM(&CoM);
    if (linear_solver != FFEA_NOMASS_CG_SOLVER) {

        // Get the centre of mass from the calculated totals
        //if (mass <= 0.0) {
        //   com_x = 0.0;
        //  com_y = 0.0;
        //com_z = 0.0;
        //} else {
        //    com_x *= 1.0 / mass;
        //    com_y *= 1.0 / mass;
        //    com_z *= 1.0 / mass;
        //}

        /* Calculate angular momentum */
        // mass matrix
        const matrix4 MM = {
            {.1, .05, .05, .05},
            {.05, .1, .05, .05},
            {.05, .05, .1, .05},
            {.05, .05, .05, .1}
        };

#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none) reduction(+:lx, ly, lz) shared(MM) private(n, r, temp1, temp2, temp3, i, j)
#endif
        for (n = 0; n < num_elements; n++) {

            // Find the separation vectors for this element
            for (i = 0; i < 4; i++) {
                r[i][0] = elem[n].n[i]->pos.x - CoM.x;
                r[i][1] = elem[n].n[i]->pos.y - CoM.y;
                r[i][2] = elem[n].n[i]->pos.z - CoM.z;
            }

            // Calculate contribution to angular momentum from this element
            temp1 = 0;
            temp2 = 0;
            temp3 = 0;
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    temp1 += MM[i][j] * (r[j][1] * elem[n].n[i]->vel.z - r[j][2] * elem[n].n[i]->vel.y);
                    temp2 += MM[i][j] * (r[j][2] * elem[n].n[i]->vel.x - r[j][0] * elem[n].n[i]->vel.z);
                    temp3 += MM[i][j] * (r[j][0] * elem[n].n[i]->vel.y - r[j][1] * elem[n].n[i]->vel.x);
                }
            }

            // Add this contribution to the sum
            lx += elem[n].mass * temp1;
            ly += elem[n].mass * temp2;
            lz += elem[n].mass * temp3;
        }
    }

    /* Calculate RMSD value for this configuration */
//#ifdef FFEA_PARALLEL_WITHIN_BLOB
//#pragma omp parallel for default(none) reduction(+:brmsd, cogx, cogy, cogz) private(i, temp1, temp2, temp3)
//#endif
    for (i = 0; i < num_nodes; i++) {

        /*
         * Center of geometry contribution (geometry better on nodes than elements)
         */
        cogx += node[i].pos.x;
        cogy += node[i].pos.y;
        cogz += node[i].pos.z;

    }
    CoG.x = cogx / num_nodes;
    CoG.y = cogy / num_nodes;
    CoG.z = cogz / num_nodes;

    // Remove translation and rotation (maybe make optional in future)
    bool remTrans = false;
    bool remRot = true;

    if(remTrans) {
        for(i = 0; i < num_nodes; ++i) {
            temp1 = node[i].pos.x - node[i].pos_0.x + CoG_0.x - CoG.x;
            temp2 = node[i].pos.y - node[i].pos_0.y + CoG_0.y - CoG.y;
            temp3 = node[i].pos.z - node[i].pos_0.z + CoG_0.z - CoG.z;
            brmsd += temp1 * temp1 + temp2 * temp2 + temp3*temp3;
        }
    } else {
        for(i = 0; i < num_nodes; ++i) {
            temp1 = node[i].pos.x - node[i].pos_0.x;
            temp2 = node[i].pos.y - node[i].pos_0.y;
            temp3 = node[i].pos.z - node[i].pos_0.z;
            brmsd += temp1 * temp1 + temp2 * temp2 + temp3*temp3;
        }
    }
    rmsd = sqrt(brmsd / num_nodes);
    L.x = lx;
    L.y = ly;
    L.z = lz;
}

void Blob::write_measurements_to_file(FILE *fout) {

    // White space for blob index bit
    fprintf(fout, "     ");
    if(there_is_mass()) {
        fprintf(fout, "%-14.6e", kineticenergy * mesoDimensions::Energy);
    }
    fprintf(fout, "%-14.6e", strainenergy * mesoDimensions::Energy);
    fprintf(fout, "%-14.6e%-14.6e%-14.6e%-14.6e", CoG.x * mesoDimensions::length, CoG.y * mesoDimensions::length, CoG.z * mesoDimensions::length, rmsd * mesoDimensions::length);
    fflush(fout);
}

void Blob::make_stress_measurements(FILE *stress_out, int blob_number) {
    int n;
    if (stress_out != NULL) {
        fprintf(stress_out, "blob\t%d\n", blob_number);
        for (n = 0; n < num_elements; n++) {
            fprintf(stress_out, "%e\n", elem[n].internal_stress_mag);
        }
        fprintf(stress_out, "\n");
    }
}

void Blob::calc_centroids_and_normals_of_all_faces() {
    int i;
    for (i = 0; i < num_surface_faces; i++)
        surface[i].calc_area_normal_centroid();
}

void Blob::calc_all_centroids() {

    int i;
    for (i = 0; i < num_surface_faces; i++) {
        surface[i].calc_area_normal_centroid();
    }
    for(i = 0; i < num_elements; ++i) {
        elem[i].calc_centroid();
    }

}

/*
 *
 */
int Blob::get_num_faces() {
    return num_surface_faces;
}

/*
 *
 */
int Blob::get_num_beads() {
    return num_beads;
}

bool Blob::is_using_beads() {
    if (num_beads > 0) return true;
    else return false;
}

/*
 *
 */
scalar Blob::get_rmsd() {
    return rmsd;
}

/*
 *
 */
int Blob::get_linear_solver() {
    return linear_solver;
};

/*
 *
arr3 Blob::get_CoG() {
    return CoG.data;
}
 */
void Blob::get_stored_centroid(arr3 &cog){
    arr3Store<scalar,arr3>(CoG.data, cog);
}

/*
 *
 */
Face * Blob::get_face(int i) {
    if (surface[i].is_ssint_active() == true) {
        return &surface[i];
    } else {
        return NULL;
    }
}

Face * Blob::absolutely_get_face(int i) {
    return &surface[i];
}

/*
 *
 */
tetra_element_linear *Blob::get_element(int i) {
    return &elem[i];
}

/**
 * @brief returns the position of bead i.
 *
 * @ingroup FMM
 **/
void Blob::get_bead_position(int i, arr3 &v) {

    arr3Store<scalar,arr3>( arr3_view<scalar,arr3>(bead_position+3*i,3), v);
}

/**
 * @brief returns the bead_type pointer.
 *
 * @ingroup FMM
 **/
int *Blob::get_bead_type_ptr() {
    return bead_type;

}

/**
 * @brief returns the bead_type of bead i.
 *
 * @ingroup FMM
 **/
int Blob::get_bead_type(int i) {
    return bead_type[i];

}

/**
 * @brief returns the list of nodes where bead i should be assigned to.
 *
 * @ingroup FMM
 **/
vector<int> Blob::get_bead_assignment(int i) {
    return bead_assignment[i];
}

scalar Blob::get_ssint_area() {
    scalar total_ssint_area = 0.0;
    for (int i = 0; i < get_num_faces(); ++i) {
        if (surface[i].is_ssint_active() == true) {
            total_ssint_area += surface[i].area;
        }
    }
    return total_ssint_area;
}
//		/*
//		 *
//		 */

int Blob::build_linear_node_elasticity_matrix(Eigen::SparseMatrix<scalar> *A) {

    int elem_index, a, b, i, j, global_a, global_a_lin, global_b, global_b_lin;
    int row, column;
    scalar dx, val; // dxtemp;
    // matrix3 J, stress;
    vector12 elastic_force[2];
    vector<Eigen::Triplet<scalar> > components;

    // Firstly, get a mapping from all node indices to just linear node indices
    // int num_linear_nodes = get_num_linear_nodes();
    int map[num_nodes];
    j = 0;
    for(i = 0; i < num_nodes; ++i) {
        if(node[i].am_I_linear()) {
            map[i] = j;
            j++;
        } else {
            map[i] = -1;
        }
    }

    // For each element
    for(elem_index = 0; elem_index < num_elements; ++elem_index) {

        // Calculate dx, how far each node should be moved for a linearisataion, as well as unstrained parameter
        dx = cbrt(elem[elem_index].calc_volume()) / 1000.0;

        // For every node a in every direction i
        for(a = 0; a < 4; ++a) {

            // Get global index for node a
            global_a = elem[elem_index].n[a]->index;
            global_a_lin = map[elem[elem_index].n[a]->index];

            for(i = 0; i < 3; ++i) {

                // Move node a in direction i and calculate the change

                // Move
                node[global_a].move(i, dx);

                // Calculate
                elem[elem_index].calc_elastic_force_vector(elastic_force[1]);

                // Other side
                node[global_a].move(i, -2 * dx);

                // Recalculate
                elem[elem_index].calc_elastic_force_vector(elastic_force[0]);

                // Move back to start
                node[global_a].move(i, dx);

                // Now, how has each component changed because of this change?
                // For the component representing node b in direction j
                for(b = 0; b < 4; ++b) {

                    // Get global index for node b
                    global_b = elem[elem_index].n[b]->index;
                    global_b_lin = map[elem[elem_index].n[b]->index];

                    for(j = 0; j < 3; ++j) {
                        val = (1.0 / (2 * dx)) * (elastic_force[1][4 * j + b] - elastic_force[0][4 * j + b]);

                        // Row is dE_p, column dx_q. Not that it should matter! Directions then nodes i.e. x0,y0,z0,x1,y1,z1...xn,yn,zn
                        //row = num_linear_nodes * i + global_a;
                        //column = num_linear_nodes * j + global_b;
                        row = 3 * global_a_lin + i;
                        column = 3 * global_b_lin + j;
                        components.push_back(Eigen::Triplet<scalar>(row, column, val));
                    }
                }
            }
        }
    }

    // Now build the matrix
    A->setFromTriplets(components.begin(), components.end());
    return FFEA_OK;
}

int Blob::build_linear_node_viscosity_matrix(Eigen::SparseMatrix<scalar> *K) {

    int i, j, a, b, global_a, global_b, row, column;
    int elem_index, num_linear_nodes;
    scalar val;
    matrix3 J;
    vector<Eigen::Triplet<scalar> > components;

    // Firstly, get a mapping from all node indices to just linear node indices
    int map[num_nodes];
    int offset = 0;
    for(int i = 0; i < num_nodes; ++i) {
        if(node[i].am_I_linear()) {
            map[i] = i - offset;
        } else {
            offset += 1;
            map[i] = -1;
        }
    }

    num_linear_nodes = get_num_linear_nodes();

    // For each element
    for(elem_index = 0; elem_index < num_elements; ++elem_index) {

        // Calculate a local viscosity matrix
        elem[elem_index].calculate_jacobian(J);
        elem[elem_index].calc_shape_function_derivatives_and_volume(J);
        elem[elem_index].create_viscosity_matrix();

        // Add each component of the local matrix to the global matrix

        //For each node
        for(a = 0; a < 4; ++a) {
            global_a = map[elem[elem_index].n[a]->index];

            for(b = 0; b < 4; ++b) {
                global_b = map[elem[elem_index].n[b]->index];

                // And each direction
                for(i = 0; i < 3; ++i) {
                    for(j = 0; j < 3; ++j) {
                        val = elem[elem_index].viscosity_matrix[4 * i + a][4 * j + b];
                        //row = num_linear_nodes * i + global_a;
                        //column = num_linear_nodes * j + global_b;
                        row = 3 * global_a + i;
                        column = 3 * global_b + j;
                        components.push_back(Eigen::Triplet<scalar>(row, column, val));
                    }
                }
            }
        }
    }

    // For each node, add stokes if necessary
    if(params->calc_stokes == 1) {
        for(global_a = 0; global_a < num_nodes; ++global_a) {
            if(map[global_a] == -1) {
                continue;
            } else {
                for(i = 0; i < 3; ++i) {
                    val = node[map[global_a]].stokes_drag;
                    //row = 4 * i + map[global_a];
                    row = 3 * map[global_a] + i;
                    column = row;
                    components.push_back(Eigen::Triplet<scalar>(row, column, val));
                }
            }
        }
    }

    // Now build the matrix and get the symmetric part (just in case)
    K->setFromTriplets(components.begin(), components.end());
    components.clear();
    for(j = 0; j < K->outerSize(); ++j) {
        for(Eigen::SparseMatrix<scalar>::InnerIterator it(*K,j); it; ++it) {
            components.push_back(Eigen::Triplet<scalar>(it.row(), it.col(), 0.5 * it.value()));
            components.push_back(Eigen::Triplet<scalar>(it.col(), it.row(), 0.5 * it.value()));
        }
    }
    K->setFromTriplets(components.begin(), components.end());
    return FFEA_OK;
}

int Blob::build_linear_node_mass_matrix(Eigen::SparseMatrix<scalar> *M) {

    int i, j, a, b, global_a, global_b, row, column;
    int elem_index, num_linear_nodes;
    scalar val;
    //MassMatrixQuadratic M_el;
    MassMatrixLinear M_el;
    vector<Eigen::Triplet<scalar> > components;

    // Firstly, get a mapping from all node indices to just linear node indices

    int map[num_nodes];
    int offset = 0;
    for(i = 0; i < num_nodes; ++i) {

        if(node[i].am_I_linear()) {
            //if(true) {
            map[i] = i - offset;
        } else {
            offset += 1;
            map[i] = -1;
        }
    }

    num_linear_nodes = get_num_linear_nodes();

    // For each element
    scalar sum = 0.0, sum2 = 0.0;
    for(elem_index = 0; elem_index < num_elements; ++elem_index) {

        // Build mass matrix
        elem[elem_index].construct_element_mass_matrix(&M_el);

        // Add each component of the local matrix to the global matrix

        // For each node
        for(a = 0; a < 4; ++a) {
            global_a = map[elem[elem_index].n[a]->index];

            for(b = 0; b < 4; ++b) {
                global_b = map[elem[elem_index].n[b]->index];

                // From another function that get's actual memory location. Dereference for value
                // Val same for each direction
                val = M_el.get_M_alpha_value(a, b);
                //sum2 += val;
                //cout << a << " " << b << " " << val * mesoDimensions::mass << endl;
                // And each direction
                for(i = 0; i < 3; ++i) {

                    row = 3 * global_a + i;
                    column = 3 * global_b + i;
                    components.push_back(Eigen::Triplet<scalar>(row, column, val));
                }
            }
        }
    }

    // Now build the matrix and get the symmetric part (just in case)
    M->setFromTriplets(components.begin(), components.end());
    components.clear();
    for(j = 0; j < M->outerSize(); ++j) {
        for(Eigen::SparseMatrix<scalar>::InnerIterator it(*M,j); it; ++it) {
            components.push_back(Eigen::Triplet<scalar>(it.row(), it.col(), 0.5 * it.value()));
            components.push_back(Eigen::Triplet<scalar>(it.col(), it.row(), 0.5 * it.value()));
            //sum += it.value();
        }
    }
    //cout << "Blob mass matrix sums: " << mesoDimensions::mass * sum2 << " " << mesoDimensions::mass * sum << endl;
    M->setFromTriplets(components.begin(), components.end());
    return FFEA_OK;

}

int Blob::build_linear_node_rp_diffusion_matrix(Eigen_MatrixX *D) {

    int i, j, n, m, num_linear_nodes;
    scalar mod, mod2, a;
    Eigen_Matrix3 block, rr;
    Eigen_Vector3 sep;

    // Firstly, get a mapping from all node indices to just linear node indices
    int map[num_nodes];
    int offset = 0;
    for(i = 0; i < num_nodes; ++i) {
        if(node[i].am_I_linear()) {
            map[i] = i - offset;
        } else {
            offset += 1;
            map[i] = -1;
        }
    }

    num_linear_nodes = get_num_linear_nodes();

    // Now build the matrix, one element at a time...
    D->setZero();

    // For each pair of nodes (upper triangle only)
    for(n = 0; n < num_nodes; ++n) {

        // If secondary node, continue
        if(map[n] == -1) {
            continue;
        } else {
            i = map[n];
        }
        for(m = n; m < num_nodes; ++m) {

            // If secondary node, continue
            if(map[m] == -1) {
                continue;
            } else {
                j = map[m];
            }

            // Initialise directional block matrix
            block.setZero();

            // If we are on the block diagonal
            if(n == m) {
                block = Eigen::Matrix3d::Identity() *  params->kT / node[n].stokes_drag;
            } else {

                // Get distance between nodes and the outer product of the separation
                Eigen::Vector3d vecn(node[n].pos.x, node[n].pos.y, node[n].pos.z);
                Eigen::Vector3d vecm(node[m].pos.x, node[m].pos.y, node[m].pos.z);

                sep = vecn - vecm;
                mod = sep.norm();
                mod2 = mod * mod;
                rr = sep * sep.transpose();

                // Effective stokes radius. Will usually equal a anyway

                a = (node[n].stokes_radius + node[m].stokes_radius) / 2.0;

                // Condition for positive-definiteness
                if(mod > 2 * node[n].stokes_radius && mod > 2 * node[m].stokes_radius) {
                    block = Eigen::Matrix3d::Identity() * mod2/ 3.0;

                    block -= rr;
                    block *= 2 * a * a / mod2;
                    block += Eigen_Matrix3::Identity() * mod2;
                    block += rr;
                    block *= params->kT / (8 * ffea_const::pi * params->stokes_visc * mod2 * mod);

                } else {
                    block = Eigen::Matrix3d::Identity() * (1 - ((9 * mod) / (32.0 * a)));
                    block += rr * (3 / (32.0 * a * mod));
                    block *= params->kT / (6 * ffea_const::pi * params->stokes_visc * a);
                }
            }

            // Linear node position from global nodes
            D->block<3,3>(3 * i, 3 * j) = block;
        }
    }

    // Fill in the lower trangle
    for(i = 0; i < 3 * num_linear_nodes; ++i) {
        for(j = i; j < 3 * num_linear_nodes; ++j) {
            (*D)(j, i) = (*D)(i, j);
        }
    }
    return FFEA_OK;
}

int Blob::solve_poisson(scalar *phi_gamma_IN, scalar *J_Gamma_OUT) {
    if (num_interior_nodes > 0) {
        /* Convert the given potential on the surface faces into the potential on each node */
        for (int i = 0; i < num_surface_nodes; i++) {
            phi_Gamma[i] = 0;
        }

        for (int i = 0; i < num_surface_faces; i++) {
            for (int j = 0; j < 3; j++) {
                phi_Gamma[surface[i].n[j]->index] += phi_gamma_IN[i];
            }
        }

        for (int i = 0; i < num_surface_nodes; i++) {
            phi_Gamma[i] /= num_contributing_faces[i];
            //			printf("phi_Gamma[%d] = %e\n", i, phi_Gamma[i]);
        }

        /* Calculate the RHS vector */
        poisson_surface_matrix->apply(phi_Gamma, poisson_rhs);
        for (int i = 0; i < num_interior_nodes; i++) {
            poisson_rhs[i] = q[i + num_surface_nodes] - poisson_rhs[i];
            //			printf("poisson_rhs[%d] = %e\n", i, poisson_rhs[i]);
        }

        poisson_solver->solve(poisson_interior_matrix, phi_Omega, poisson_rhs);

        //		printf("int:\n");
        //		poisson_interior_matrix->print_dense();

        for (int n = 0; n < num_surface_nodes; n++)
            node[n].phi = phi_Gamma[n];

        for (int n = 0; n < num_interior_nodes; n++) {
            node[n + num_surface_nodes].phi = phi_Omega[n];
            //			printf("phi_Omega[%d] = %e\n", n, phi_Omega[n]);
        }

        for (int i = 0; i < num_surface_faces; i++) {
            J_Gamma_OUT[i] = surface[i].get_normal_flux();
        }
    }

    return FFEA_OK;
}


int Blob::apply_ctforces() {

    // CTF 1 - Add the linear constant forces, stored in ctf_force.
    // If there are no forces, num_l_ctf will be zero.
    for (int i=0; i<num_l_ctf; i++) {
        force[ctf_l_nodes[i]].x += ctf_l_forces[3*i  ];
        force[ctf_l_nodes[i]].y += ctf_l_forces[3*i+1];
        force[ctf_l_nodes[i]].z += ctf_l_forces[3*i+2];
    }
    // CTF 2 - Add the rotational forces:
    for (int i=0; i<num_r_ctf; i++) {
        // 2.1 - get the axis:
        arr3 r, f;
        scalar rsize;
        // 2.1.1 - if it is defined by points:
        if (ctf_r_type[2*i] == 'p') {
            // get r, and its length rsize
            intersectingPointToLine<scalar,arr3>(node[ctf_r_nodes[i]].pos,// point out of the line
                                                 arr3_view<scalar,arr3>(ctf_r_axis+6*i, 3), //  point
                                                 arr3_view<scalar,arr3>(ctf_r_axis+(6*i+3), 3), r); // vector, axis
            arr3arr3Substract<scalar,arr3>(r, node[ctf_r_nodes[i]].pos.data, r);
        } /* else if (ctf_l_type[2*i] == 'n') {
        // not working yet!
        that should read:
        axis[0] = blob[ctf_r_axis[6*i + 3]].conf[ctf_r_axis[6*i + 4]].node[ctf_r_axis[6*i + 5]]->pos.x -
                  blob[ctf_r_axis[6*i + 0]].conf[ctf_r_axis[6*i + 1]].node[ctf_r_axis[6*i + 2]]->pos.x
        axis[1] =  ...
        axis[2] =  ...
        arr3Normalize(axis);
        // We do not have this information in "Blob"!!
      } */
        arr3arr3VectorProduct<scalar,arr3>(r, arr3_view<scalar,arr3>(ctf_r_axis+(6*i+3), 3), f);
        // 2.2.a - apply a constant circular force:
        if ( ctf_r_type[2*i + 1] == 'f' ) {
            arr3Resize<scalar,arr3>(ctf_r_forces[i], f);
            // 2.2.b - apply a constant torque:
        } else if (ctf_r_type[2*i + 1] == 't') {
            rsize = mag<scalar,arr3>(r);
            arr3Resize<scalar,arr3>(ctf_r_forces[i]/rsize, f);
        }
        force[ctf_r_nodes[i]].x += f[0];
        force[ctf_r_nodes[i]].y += f[1];
        force[ctf_r_nodes[i]].z += f[2];
    }

    // CTF 3 - Add the linear surface forces:
    scalar totalArea;
    arr3 traction;
    int auxndx = 0;
    for (int i=0; i<num_slsets_ctf; i++) {
        totalArea = 0;
        scalar *faceAreas = new scalar[ctf_sl_surfsize[i]];
        for (int j=0; j<ctf_sl_surfsize[i]; j++) {
            faceAreas[j] = surface[ctf_sl_faces[auxndx + j]].get_area();
            totalArea += faceAreas[j];
        }
        for (int j=0; j<ctf_sl_surfsize[i]; j++) {
            arr3Resize2<scalar,arr3>(faceAreas[j]/(3*totalArea), arr3_view<scalar,arr3>(ctf_sl_forces+3*i, 3), traction);
            surface[ctf_sl_faces[auxndx+j]].add_force_to_node(0, traction);
            surface[ctf_sl_faces[auxndx+j]].add_force_to_node(1, traction);
            surface[ctf_sl_faces[auxndx+j]].add_force_to_node(2, traction);
        }
        auxndx += ctf_sl_surfsize[i];
        delete[] faceAreas;
    }

	return FFEA_OK;

}

/*
 */
void Blob::zero_force() {
    for (int i = 0; i < num_elements; i++) {
        elem[i].zero_force();
    }
    for (int i = 0; i < num_surface_faces; i++) {
        surface[i].zero_force();
    }
}

void Blob::set_forces_to_zero() {
    for (int i = 0; i < num_nodes; ++i) {
        force[i][0] = 0;
        force[i][1] = 0;
        force[i][2] = 0;
    }
}

void Blob::get_node(int index, arr3 &v) {

    arr3Store<scalar,arr3>(node[index].pos.data, v);

}

void Blob::get_node_0(int index, arr3 &v) {

    arr3Store<scalar,arr3>(node[index].pos_0.data, v);

}

void Blob::add_force_to_node(vector3 f, int index) {
    force[index].x += f.x;
    force[index].y += f.y;
    force[index].z += f.z;
}

void Blob::velocity_all(scalar vel_x, scalar vel_y, scalar vel_z) {
    int i;
    for (i = 0; i < num_nodes; i++) {
        node[i].vel.x = vel_x;
        node[i].vel.y = vel_y;
        node[i].vel.z = vel_z;
    }
}

void Blob::build_poisson_matrices() {
    if (num_interior_nodes > 0) {
        /*
                        scalar K[num_nodes*num_nodes];
                        for(int i = 0; i < num_nodes * num_nodes; i++) {
                                K[i] = 0;
                        }
         */

        /* Calculate the K_alpha matrices for each element (the diffusion matrix * epsilon * volume) */
        for (int n = 0; n < num_elements; n++) {
            elem[n].calculate_K_alpha();
            //			elem[n].add_K_alpha(K, num_nodes);
        }


        /*
                        printf("BLOB Full K:\n");
                        int c = 0;
                        for(int i = 0; i < num_nodes; i++) {
                                for(int j = 0; j < num_nodes; j++) {
                                        printf("%+e ", K[c]);
                                        c++;
                                }
                                printf("\n");
                        }
         */
        /* Construct the poisson matrices for the current blob (based on diffusion matrices of elements) */
        poisson_surface_matrix->build();
        poisson_interior_matrix->build();
    }
}

//		int elements_are_connected(int e1, int e2);

/*

 */
scalar Blob::get_mass() {
    return mass;
}

/*
 */
void Blob::enforce_box_boundaries(vector3 *box_dim) {
    if (params->wall_x_1 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos.x < 0) {
                if (node[i].vel.x < 0) {
                    node[i].vel.x = 0;
                }
            }
        }
    }
    if (params->wall_x_2 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos.x > box_dim->x) {
                if (node[i].vel.x > 0) {
                    node[i].vel.x = 0;
                }
            }
        }
    }
    if (params->wall_y_1 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos.y < 0) {
                if (node[i].vel.y < 0) {
                    node[i].vel.y = 0;
                }
            }
        }
    }
    if (params->wall_y_2 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos.y > box_dim->y) {
                if (node[i].vel.y > 0) {
                    node[i].vel.y = 0;
                }
            }
        }
    }
    if (params->wall_z_1 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos.z < 0) {
                if (node[i].vel.z < 0) {
                    node[i].vel.z = 0;
                }
            }
        }
    }
    if (params->wall_z_2 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos.z > box_dim->z) {
                if (node[i].vel.z > 0) {
                    node[i].vel.z = 0;
                }
            }
        }
    }
}

int Blob::get_num_nodes() {
    return num_nodes;
}

int Blob::get_num_elements() {
    return num_elements;
}

int Blob::get_motion_state() {
    return blob_state;
}

scalar Blob::get_scale() {
    return scale;
}

scalar Blob::get_RandU01() {
    return rng->RandU01();
}

int Blob::get_num_linear_nodes() {

    int n, i;
    set<int> node_indices;
    for(n = 0; n < num_elements; ++n) {
        for(i = 0; i < 4; ++i) {
            node_indices.insert(elem[n].n[i]->index);
        }
    }

    return node_indices.size();
}

scalar Blob::get_kinetic_energy() {
    return kineticenergy;
}

scalar Blob::get_strain_energy() {
    //std::cout << "Blob " << blob_index << " strain energy: " << strainenergy << "\n";
    return strainenergy;
}

void Blob::get_min_max(vector3 *blob_min, vector3 *blob_max) {

    blob_min->x = INFINITY;
    blob_max->x = -1 * INFINITY;
    blob_min->y = INFINITY;
    blob_max->y = -1 * INFINITY;
    blob_min->z = INFINITY;
    blob_max->z = -1 * INFINITY;

    for(int i = 0; i < num_nodes; ++i) {
        if(node[i].pos.x > blob_max->x) {
            blob_max->x = node[i].pos.x;
        } else if (node[i].pos.x < blob_min->x) {
            blob_min->x = node[i].pos.x;
        }

        if(node[i].pos.y > blob_max->y) {
            blob_max->y = node[i].pos.y;
        } else if (node[i].pos.y < blob_min->y) {
            blob_min->y = node[i].pos.y;
        }

        if(node[i].pos.z > blob_max->z) {
            blob_max->z = node[i].pos.z;
        } else if (node[i].pos.z < blob_min->z) {
            blob_min->z = node[i].pos.z;
        }
    }
}

/*
 */
int Blob::load_nodes(const char *node_filename, scalar scale) {
    FILE *in = NULL;
    int i;
    double x, y, z;
    const int max_line_size = 50;
    char line[max_line_size];

    // open the node file
    if ((in = fopen(node_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(node_filename)
    }
    printf("\t\tReading in nodes file: %s\n", node_filename);

    // first line should be the file type "ffea node file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of node file\n")
    }
    if (strcmp(line, "walrus node file\n") != 0 && strcmp(line, "ffea node file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea node file' (read '%s') \n", line)
    }

    // read in the number of nodes in the file
    if (fscanf(in, "num_nodes %d\n", &num_nodes) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of nodes\n")
    }
    printf("\t\t\tNumber of nodes = %d\n", num_nodes);

    // read in the number of surface nodes in the file
    if (fscanf(in, "num_surface_nodes %d\n", &num_surface_nodes) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of surface nodes\n")
    }
    printf("\t\t\tNumber of surface nodes = %d\n", num_surface_nodes);

    // read in the number of interior nodes in the file
    if (fscanf(in, "num_interior_nodes %d\n", &num_interior_nodes) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of interior nodes\n")
    }
    printf("\t\t\tNumber of interior nodes = %d\n", num_interior_nodes);

    // Allocate the memory for all these nodes
    node = new(std::nothrow) mesh_node[num_nodes];
    node_position = new(std::nothrow) vector3*[num_nodes];

    if (node == NULL || node_position == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Unable to allocate memory for nodes array.\n")
    }

    // Check for "surface nodes:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error when looking for 'surface nodes:' line\n")
    }
    if (strcmp(line, "surface nodes:\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("Could not find 'surface nodes:' line (found '%s' instead)\n", line)
    }

    // Read in all the surface nodes from file
    for (i = 0; i < num_nodes; i++) {

        if (i == num_surface_nodes) {
            // Check for "interior nodes:" line
            if (fgets(line, max_line_size, in) == NULL) {
                fclose(in);
                FFEA_ERROR_MESSG("Error when looking for 'interior nodes:' line\n")
            }
            if (strcmp(line, "interior nodes:\n") != 0) {
                fclose(in);
                FFEA_ERROR_MESSG("Could not find 'interior nodes:' line (found '%s' instead)\n", line)
            }
        }

        if (fscanf(in, "%le %le %le\n", &x, &y, &z) != 3) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from nodes file at node %d\n", i)
        } else {
            node[i].pos.x = scale * (scalar) x;
            node[i].pos.y = scale * (scalar) y;
            node[i].pos.z = scale * (scalar) z;
            node_position[i] = &node[i].pos;
            node[i].vel.x = 0;
            node[i].vel.y = 0;
            node[i].vel.z = 0;

            node[i].index = i;
        }
    }

    fclose(in);
    printf("\t\t\tRead %d nodes from %s\n", i, node_filename);

    return FFEA_OK;
}

/*
 */
int Blob::load_topology(const char *topology_filename) {
    FILE *in = NULL;
    int i, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
    const int max_line_size = 50;
    char line[max_line_size];

    // Now open the topology file
    if ((in = fopen(topology_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(topology_filename)
    }
    printf("\t\tReading in topology file: %s\n", topology_filename);

    // first line should be the file type "ffea topology file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of topology file\n")
    }
    if (strcmp(line, "walrus topology file\n") != 0 && strcmp(line, "ffea topology file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea topology file' (read '%s') \n", line)
    }

    // read in the total number of elements in the file
    if (fscanf(in, "num_elements %d\n", &num_elements) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of elements\n")
    }
    printf("\t\t\tNumber of elements = %d\n", num_elements);

    // read in the number of surface elements in the file
    if (fscanf(in, "num_surface_elements %d\n", &num_surface_elements) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of surface elements\n")
    }
    printf("\t\t\tNumber of surface elements = %d\n", num_surface_elements);

    // read in the number of interior elements in the file
    if (fscanf(in, "num_interior_elements %d\n", &num_interior_elements) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of interior elements\n")
    }
    printf("\t\t\tNumber of interior elements = %d\n", num_interior_elements);

    // Allocate the memory for all these elements
    elem = new(std::nothrow) tetra_element_linear[num_elements];
    if (elem == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Unable to allocate memory for element array.\n")
    }

    // Check for "surface elements:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error when looking for 'surface elements:' line\n")
    }
    if (strcmp(line, "surface elements:\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("Could not find 'surface elements:' line (found '%s' instead)\n", line)
    }

    // Read in all the elements from file
    for (i = 0; i < num_elements; i++) {
        if (i == num_surface_elements) {
            // Check for "interior elements:" line
            if (fgets(line, max_line_size, in) == NULL) {
                fclose(in);
                FFEA_ERROR_MESSG("Error when looking for 'interior elements:' line\n")
            }
            if (strcmp(line, "interior elements:\n") != 0) {
                fclose(in);
                FFEA_ERROR_MESSG("Could not find 'interior elements:' line (found '%s' instead)\n", line)
            }
        }

        if (fscanf(in, "%d %d %d %d %d %d %d %d %d %d\n", &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &n9, &n10) != 10) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from elements file at element %d\n", i)
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= num_nodes ||
                    n2 < 0 || n2 >= num_nodes ||
                    n3 < 0 || n3 >= num_nodes ||
                    n4 < 0 || n4 >= num_nodes ||
                    n5 < 0 || n5 >= num_nodes ||
                    n6 < 0 || n6 >= num_nodes ||
                    n7 < 0 || n7 >= num_nodes ||
                    n8 < 0 || n8 >= num_nodes ||
                    n9 < 0 || n9 >= num_nodes ||
                    n10 < 0 || n10 >= num_nodes) {
                fclose(in);
                FFEA_ERROR_MESSG("Error: Element %d references an out of bounds node index\n", i)
            }

            // Link element nodes to actual nodes
            elem[i].n[0] = &node[n1];
            elem[i].n[1] = &node[n2];
            elem[i].n[2] = &node[n3];
            elem[i].n[3] = &node[n4];
            elem[i].n[4] = &node[n5];
            elem[i].n[5] = &node[n6];
            elem[i].n[6] = &node[n7];
            elem[i].n[7] = &node[n8];
            elem[i].n[8] = &node[n9];
            elem[i].n[9] = &node[n10];

            // Assign bool to linear nodes
            node[n1].set_linear();
            node[n2].set_linear();
            node[n3].set_linear();
            node[n4].set_linear();
            elem[i].daddy_blob = this;
            elem[i].index = i;
        }
    }

    fclose(in);

    if (i == 1)
        printf("\t\t\tRead 1 element from %s\n", topology_filename);
    else
        printf("\t\t\tRead %d elements from %s\n", i, topology_filename);

    return FFEA_OK;
}

/*
 */
int Blob::load_surface(const char *surface_filename, SimulationParams* params) {
    FILE *in = NULL;
    int i, n1, n2, n3, element;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(surface_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(surface_filename)
    }
    printf("\t\tReading in surface file: %s\n", surface_filename);

    // first line should be the file type "ffea surface file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of surface file\n")
    }
    if (strcmp(line, "walrus surface file\n") != 0 && strcmp(line, "ffea surface file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea surface file' (read '%s') \n", line)
    }

    // read in the number of faces in the file
    if (fscanf(in, "num_surface_faces %d\n", &num_surface_faces) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of faces\n")
    }
    printf("\t\t\tNumber of faces = %d\n", num_surface_faces);

    // Allocate the memory for all these faces
    surface = new(std::nothrow) Face[num_surface_faces];
    if (surface == NULL) FFEA_ERROR_MESSG("Failed to allocate memory for the faces\n");

    // Check for "faces:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error when looking for 'faces:' line\n")
    }
    if (strcmp(line, "faces:\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("Could not find 'faces:' line (found '%s' instead)\n", line)
    }

    // Read in all the faces from file
    scalar smallest_A = INFINITY;
    for (i = 0; i < num_surface_faces; i++) {
        if (fscanf(in, "%d %d %d %d\n", &element, &n1, &n2, &n3) != 4) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from surface file at face %d. There should be 4 space separated integers. \n", i);
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= num_nodes ||
                    n2 < 0 || n2 >= num_nodes ||
                    n3 < 0 || n3 >= num_nodes) {
                fclose(in);
                FFEA_ERROR_MESSG("Error: Surface face %d references an out of bounds node index\n", i);
            } else if (element < 0 || element >= num_elements) {
                fclose(in);
                FFEA_ERROR_MESSG("Error: Surface face %d references an out of bounds element index\n", i);
            }


            int n1_el = elem[element].what_node_is_this(n1), n2_el = elem[element].what_node_is_this(n2), n3_el = elem[element].what_node_is_this(n3);

            SecondOrderFunctions::stu n1_stu = {SecondOrderFunctions::stu_lookup[n1_el].s, SecondOrderFunctions::stu_lookup[n1_el].t, SecondOrderFunctions::stu_lookup[n1_el].u};
            SecondOrderFunctions::stu n2_stu = {SecondOrderFunctions::stu_lookup[n2_el].s, SecondOrderFunctions::stu_lookup[n2_el].t, SecondOrderFunctions::stu_lookup[n2_el].u};
            SecondOrderFunctions::stu n3_stu = {SecondOrderFunctions::stu_lookup[n3_el].s, SecondOrderFunctions::stu_lookup[n3_el].t, SecondOrderFunctions::stu_lookup[n3_el].u};

            SecondOrderFunctions::stu centroid_stu = {
                (n1_stu.s + n2_stu.s + n3_stu.s) / 3.0,
                (n1_stu.t + n2_stu.t + n3_stu.t) / 3.0,
                (n1_stu.u + n2_stu.u + n3_stu.u) / 3.0
            };

            int n_op = elem[element].get_opposite_node(n1_el, n2_el, n3_el);
            if (n_op == -1)
                FFEA_ERROR_MESSG("Error: Could not find the opposite node\n");
            // now the node that we can pass is: elem[element].n[n_op]
            // elem[element].n[n1_el]->print()  =  node[n1].print();


            surface[i].init(i, &elem[element], &node[n1], &node[n2], &node[n3], elem[element].n[n_op], centroid_stu, this, params);
            if (surface[i].area_0 < smallest_A) {
                smallest_A = surface[i].area_0;
            }
        }
    }
    fclose(in);

    printf("\t\t\tSmallest Face Area = %e\n", smallest_A * mesoDimensions::area);
    if (i == 1)
        printf("\t\t\tRead 1 surface face from %s\n", surface_filename);
    else
        printf("\t\t\tRead %d surface faces from %s\n", i, surface_filename);

    return FFEA_OK;
}

/*
 */
int Blob::load_surface_no_topology(const char *surface_filename, SimulationParams *params) {
    FILE *in = NULL;
    int i, n1, n2, n3, element;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(surface_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(surface_filename)
    }
    printf("\t\tReading in surface file: %s\n", surface_filename);

    // first line should be the file type "ffea surface file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of surface file\n")
    }
    if (strcmp(line, "walrus surface file\n") != 0 && strcmp(line, "ffea surface file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea surface file' (read '%s') \n", line)
    }

    // read in the number of faces in the file
    if (fscanf(in, "num_surface_faces %d\n", &num_surface_faces) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of faces\n")
    }
    printf("\t\t\tNumber of faces = %d\n", num_surface_faces);

    // Allocate the memory for all these faces
    surface = new(std::nothrow) Face[num_surface_faces];
    if (surface == NULL) FFEA_ERROR_MESSG("Failed to allocate memory for the faces\n");

    // Check for "faces:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error when looking for 'faces:' line\n")
    }
    if (strcmp(line, "faces:\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("Could not find 'faces:' line (found '%s' instead)\n", line)
    }

    // Read in all the faces from file (element will always be zero here, because no internal structure exists)
    scalar smallest_A = INFINITY;
    for (i = 0; i < num_surface_faces; i++) {
        if (fscanf(in, "%d %d %d %d\n", &element, &n1, &n2, &n3) != 4) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from surface file at face %d. There should be 4 space separated integers. \n", i);
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= num_nodes ||
                    n2 < 0 || n2 >= num_nodes ||
                    n3 < 0 || n3 >= num_nodes) {
                fclose(in);
                FFEA_ERROR_MESSG("Error: Surface face %d references an out of bounds node index\n", i);
            }

            surface[i].init(i, &node[n1], &node[n2], &node[n3], NULL, this, params);

            if (surface[i].area_0 < smallest_A) {
                smallest_A = surface[i].area_0;
            }
        }
    }

    fclose(in);

    printf("\t\t\tSmallest Face Area = %e\n", smallest_A);
    if (i == 1)
        printf("\t\t\tRead 1 surface face from %s\n", surface_filename);
    else
        printf("\t\t\tRead %d surface faces from %s\n", i, surface_filename);

    return FFEA_OK;
}

/*
 */
int Blob::load_material_params(const char *material_params_filename) {
    FILE *in = NULL;
    int num_material_elements;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(material_params_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(material_params_filename)
    }
    printf("\t\tReading in material parameters file: %s\n", material_params_filename);

    // first line should be the file type "ffea material params file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of material params file\n")
    }
    if (strcmp(line, "walrus material params file\n") != 0 && strcmp(line, "ffea material params file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea material params file' (read '%s') \n", line)
    }

    // read in the number of elements in the file
    if (fscanf(in, "num_elements %d\n", &num_material_elements) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of elements\n")
    }
    printf("\t\t\tNumber of elements in material params file = %d\n", num_material_elements);

    // Check that we have same number of elements in material params file as in topology file
    if (num_material_elements != num_elements) {
        fclose(in);
        FFEA_ERROR_MESSG("Number of elements in material params file (%d) does not match number of elements in topology file (%d)\n", num_material_elements, num_elements)
    }

    // Set the material parameters for each element in the Blob
    scalar density = 0.0, shear_visc = 0.0, bulk_visc = 0.0, shear_mod = 0.0, bulk_mod = 0.0, dielectric = 0.0;
    int i;
    for (i = 0; i < num_elements; i++) {
        if (fscanf(in, "%le %le %le %le %le %le\n", &density, &shear_visc, &bulk_visc, &shear_mod, &bulk_mod, &dielectric) != 6) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from material params file at element %d. There should be 6 space separated real values (density, shear_visc, bulk_visc, shear_mod, bulk_mod, dielectric).\n", i);
        }
        elem[i].rho = density * mesoDimensions::volume / mesoDimensions::mass ;
        elem[i].A = shear_visc / (mesoDimensions::pressure * mesoDimensions::time);
        elem[i].B = bulk_visc / (mesoDimensions::pressure * mesoDimensions::time) - (2.0 / 3.0) * shear_visc; // Code uses second coefficient of viscosity
        elem[i].G = shear_mod / mesoDimensions::pressure;
        elem[i].E = bulk_mod / mesoDimensions::pressure;
        elem[i].dielectric = dielectric; // relative permittivity.
    }

    fclose(in);

    printf("\t\t\tRead %d element material params from %s\n", i, material_params_filename);

    return FFEA_OK;
}

/*
 */
int Blob::load_stokes_params(const char *stokes_filename, scalar scale) {
    FILE *in = NULL;
    int num_stokes_nodes;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(stokes_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(stokes_filename)
    }
    printf("\t\tReading in material parameters file: %s\n", stokes_filename);

    // first line should be the file type "ffea stokes radii file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of stokes radii file\n")
    }
    if (strcmp(line, "walrus stokes radii file\n") != 0 && strcmp(line, "ffea stokes radii file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea stokes radii file' (read '%s') \n", line)
    }

    // read in the number of nodes in the file
    if (fscanf(in, "num_nodes %d\n", &num_stokes_nodes) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of nodes\n")
    }
    printf("\t\t\tNumber of nodes in stokes radii file = %d\n", num_stokes_nodes);

    // Check that we have same number of nodes in stokes radii file as in nodes file
    if (num_stokes_nodes != num_nodes) {
        fclose(in);
        FFEA_ERROR_MESSG("Number of nodes in stokes radii file (%d) does not match number of nodes in nodes file (%d)\n", num_stokes_nodes, num_nodes)
    }

    // Set the stokes radius for each node in the Blob
    scalar stokes_radius = 0.0;
    int crap, i, check = 0;
    for (i = 0; i < num_nodes; i++) {
        if (fscanf(in, "%le\n", &stokes_radius) != 1) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from stokes radii file at node %d. There should be 1 real value (stokes radius).\n", i);
        }
        node[i].stokes_radius = stokes_radius * scale;
        if (node[i].stokes_radius < 1e-12 && node[i].stokes_radius > 0.0 && check == 0) {
            char finish[1];
            int done = 0;
            printf("WARNING. Stokes Radius on node %d in this Blob is very small, %e.\nStokes Radius is scaled by same factor as node positions, so specify in same units.\n", i, node[i].stokes_radius);
            printf("Would you like to continue (y or n)?:");
            while (done == 0) {
                crap = scanf("%s", finish);
                if (strcmp(finish, "y") == 0) {
                    done = 1;
                    check = 1;
                } else if (strcmp(finish, "n") == 0) {
                    return FFEA_ERROR;
                } else {
                    printf("Please enter y or n:");
                }

            }
        }
    }

    fclose(in);

    printf("\t\t\tRead %d stokes radii from %s\n", i, stokes_filename);
    return FFEA_OK;
}

/**
 * @brief Read the beads_filename, loading beads position and types.
 *
 * @ingroup FMM
 * @details
 */
int Blob::load_beads(const char *beads_filename, PreComp_params *pc_params, scalar scale) {

    ifstream fin;
    string line;
    vector<string> vec_line;
    // int typeBead;

    fin.open(beads_filename, std::ifstream::in);
    if (fin.fail()) {
        FFEA_FILE_ERROR_MESSG(beads_filename);
        return FFEA_ERROR;
    }
    printf("\t\tReading in Beads file: %s\n", beads_filename);
    printf("\t\tScaling beads positions using scale: %e\n", scale);


    vector<string> stypes;
    vector<scalar> positions;
    string type;
    scalar x, y, z;

    // a set of constant strings, and a number of temporary vectors
    //    to parse the lines in search of " < nodes = ... > ".
    const string nodesKeyword = "nodes";
    const char *openField = "<";
    const char *closeField = ">";
    const char *splitNodes = ",";
    const char *defineRange = "-";
    const string defineField = "=";
    vector<string> v1, v2, v3, v4;

    // 1 - read the data, positions and bead-types to memory before storing,
    //       as well as the set of nodes where every bead will be associated to.
    int cnt = 0;
    while (getline(fin, line)) {
        // 1.1 - ignore those lines that do not start with "ATOM"
        if (line.find("ATOM",0,4) != 0)
            continue;

        // cout << line << endl;
        // 1.2 - get bead type and position from its positioning within the line:
        type = line.substr(11,5);
        boost::trim (type);
        stypes.push_back(type);
        type.clear();

        x = stod( line.substr(28,10) ) * scale;
        y = stod( line.substr(38,8) ) * scale;
        z = stod( line.substr(46,8) ) * scale;
        positions.push_back(x);
        positions.push_back(y);
        positions.push_back(z);


        // 1.3 - look for node restrictions within "< nodes = ...  >"
        // 1.3.1 - split the line in a vector using "<":
        boost::split(v1, line, boost::is_any_of(openField));
        bead_assignment.push_back(vector<int>()); // add a row.
        // 1.3.2 - for each of them:
        for (unsigned int j=0; j<v1.size(); j++) {
            // 1.3.3 - remove the closing bracket ">"
            v1[j] = boost::erase_last_copy(v1[j], closeField);
            boost::trim(v1[j]);
            // 1.3.4 - split using "=":
            boost::split(v2, v1[j], boost::is_any_of(defineField));
            boost::trim(v2[0]);
            // 1.3.5 - see if we have "nodes" there:
            if (boost::iequals(nodesKeyword, v2[0])) {
                // cout << " appending nodes " << v2[1] << endl;
                // 1.3.6 - the nodes are comma sepparated:
                boost::split(v3, v2[1], boost::is_any_of(splitNodes));
                // 1.3.7 - and there could be ranges:
                for (unsigned int k=0; k<v3.size(); k++) {
                    boost::split(v4, v3[k], boost::is_any_of(defineRange));
                    if (v4.size() == 1) { // append a single integer:
                        bead_assignment[cnt].push_back(stoi(v4[0]));
                    } else if (v4.size() == 2) {
                        if (stoi(v4[0]) > stoi(v4[1])) {
                            swap(v4[0], v4[1]);
                        }
                        for (int ik = stoi(v4[0]); ik < stoi(v4[1]) + 1; ik++) {
                            bead_assignment[cnt].push_back(ik);
                        }
                        // 1.3.8 - and typos:
                    } else {
                        FFEA_error_text();
                        cout << " failed to parse this set of nodes: " << v2[1] << endl;
                        return FFEA_ERROR;
                    }
                    v4.clear(); // clear temporary vector v4 for the next round.
                }
                v3.clear(); // clear temporary vector v3 for the next round.
            }
            v2.clear(); // clear temporary vector v2 for the next round.
        }
        cnt += 1;
        v1.clear(); // clear temporary vector v1 for the next round.
    }

    // 2 - store the data efficiently:
    // 2.1 - positions:
    bead_position = new(std::nothrow) scalar[positions.size()];
    if (bead_position == NULL) FFEA_ERROR_MESSG("Failed to allocate memory for bead positions\n")
        for (unsigned int i=0; i<positions.size(); i++) {
            bead_position[i] = positions[i];
        }

    // 2.2 - bead types are integers starting from zero:
    vector<string>::iterator it;
    bead_type = new(std::nothrow) int[stypes.size()];
    if (bead_type == NULL) FFEA_ERROR_MESSG("Failed to allocate memory for bead types\n")
        int index;
    for (unsigned int i=0; i<stypes.size(); i++) {
        it = std::find(pc_params->types.begin(), pc_params->types.end(), stypes[i]);
        if (it == pc_params->types.end()) { // type in beads file not matching the types in .ffea file!!
            FFEA_ERROR_MESSG("Type '%s' read in beads file does not match any of the bead types specified in the .ffea file\n", stypes[i].c_str());
        }
        index = std::distance(pc_params->types.begin(), it);
        bead_type[i] = index;
    }

    // 2.3 - num_beads:
    num_beads = stypes.size();
    beads_on_blob = true;

    return FFEA_OK;

}

int Blob::load_ctforces(string ctforces_fname) {

    FFEA_input_reader reader;
    vector<string> ctforces_lines, line_split;
    int n_ctforces = 0;  // total number of forces in the input file.
    int n_ct_lforces = 0; // number of linear forces in the input file
    int n_ct_rforces = 0; // number of rotational forces in the input file
    int n_cts_lforces = 0; // number of forces to be spread on contiguous surface faces.

    // read the input file, taking out comments as <!-- -->
    //  and put it into the ctforces_lines vector of strings:
    reader.file_to_lines(ctforces_fname, &ctforces_lines);


    // 1 - READ AND CHECK THE HEADER:
    //   Check 1:
    if (ctforces_lines.size() < 4) {
        FFEA_ERROR_MESSG("ctforces_fname '%s' is too short, something is wrong\n", ctforces_fname.c_str());
        return FFEA_ERROR;
    }
    // 1.2 - Now read num_ctforces, num_linear_forces, num_rot_forces, num_linear_surface_forces
    int now_reading = 1;
    // WRITE A FOR LOOP TO READ THE HEADER OUT OF ORDER;
    for (now_reading; now_reading < 5; now_reading++) {
        //   Check 2: it ":" appears -> the header is over -> break out.
        if (now_reading == ctforces_lines.size()) break;
        if (ctforces_lines[now_reading].find_first_of(":") != string::npos) break;

        boost::split(line_split, ctforces_lines[now_reading], boost::is_any_of(" \t"), boost::token_compress_on);

        if (line_split[0] == "num_ctforces") {
            try {
                n_ctforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                FFEA_ERROR_MESSG("Invalid number of ctforces read: %s\n", line_split[1].c_str());
                return FFEA_ERROR;
            }
        }

        if (line_split[0] == "num_linear_forces") {
            try {
                n_ct_lforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                FFEA_ERROR_MESSG("Invalid number of constant linear forces read: %s\n", line_split[1].c_str());
                return FFEA_ERROR;
            }
        }

        if (line_split[0] == "num_rot_forces") {
            try {
                n_ct_rforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                FFEA_ERROR_MESSG("Invalid number of circular forces read: %s\n", line_split[1].c_str());
                return FFEA_ERROR;
            }
        }

        if (line_split[0] == "num_linear_surface_forces") {
            try {
                n_cts_lforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                FFEA_ERROR_MESSG("Invalid number of n_cts_lforces read: %s\n", line_split[1].c_str());
                return FFEA_ERROR;
            }
        }
    }
    // 1.2 - Now check that the values are consistent:
    bool Err = false;
    if (n_ctforces != n_ct_lforces + n_ct_rforces + n_cts_lforces) {
        cout << "--- n_ctforces != n_ct_lforces + n_ct_rforces + n_cts_lforces" << endl;
        Err = true;
    } else if (n_ctforces == 0) Err = true;
    int calc_n_lines = n_ctforces + now_reading;
    if (n_ct_lforces) calc_n_lines += 1;
    if (n_ct_rforces) calc_n_lines += 1;
    if (n_cts_lforces) calc_n_lines += 1;
    if (calc_n_lines != ctforces_lines.size()) {
        Err = true;
        cout << "--- The total number of lines in the ctforces file should be " << calc_n_lines
             << " but is actually " << ctforces_lines.size() << endl;
    }

    if (Err) {
        cout << "--- ABORTING. Something went wrong when parsing the header of the ctforces file" << endl;
        cout << "--- Check it or turn calc_ctforces to 0 in the FFEA input file" << endl;
        return FFEA_ERROR;
    }


    // 2 - READ AND STORE THE LINEAR PART:
    // 2.1 - First the header:
    // check that we have a line saying "linear forces:"
    if (ctforces_lines[now_reading].compare(0, 14, "linear forces:") != 0) {
        // if not, and needed: ABORT.
        if (n_ct_lforces > 0) {
            FFEA_ERROR_MESSG("Error reading the ctforces file; it should announce the start of linear ctforces, but instead read: %s\n", ctforces_lines[now_reading].c_str());
            return FFEA_ERROR;
        }
    } else {
        now_reading += 1;
    }

    // 2.2 - Now put the lines referring to this blob/conf into a vector,
    //   and calculate the number of constant forces we're adding:
    vector<string> my_lines;
    for (int i=now_reading; i<now_reading+n_ct_lforces; i++) {
        boost::split(line_split, ctforces_lines[i], boost::is_any_of("  \t"), boost::token_compress_on);
        // at least check the length of the line:
        if (line_split.size() != 7) {
            FFEA_ERROR_MESSG("Invalid line in the ctforces file:\n \t %s\n", ctforces_lines[i].c_str());
            return FFEA_ERROR;
        }
        int b_i = boost::lexical_cast<int>(line_split[4]); // read the blob index
        if (b_i != blob_index) continue;
        int c_i = boost::lexical_cast<int>(line_split[5]); // read the conformation index
        if (c_i != conformation_index) continue;
        if (line_split[6].compare("all") == 0) {
            num_l_ctf += num_nodes;
        } else {
            num_l_ctf += 1;
        }
        my_lines.push_back(ctforces_lines[i]);
    }

    if (num_l_ctf > 0) cout << num_l_ctf << " linear forces were loaded for blob " << blob_index << ":" << conformation_index << endl;

    // 2.3 - And finally store the stuff properly:
    ctf_l_nodes = new(std::nothrow) int[num_l_ctf];       // allocate nodes
    ctf_l_forces = new(std::nothrow) scalar[3*num_l_ctf]; // allocate forces
    if (ctf_l_nodes == NULL || ctf_l_forces == NULL) FFEA_ERROR_MESSG("Failed to allocate ctf_l relevant arrays\n");
    arr3 ctf_d; // direction of the force
    int cnt = 0;
    for (int i=0; i<my_lines.size(); i++) {
        boost::split(line_split, my_lines[i], boost::is_any_of(" \t"), boost::token_compress_on);
        scalar F = boost::lexical_cast<scalar>(line_split[0]);
        ctf_d[0] = boost::lexical_cast<scalar>(line_split[1]);
        ctf_d[1] = boost::lexical_cast<scalar>(line_split[2]);
        ctf_d[2] = boost::lexical_cast<scalar>(line_split[3]);
        arr3Normalise<scalar,arr3>(ctf_d);  // normalise the direction of the force.
        scalar Fx = F * ctf_d[0] / mesoDimensions::force;
        scalar Fy = F * ctf_d[1] / mesoDimensions::force;
        scalar Fz = F * ctf_d[2] / mesoDimensions::force;
        if (line_split[6].compare("all") != 0) {
            ctf_l_nodes[cnt] = boost::lexical_cast<int>(line_split[6]);
            ctf_l_forces[3*cnt   ]  = Fx;
            ctf_l_forces[3*cnt +1]  = Fy;
            ctf_l_forces[3*cnt +2]  = Fz;
            cnt += 1;
        } else {
            for (int j=0; j<num_nodes; j++) {
                ctf_l_nodes[cnt] = j;
                ctf_l_forces[3*cnt   ]  = Fx;
                ctf_l_forces[3*cnt +1]  = Fy;
                ctf_l_forces[3*cnt +2]  = Fz;
                cnt += 1;
            }
        }
    }



    // 3 - READ AND STORE THE ROTATIONAL PART:
    // 3.1 - First the header:
    now_reading = now_reading + n_ct_lforces;
    // check that we're still within bounds:
    if (now_reading >= ctforces_lines.size()) { // Out of bounds!
        if (n_ct_rforces + n_cts_lforces > 0) { // And we should be in!!!
            FFEA_ERROR_MESSG("Wrong number of lines in ctforces. ABORTING");
            return FFEA_ERROR;
        } else return FFEA_OK; // Oh, the work was over.
    }
    // check that we have a line saying "rotational forces:"
    if (ctforces_lines[now_reading].compare(0, 18, "rotational forces:") != 0) {
        // if not, and needed: ABORT.
        if (n_ct_rforces > 0) {
            FFEA_ERROR_MESSG("Wrong header; it should announce the start of rotational ctforces, but instead read: %s\n", ctforces_lines[now_reading].c_str());
            return FFEA_ERROR;
        }
    } else {
        now_reading += 1;
        cout << " loading rotational forces for blob " << blob_index << ":" << conformation_index << endl;
    }

    // 3.2 - Now put the lines referring to this blob/conf into a vector,
    //   and calculate the number of rot constant forces we're adding:
    my_lines.clear();
    for (int i=now_reading; i<now_reading+n_ct_rforces; i++) {
        boost::split(line_split, ctforces_lines[i], boost::is_any_of("  \t"), boost::token_compress_on);
        // at least check the length of the line:
        if (line_split.size() != 11) { // check length:
            FFEA_ERROR_MESSG("Invalid line in the ctforces file:\n \t %s\n", ctforces_lines[i].c_str());
            return FFEA_ERROR;
        }
        int b_i = boost::lexical_cast<int>(line_split[8]); // read the blob index
        if (b_i != blob_index) continue;
        int c_i = boost::lexical_cast<int>(line_split[9]); // read the conformation index
        if (c_i != conformation_index) continue;
        if (line_split[10].compare("all") == 0) {
            num_r_ctf += num_nodes;
        } else {
            num_r_ctf += 1;
        }
        my_lines.push_back(ctforces_lines[i]);
    }

    // 3.3 - And finally store the stuff properly:
    arr3 ctf_p; // temporary point in the axis.
    ctf_r_nodes = new(std::nothrow) int[num_r_ctf];       // allocate nodes
    ctf_r_forces = new(std::nothrow) scalar[num_r_ctf]; // allocate forces
    ctf_r_axis = new(std::nothrow) scalar[6*num_r_ctf]; // allocate axis
    ctf_r_type = new(std::nothrow) char[2*num_r_ctf]; // allocate type of rotational force
    if (ctf_r_nodes == NULL || ctf_r_forces == NULL || ctf_r_axis == NULL || ctf_r_type == NULL) FFEA_ERROR_MESSG("Failed to allocate ctf_r relevant arrays\n");
    const scalar mdfm1 = 1./mesoDimensions::force;
    const scalar mdlm1 = 1./mesoDimensions::length;
    cnt = 0; // reinitialise cnt
    for (int i=0; i<my_lines.size(); i++) {
        boost::split(line_split, my_lines[i], boost::is_any_of(" \t"), boost::token_compress_on);
        scalar F = boost::lexical_cast<scalar>(line_split[0]) * mdfm1;
        string type = line_split[1];
        ctf_p[0] = boost::lexical_cast<scalar>(line_split[2]);
        ctf_p[1] = boost::lexical_cast<scalar>(line_split[3]);
        ctf_p[2] = boost::lexical_cast<scalar>(line_split[4]);
        ctf_d[0] = boost::lexical_cast<scalar>(line_split[5]);
        ctf_d[1] = boost::lexical_cast<scalar>(line_split[6]);
        ctf_d[2] = boost::lexical_cast<scalar>(line_split[7]);
        if (type.compare(0,1,"p") == 0) {    // store as point + direction:
            arr3Normalise<scalar,arr3>(ctf_d);  //  and thus normalise.
            arr3Resize<scalar,arr3>(mdlm1, ctf_p);  // and rescale CTFPENDING: check!!!
        } else if (type.compare(0,1,"n")) { // otherwise store as pairs of nodes, or complain.
            FFEA_ERROR_MESSG("Invalid rotational force: %s, in line read: %s\n", type.substr(0).c_str(), my_lines[i].c_str());
            return FFEA_ERROR;
        }
        if ((type.compare(1,1,"f")) and (type.compare(1,1,"t"))) { // check and store type force or torque
            FFEA_ERROR_MESSG("Invalid rotational force type: %s, in line read: %s\n", type.substr(1).c_str(), my_lines[i].c_str());
            return FFEA_ERROR;
        }
        if (line_split[10].compare("all") != 0) {
            ctf_r_nodes[cnt] = boost::lexical_cast<int>(line_split[10]);
            ctf_r_forces[cnt]  = F;
            ctf_r_axis[6*cnt    ] = ctf_p[0];
            ctf_r_axis[6*cnt + 1] = ctf_p[1];
            ctf_r_axis[6*cnt + 2] = ctf_p[2];
            ctf_r_axis[6*cnt + 3] = ctf_d[0];
            ctf_r_axis[6*cnt + 4] = ctf_d[1];
            ctf_r_axis[6*cnt + 5] = ctf_d[2];
            ctf_r_type[2*cnt    ] = type[0];
            ctf_r_type[2*cnt + 1] = type[1];
            cnt += 1;
        } else {
            for (int j=0; j<num_nodes; j++) {
                ctf_r_nodes[cnt] = j;
                ctf_r_forces[cnt]  = F;
                ctf_r_axis[6*cnt    ] = ctf_p[0];
                ctf_r_axis[6*cnt + 1] = ctf_p[1];
                ctf_r_axis[6*cnt + 2] = ctf_p[2];
                ctf_r_axis[6*cnt + 3] = ctf_d[0];
                ctf_r_axis[6*cnt + 4] = ctf_d[1];
                ctf_r_axis[6*cnt + 5] = ctf_d[2];
                ctf_r_type[2*cnt    ] = type[0];
                ctf_r_type[2*cnt + 1] = type[1];
                cnt += 1;
            }
        }
    }



    // 4 - READ AND STORE THE SURFACE PART:
    // 4.1 - First the header:
    now_reading = now_reading + n_ct_rforces;
    // check that we're still within bounds:
    if (now_reading >= ctforces_lines.size()) { // Out of bounds!
        if (n_cts_lforces > 0) { // And we should be in!!!
            FFEA_ERROR_MESSG("Wrong number of lines in ctforces. ABORTING");
            return FFEA_ERROR;
        } else return FFEA_OK; // Oh, the work was over.
    }
    // check that we have a line saying "linear surface forces:"
    if (ctforces_lines[now_reading].compare(0, 22, "linear surface forces:") != 0) {
        // if not, and needed: ABORT.
        if (n_cts_lforces > 0) {
            FFEA_ERROR_MESSG("Wrong header; it should announce the start of the linear surface ctforces, but instead read: %s\n", ctforces_lines[now_reading].c_str());
            return FFEA_ERROR;
        }
    } else {
        now_reading += 1;
        cout << " loading linear surface forces for blob " << blob_index << ":" << conformation_index << endl;
    }

    // 4.2 - Now put the lines referring to this blob/conf into a vector,
    //   and calculate the number of rot constant forces we're adding:
    my_lines.clear();
    for (int i=now_reading; i<now_reading+n_cts_lforces; i++) {
        boost::split(line_split, ctforces_lines[i], boost::is_any_of("  \t"), boost::token_compress_on);
        // at least check the minimum length of the line:
        if (line_split.size() < 7) {
            FFEA_ERROR_MESSG("Invalid line in the ctforces file:\n \t %s\n", ctforces_lines[i].c_str());
            return FFEA_ERROR;
        }
        int b_i = boost::lexical_cast<int>(line_split[4]); // read the blob index
        if (b_i != blob_index) continue;
        int c_i = boost::lexical_cast<int>(line_split[5]); // read the conformation index
        if (c_i != conformation_index) continue;
        if (line_split[6].compare("all") == 0) {
            num_sltotal_ctf += num_surface_faces;
            return FFEA_ERROR;
        } else {
            num_sltotal_ctf += line_split.size() - 6;  // there may be many faces in every 'surface'!
        }
        my_lines.push_back(ctforces_lines[i]);
    }

    // 4.3 - And finally store the stuff properly:
    //     - If a face were appearing different lines, force will be applied "serially",
    //            and not summed up. They may be part of different "surfaces".
    ctf_sl_faces = new(std::nothrow) int[num_sltotal_ctf];       // allocate faces
    num_slsets_ctf = my_lines.size();
    ctf_sl_surfsize = new(std::nothrow) int[num_slsets_ctf];
    ctf_sl_forces = new(std::nothrow) scalar[3*num_slsets_ctf]; // allocate forces
    if (ctf_sl_faces == NULL || ctf_sl_surfsize == NULL || ctf_sl_forces == NULL) FFEA_ERROR_MESSG("Failed to allocate ctf_sl relevant arrays\n");
    ctf_d; // direction of the force
    cnt = 0;
    for (int i=0; i<my_lines.size(); i++) {
        boost::split(line_split, my_lines[i], boost::is_any_of(" \t"), boost::token_compress_on);
        scalar F = boost::lexical_cast<scalar>(line_split[0]);
        ctf_d[0] = boost::lexical_cast<scalar>(line_split[1]);
        ctf_d[1] = boost::lexical_cast<scalar>(line_split[2]);
        ctf_d[2] = boost::lexical_cast<scalar>(line_split[3]);
        arr3Normalise<scalar,arr3>(ctf_d);  // normalise the direction of the force.
        ctf_sl_forces[3*i   ]  = F * ctf_d[0] / mesoDimensions::force;
        ctf_sl_forces[3*i +1]  = F * ctf_d[1] / mesoDimensions::force;
        ctf_sl_forces[3*i +2]  = F * ctf_d[2] / mesoDimensions::force;
        if (line_split[6].compare("all") != 0) {
            ctf_sl_surfsize[i] = line_split.size() - 6;
            for (int j=6; j<line_split.size(); j++) {
                ctf_sl_faces[cnt] = boost::lexical_cast<int>(line_split[j]);
                cnt += 1;
            }
        } else {
            ctf_sl_surfsize[i] = num_surface_faces;
            for (int j=0; j<num_surface_faces; j++) {
                ctf_sl_faces[cnt] = j;
                cnt += 1;
            }
        }
    }

    /* // SELF-EXPLANATORY CHECK //
    int auxndx = 0;
    for (int i=0; i<num_slsets_ctf; i++) {
      cout << "set " << i << ": ";
      cout << "force: " << ctf_sl_forces[3*i] << ", " << ctf_sl_forces[3*i+1] << ", " << ctf_sl_forces[3*i+2] << " affecting " << ctf_sl_surfsize[i] << " faces. The faces being: ";
      for (int j=0; j<ctf_sl_surfsize[i]; j++) {
        cout << ctf_sl_faces[auxndx + j] << " " ;
      }
      cout << endl;
      auxndx += ctf_sl_surfsize[i];
    } */

    return FFEA_OK;
}


void Blob::add_steric_nodes() {

    int i;
    for(i = 0; i < num_surface_faces; ++i) {
        surface[i].build_opposite_node();
    }
}

/**
 * @brief num_beads = 0; delete bead_type; delete bead_position.
 *
 * @ingroup FMM
 * @details Beads are only useful before PreComp_solver.init is called.
 *      * They can be removed later on.
 */
int Blob::forget_beads() {

    num_beads = 0;
    delete[] bead_position;
    bead_position = NULL;
    delete[] bead_type;
    bead_type = NULL;
    return FFEA_OK;

}


void Blob::print_node_positions() {

    for (int n = 0; n < num_nodes; n++) {
        cout << "---n: " << node[n].pos.x << " " << node[n].pos.y << "  " << node[n].pos.z << endl;
    }

}

void Blob::print_bead_positions() {

    for (int i=0; i<num_beads; i++) {
        cout << "---b: " << bead_position[3*i] << " "
             << bead_position[3*i+1] << " "
             << bead_position[3*i+2] << endl;
    }

}

/*
 */
int Blob::load_ssint(const char *ssint_filename, int num_ssint_face_types, string ssint_method) {
    FILE *in = NULL;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(ssint_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(ssint_filename)
    }
    printf("\t\tReading in VDW file: %s\n", ssint_filename);

    // first line should be the file type "ffea vdw file" (.vdw) or "ffea ssint file" (.ssint)
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of VDW file\n")
    }
    if (strcmp(line, "walrus vdw file\n") != 0 && strcmp(line, "ffea vdw file\n") != 0 && strcmp(line, "ffea ssint file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea vdw file' (read '%s') \n", line)
    }

    // read in the number of faces in the file
    int num_ssint_faces = 0;
    if (fscanf(in, "num_faces %d\n", &num_ssint_faces) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of faces\n")
    }
    printf("\t\t\tNumber of faces = %d\n", num_ssint_faces);

    if (num_ssint_faces != num_surface_faces) {
        FFEA_ERROR_MESSG("Number of faces specified in VDW file (%d) does not agree with number in surface file (%d)\n", num_ssint_faces, num_surface_faces)
    }

    // Check for "ssint params:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error when looking for 'vdw params:' line\n")
    }
    if (strcmp(line, "vdw params:\n") != 0 && strcmp(line, "ssint params:\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("Could not find 'vdw params:' line (found '%s' instead)\n", line)
    }

    // Read in all the ssint parameters from the file, assigning them to the appropriate faces
    int i;
    int ssint_type = 0;

    // If steric only, set all to type 0
    if (params->calc_ssint == 0) {
        for(i = 0; i < num_surface_faces; ++i) {
            if (fscanf(in, "%d\n", &ssint_type) != 1) {
                fclose(in);
                FFEA_ERROR_MESSG("Error reading from VDW file at face %d. There should be 1 integer denoting ssint face species (-1 - unreactive). \n", i);
            } else {
                if (ssint_type > num_ssint_face_types - 1) {
                    ssint_type = 0;
                }
                surface[i].set_ssint_interaction_type(ssint_type);
            }
        }
    } else {
        for (i = 0; i < num_surface_faces; i++) {
            if (fscanf(in, "%d\n", &ssint_type) != 1) {
                fclose(in);
                FFEA_ERROR_MESSG("Error reading from VDW file at face %d. There should be 1 integer denoting ssint face species (-1 - unreactive). \n", i);
            } else {
                if (ssint_type > num_ssint_face_types - 1) {
                    FFEA_ERROR_MESSG("Error reading from VDW file at face %d. The given VDW face type (%d) is higher than that allowed by the ssint forcefield params file (%d). \n", i, ssint_type, num_ssint_face_types - 1);
                }
                surface[i].set_ssint_interaction_type(ssint_type);
            }
        }
    }

    // Set whether vdw is active on blob
    for(i = 0; i < num_surface_faces; ++i) {
        if(surface[i].ssint_interaction_type != -1) {
            ssint_on_blob = true;
            break;
        }
    }
    fclose(in);

    printf("\t\t\tRead %d VDW faces from %s\n", i, ssint_filename);

    return FFEA_OK;
}

/*
 */
int Blob::load_binding_sites() {

    int num_binding_site_types = binding_matrix->get_num_interaction_types();

    // Return successful as params.calc_kinetics == 0 or no sites are required
    if (s_binding_filename.empty()) return FFEA_OK;
    // Open file
    ifstream fin;
    fin.open(s_binding_filename.c_str());
    if(fin.fail()) {
        FFEA_ERROR_MESSG("'binding_params_fname' %s not found\n", s_binding_filename.c_str())
    }

    cout << "\t\tReading in Binding Sites file: " << s_binding_filename << endl;

    // Check if correct file
    int MAX_BUF_SIZE = 255;
    char buf[MAX_BUF_SIZE];
    string buf_string;
    vector<string> string_vec;
    fin.getline(buf, MAX_BUF_SIZE);
    buf_string = string(buf);
    boost::trim(buf_string);

    if(buf_string != "ffea binding sites file") {
        FFEA_ERROR_MESSG("This is not a 'ffea binding site file' (read '%s') \n", buf)
    }

    // read in the number of binding sites in the file
    fin >> buf_string >> num_binding_sites;
    cout << "\t\t\tNumber of binding sites = " << num_binding_sites << endl;

    if (num_binding_sites > num_surface_faces) {
        FFEA_ERROR_MESSG("Number of binding sites specified in binding sites file (%d) cannot exceed number of surface faces (%d)\n", num_binding_sites, num_surface_faces)
    }

    if (num_binding_sites == 0) {
        return FFEA_OK;
    }

    // Create binding sites
    binding_site = new(std::nothrow) BindingSite[num_binding_sites];
    if (binding_site == NULL) FFEA_ERROR_MESSG("Failed to allocate array of binding sites\n");

    // Check for "binding sites:" line
    fin.getline(buf, MAX_BUF_SIZE);
    fin.getline(buf, MAX_BUF_SIZE);
    buf_string = string(buf);
    boost::trim(buf_string);
    if(buf_string != "binding sites:") {
        FFEA_ERROR_MESSG("Could not find 'binding sites:' line (found '%s' instead)\n", buf)
    }

    // Get all binding sites
    unsigned int num_faces = 0;
    int bind_type = -1, face_index;
    for(int i = 0; i < num_binding_sites; ++i) {

        // Get structural details first
        fin.getline(buf, MAX_BUF_SIZE);
        buf_string = string(buf);
        boost::trim(buf_string);
        boost::split(string_vec, buf_string, boost::is_space(), boost::token_compress_on);
        try {
            bind_type = atoi(string_vec.at(1).c_str());
            num_faces = atoi(string_vec.at(3).c_str());

        } catch (...) {
            FFEA_ERROR_MESSG("Unable to read type %%d num_faces %%d line for binding site %d in %s.\n", i, s_binding_filename.c_str())
        }

        if(bind_type >= num_binding_site_types) {
            FFEA_ERROR_MESSG("Binding site %d specifies site type %d, which is outside range of types allowed by the 'bsite_in_fname' matrix (%d types allowed)\n", i, bind_type, num_binding_site_types)
            return FFEA_ERROR;
        }

        binding_site[i].set_type(bind_type);
        binding_site[i].set_num_faces(num_faces);

        // Now build list of faces
        fin.getline(buf, MAX_BUF_SIZE);
        buf_string = string(buf);
        boost::trim(buf_string);
        boost::split(string_vec, buf_string, boost::is_space(), boost::token_compress_on);
        if(string_vec.size() != num_faces + 1) {
            FFEA_ERROR_MESSG("In %s, num_faces specified, %d, != num_faces in following line, %zd.\n", s_binding_filename.c_str(), num_faces, string_vec.size() - 1)
        }

        for(unsigned int j = 0; j < num_faces; ++j) {
            face_index = atoi(string_vec.at(j + 1).c_str());
            if(face_index >= num_surface_faces) {
                FFEA_ERROR_MESSG("Face index %d specifies face outside range of surface faces defined in surface file (%d)\n", face_index, num_surface_faces)
                return FFEA_ERROR;
            } else {
                binding_site[i].add_face(&surface[face_index]);

            }
        }

        // Properties continually change so no need to calculate stuff unless about to be used
    }

    fin.close();
    return FFEA_OK;
}

/*
 */
int Blob::load_pinned_nodes(const char *pin_filename) {
    FILE *in = NULL;
    int i, pn_index;
    const int max_line_size = 50;
    char line[max_line_size];

    if ((in = fopen(pin_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(pin_filename)
    }
    printf("\t\tReading in pinned nodes file: %s\n", pin_filename);

    // first line should be the file type "ffea pinned nodes file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of pin file\n")
    }
    if (strcmp(line, "walrus pinned nodes file\n") != 0 && strcmp(line, "ffea pinned nodes file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea pinned nodes file' (read '%s') \n", line)
    }

    // read in the number of pinned node indices in the file
    if (fscanf(in, "num_pinned_nodes %d\n", &num_pinned_nodes) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of pinned nodes\n")
    }
    printf("\t\t\tNumber of pinned nodes = %d\n", num_pinned_nodes);

    // Allocate the memory for the list of pinned node indices
    pinned_nodes_list = new int[num_pinned_nodes];
    if (pinned_nodes_list == NULL) FFEA_ERROR_MESSG("Failed to allocate pinned_nodes_list\n");

    // Check for "pinned nodes:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error when looking for 'pinned nodes:' line\n")
    }
    if (strcmp(line, "pinned nodes:\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("Could not find 'pinned nodes:' line (found '%s' instead)\n", line)
    }

    // Read in all the pinned node indices from file
    for (i = 0; i < num_pinned_nodes; i++) {
        if (fscanf(in, "%d\n", &pn_index) != 1) {
            fclose(in);
            FFEA_ERROR_MESSG("Error reading from pinned node file at face %d. There should be 1 integer per line. \n", i);
        } else {
            // check that this does not reference nodes outside of the node array
            if (pn_index < 0 || pn_index >= num_nodes) {
                fclose(in);
                FFEA_ERROR_MESSG("Error: Pinned node %d references an out of bounds node index\n", i);
            }
            pinned_nodes_list[i] = pn_index;
        }
    }

    fclose(in);

    if (i == 1)
        printf("\t\t\tRead 1 pinned node index from %s\n", pin_filename);
    else
        printf("\t\t\tRead %d pinned node indices from %s\n", i, pin_filename);

    return FFEA_OK;
}

/*
 */
void Blob::calc_rest_state_info() {
    // Calculate some rest state information for each element:
    // Calculate the inverse jacobian (needed for gradient deformation tensor
    // calc.) and use the determinant to calculate the element's rest volume
    // (needed for volumetric spring calc.). Also, get the diffusion matrix info
    // for constructing the preliminary poisson solver matrix.
    matrix3 J;
    scalar min_vol = INFINITY, temp;
    scalar longest_edge = 0;
    scalar longest_surface_edge = 0;
    int min_vol_elem = 0;
    mass = 0;
    scalar total_vol = 0;
    for (int i = 0; i < num_elements; i++) {
        // Get jacobian matrix for this element
        elem[i].calculate_jacobian(J);

        // Get the inverse jacobian matrix
        mat3_invert(J, elem[i].J_inv_0, &temp);

        // get the 12 derivatives of the shape functions (by inverting the jacobian)
        // and also get the element volume
        elem[i].calc_shape_function_derivatives_and_volume(J);
        elem[i].vol_0 = elem[i].vol;

        total_vol += elem[i].vol_0;

        // Keep a check on which element has the smallest volume
        if (elem[i].vol_0 < min_vol) {
            min_vol = elem[i].vol_0;
            min_vol_elem = i;
        }

        scalar longest_edge_i = elem[i].length_of_longest_edge();
        if (longest_edge_i > longest_edge) longest_edge = longest_edge_i;

        // Calc the mass of the element
        elem[i].mass = elem[i].vol_0 * elem[i].rho;
    }

    for (int i=0; i<num_surface_faces; i++) {
        scalar longest_surface_edge_i = surface[i].length_of_longest_edge();
        if (longest_surface_edge < longest_surface_edge_i) longest_surface_edge = longest_surface_edge_i;
    }

    if (blob_state == FFEA_BLOB_IS_STATIC) {
        printf("\t\tBlob is static, so volume not defined within simulation.\n");
        printf("\t\tDefining Total rest volume of Blob to be 0 cubic Angstroms.\n");
        printf("\t\tAll elements have volume 0 cubic Angstroms.\n");
        min_vol = 0.0;
    } else {
        printf("\t\tTotal rest volume of Blob is %e cubic Angstroms.\n", total_vol * mesoDimensions::volume* 1e30);
        printf("\t\tSmallest element (%i) has volume %e cubic Angstroms.\n", min_vol_elem, min_vol * mesoDimensions::volume * 1e30);
        printf("\t\tLongest edge has length %e Angstroms.\n", longest_edge * mesoDimensions::length * 1e10);
        printf("\t\tLongest surface edge has length %e Angstroms.\n", longest_surface_edge * mesoDimensions::length * 1e10);
    }

    // Calc the total mass of this Blob
    for (int i = 0; i < num_elements; i++) {
        mass += elem[i].mass;
    }
    printf("\t\tTotal mass of blob = %e\n", mass*mesoDimensions::mass);
}

/*
 */
int Blob::aggregate_forces_and_solve() {
    int n, m;

    // Aggregate the forces on each node by summing the contributions from each element.
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) private(n, m) schedule(guided)
#endif
    for (n = 0; n < num_nodes; n++) {
        for (m = 0; m < node[n].num_element_contributors; m++) {
            force[n][0] += node[n].force_contributions[m]->x;
            force[n][1] += node[n].force_contributions[m]->y;
            force[n][2] += node[n].force_contributions[m]->z;
        }
    }

    // Aggregate surface forces onto nodes
    for (n = 0; n < num_surface_faces; n++) {
        for (int i = 0; i < 4; i++) {
            int sni = surface[n].n[i]->index;
            force[sni].x += surface[n].force[i].x;
            force[sni].y += surface[n].force[i].y;
            force[sni].z += surface[n].force[i].z;

            //			printf("force on %d from face %d = %e %e %e\n", sni, n, force[sni].x, force[sni].y, force[sni].z);
        }
    }
    //	printf("----\n\n");
    if (params->calc_stokes == 1) {
        if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            #pragma omp parallel default(none)
            {
#endif
#ifdef USE_OPENMP
                int thread_id = omp_get_thread_num();
#else
                int thread_id = 0;
#endif
#ifdef FFEA_PARALLEL_WITHIN_BLOB
                #pragma omp for schedule(guided)
#endif
                for (int i = 0; i < num_nodes; i++) {
                    force[i].x -= node[i].vel.x * node[i].stokes_drag;
                    force[i].y -= node[i].vel.y * node[i].stokes_drag;
                    force[i].z -= node[i].vel.z * node[i].stokes_drag;
                    if (params->calc_noise == 1) {
                        force[i].x -= RAND(-.5, .5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
                        force[i].y -= RAND(-.5, .5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
                        force[i].z -= RAND(-.5, .5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
                    }
                }
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            }
#endif
        } else {
            if (params->calc_noise == 1) {

#ifdef FFEA_PARALLEL_WITHIN_BLOB
                #pragma omp parallel default(none)
                {
#endif
#ifdef USE_OPENMP
                    int thread_id = omp_get_thread_num();
#else
                    int thread_id = 0;
#endif
#ifdef FFEA_PARALLEL_WITHIN_BLOB
                    #pragma omp for schedule(guided)
#endif
                    for (int i = 0; i < num_nodes; i++) {
                        force[i].x -= RAND(-.5, .5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
                        force[i].y -= RAND(-.5, .5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
                        force[i].z -= RAND(-.5, .5) * sqrt((24 * params->kT * node[i].stokes_drag) / (params->dt));
                    }
#ifdef FFEA_PARALLEL_WITHIN_BLOB
                }
#endif
            }
        }
    }

    // Take forces off of second order nodes onto linear ones
    linearise_force();

    // Set to zero any forces on the pinned nodes
    for (n = 0; n < num_pinned_nodes; n++) {
        int pn_index = pinned_nodes_list[n];
        force[pn_index].x = 0;
        force[pn_index].y = 0;
        force[pn_index].z = 0;
    }

    // This should have a similar way of making the solver matrix become the identity matrix as the pinned nodes do
    for(set<int>::iterator it = bsite_pinned_nodes_list.begin(); it != bsite_pinned_nodes_list.end(); ++it) {
        force[*it].x = 0;
        force[*it].y = 0;
        force[*it].z = 0;
    }

    // Use the linear solver to solve for Mx = f where M is the Blob's mass matrix,
    // or Kv = f where K is the viscosity matrix for the system
    // x/v is the (unknown) force solution and f is the force vector for the system.
    if (solver->solve(force) == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Error reported by Solver.\n");
    }
    return FFEA_OK;
}

/*
 *
 */
void Blob::euler_integrate() {
    int i;

    // Update the velocities and positions of all the nodes
    if (linear_solver == FFEA_NOMASS_CG_SOLVER) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none) private(i) schedule(static)
#endif
        for (i = 0; i < num_nodes; i++) {
            node[i].vel[0] = force[i][0];
            node[i].vel[1] = force[i][1];
            node[i].vel[2] = force[i][2];

            node[i].pos[0] += force[i][0] * params->dt; // really meaning v * dt
            node[i].pos[1] += force[i][1] * params->dt;
            node[i].pos[2] += force[i][2] * params->dt;
        }

    } else {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none) private(i) schedule(static)
#endif
        for (i = 0; i < num_nodes; i++) {
            node[i].vel[0] += force[i][0] * params->dt;
            node[i].vel[1] += force[i][1] * params->dt;
            node[i].vel[2] += force[i][2] * params->dt;

            node[i].pos[0] += node[i].vel[0] * params->dt;
            node[i].pos[1] += node[i].vel[1] * params->dt;
            node[i].pos[2] += node[i].vel[2] * params->dt;
        }
    }
}

/*
 *
 */
int Blob::calculate_node_element_connectivity() {
    int i, j;
    int *node_counter = NULL;

    // initialise num_element_contributors to zero
    for (i = 0; i < num_nodes; i++)
        node[i].num_element_contributors = 0;

    // count how many times each node is referenced in the list of elements
    for (i = 0; i < num_elements; i++)
        for (j = 0; j < NUM_NODES_QUADRATIC_TET; j++)
            node[elem[i].n[j]->index].num_element_contributors++;

    // allocate the contributions array for each node to the length given by num_element_contributors
    for (i = 0; i < num_nodes; i++) {
        node[i].force_contributions = new(std::nothrow) vector3 *[node[i].num_element_contributors];
        if (node[i].force_contributions == NULL) {
            FFEA_ERROR_MESSG("Failed to allocate memory for 'force_contributions' array (on node %d)\n", i);
        }
    }

    // create an array of counters keeping track of how full each contributions array is on each node
    node_counter = new(std::nothrow) int[num_nodes];
    if (node_counter == NULL) {
        FFEA_ERROR_MESSG("Failed to allocate node_counter array.\n");
    }

    // initialise the node_counter array to zero
    for (i = 0; i < num_nodes; i++)
        node_counter[i] = 0;

    // go back through the elements array and fill the contributions arrays with pointers to the
    // appropriate force contributions in the elements
    int node_index;
    for (i = 0; i < num_elements; i++)
        for (j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
            node_index = elem[i].n[j]->index;
            node[node_index].force_contributions[node_counter[node_index]] = &elem[i].node_force[j];
            node_counter[node_index]++;
        }

    // release node counter array
    delete[] node_counter;

    return FFEA_OK;
}

void Blob::pin_binding_site(set<int> node_indices) {

    set<int>::iterator it;
    for(it = node_indices.begin(); it != node_indices.end(); ++it) {
        bsite_pinned_nodes_list.insert(*it);
    }
}

void Blob::unpin_binding_site(set<int> node_indices) {

    set<int>::iterator it, it2;
    for(it = node_indices.begin(); it != node_indices.end(); ++it) {
        bsite_pinned_nodes_list.erase(*it);
    }
}

int Blob::create_pinned_nodes(set<int> list) {

    int i;
    set<int>::iterator it;
    delete[] pinned_nodes_list;

    pinned_nodes_list = new(std::nothrow) int[list.size()];
    if (pinned_nodes_list == NULL) FFEA_ERROR_MESSG("Could not allocate memory for pinned_nodes_list\n");
    i = 0;
    for(it = list.begin(); it != list.end(); ++it) {
        pinned_nodes_list[i++] = *it;
    }

    return FFEA_OK;
}

int Blob::get_state_index() {

    return state_index;
}

void Blob::set_state_index(int index) {

    this->state_index = index;
}

int Blob::get_previous_state_index() {

    return previous_state_index;
}

void Blob::set_previous_state_index(int index) {

    this->previous_state_index = index;
}

int Blob::get_conformation_index() {

    return conformation_index;
}

int Blob::get_previous_conformation_index() {

    return previous_conformation_index;
}

void Blob::set_previous_conformation_index(int index) {

    this->previous_conformation_index = index;
}

BindingSite* Blob::get_binding_site(int index) {

    return &binding_site[index];
}

int Blob::build_mass_matrix() {
    // Calculate the Sparsity Pattern for the Mass matrix
    printf("\t\tCalculating sparsity pattern for 2nd order Mass matrix\n");
    SparsityPattern sparsity_pattern_mass_matrix;
    sparsity_pattern_mass_matrix.init(num_nodes);

    MassMatrixQuadratic *M_alpha = new(std::nothrow) MassMatrixQuadratic[num_elements];
    if (M_alpha == NULL) FFEA_ERROR_MESSG("Failed to allocate memory for the mass matrix\n");

    scalar *mem_loc;
    int ni_index, nj_index;
    for (int el = 0; el < num_elements; el++) {

        elem[el].construct_element_mass_matrix(&M_alpha[el]);

        for (int ni = 0; ni < 10; ni++) {
            for (int nj = 0; nj < 10; nj++) {

                ni_index = elem[el].n[ni]->index;
                nj_index = elem[el].n[nj]->index;

                mem_loc = M_alpha[el].get_M_alpha_mem_loc(ni, nj);

                sparsity_pattern_mass_matrix.register_contribution(ni_index, nj_index, mem_loc);
            }
        }
    }

    printf("\t\tBuilding sparsity pattern\n");

    // Use the sparsity patterns to create a fixed pattern sparse matrix
    M = sparsity_pattern_mass_matrix.create_sparse_matrix();

    // Build the mass matrix
    M->build();

    return FFEA_OK;
}

bool Blob::there_is_mass() {
    return mass_in_blob;
}

void Blob::set_springs_on_blob(bool state) {
    springs_on_blob = state;
}

bool Blob::there_are_springs() {
    return springs_on_blob;
}

bool Blob::there_is_ssint() {
    return ssint_on_blob;
}

bool Blob::there_are_beads() {
    return beads_on_blob;
}

