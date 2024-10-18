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
    num_l_ctf = 0;
    num_r_ctf = 0;
    num_sltotal_ctf = 0;
    num_slsets_ctf = 0;
    num_surface_nodes = 0;
    num_interior_nodes = 0;
    num_surface_elements = 0;
    num_interior_elements = 0;
    ctf_l_nodes = {};
    ctf_r_nodes = {};
    ctf_sl_faces = {};
    ctf_slsurf_ndx = {};
    ctf_sl_faces = {};
    ctf_sl_surfsize = {};
    ctf_l_forces = {};
    ctf_r_forces = {};
    ctf_r_axis = {};
    ctf_r_type = {};
    ctf_sl_forces = {};
    mass = 0;
    blob_state = FFEA_BLOB_IS_STATIC;
    node = {};
    node_position = {};
    elem = {};
    surface = {};
    binding_site = {};
    solver = nullptr;
    linear_solver = 0;
    mass_in_blob = false;
    ssint_on_blob = false;
    springs_on_blob = false;
    beads_on_blob = false;
    force = {};
    poisson_solver.reset();
    phi_Omega = {};
    phi_Gamma = {};
    q = {};
    // nodal_q = {};
    poisson_rhs = {};
    pinned_nodes_list = {};
    bsite_pinned_nodes_list = {};

    num_contributing_faces = {};

    pbc_count[0] = 0;
    pbc_count[1] = 0;
    pbc_count[2] = 0;

    toBePrinted_nodes = {};

    poisson_surface_matrix = nullptr;
    poisson_interior_matrix = nullptr;
}

Blob::~Blob() {
    /* Release the node, element and surface arrays */
    node.clear();
    node_position.clear();
    elem.clear();
    surface.clear();
    
    binding_site.clear();

    /* Release the force vector */
    force.clear();

    /* Release the Solver */
    solver.reset();
    linear_solver = 0;
    mass_in_blob = false;
    ssint_on_blob = false;
    springs_on_blob = false;
    beads_on_blob = false;

    /* Release the Poisson Solver */
    poisson_solver.reset();
    phi_Omega.clear();
    phi_Gamma.clear();
    q.clear();
    // nodal_q.clear();
    poisson_rhs.clear();
    pinned_nodes_list.clear();

    /* delete precomp stuff */
    if (!bead_position.empty()) {
        bead_position.clear();
        bead_type.clear();
    }

    /* delete ctforces stuff */
    //  first linear:
    if (num_l_ctf > 0) {
        num_l_ctf = 0;
        ctf_l_nodes.clear();
        ctf_l_forces.clear();
    }
    // then rotational:
    if (num_r_ctf > 0) {
        num_r_ctf = 0;
        ctf_r_nodes.clear();
        ctf_r_forces.clear();
        ctf_r_axis.clear();
        ctf_r_type.clear();
    }
    // and then linear onto surfaces:
    if (num_sltotal_ctf > 0) {
        num_sltotal_ctf = 0;
        num_slsets_ctf = 0;
        ctf_slsurf_ndx.clear();
        ctf_sl_faces.clear();
        ctf_sl_surfsize.clear();
        ctf_sl_forces.clear();
    }

    /* Set relevant data to zero */
    conformation_index = 0;
    previous_conformation_index = 0;
    state_index = 0;
    previous_state_index = 0;
    num_surface_elements = 0;
    num_interior_elements = 0;
    num_surface_nodes = 0;
    num_interior_nodes = 0;
    blob_state = FFEA_BLOB_IS_STATIC;
    mass = 0;
    rng.reset();
    bsite_pinned_nodes_list.clear();

    toBePrinted_nodes.clear();
    
    num_contributing_faces.clear();
    
    poisson_surface_matrix.reset();
    poisson_interior_matrix.reset();
}


void Blob::config(const int blob_index, const int conformation_index, const string& node_filename,
                 const string& topology_filename, const string& surface_filename, const string& material_params_filename,
                 const string& stokes_filename, const string& ssint_filename, const string& pin_filename,
                 const string& binding_filename, const string& beads_filename, scalar scale, scalar calc_compress,
                 scalar compress, int linear_solver, int blob_state, const SimulationParams &params,
                 const PreComp_params &pc_params, SSINT_matrix *ssint_matrix,
                 BindingSite_matrix *binding_matrix, std::shared_ptr<std::vector<RngStream>> &rng){

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

    // Store the simulation parameters
    this->params = params;

    // Store the pointer to the random number generator array
    this->rng = rng;

    // Get the blob state
    this->blob_state = blob_state;

    // Need to know solver type
    this->linear_solver = linear_solver;
}

void Blob::init(){
    //Load the node, topology, surface, materials and stokes parameter files.
    load_nodes(s_node_filename.c_str(), scale);

    if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        load_topology(s_topology_filename.c_str());
        load_surface(s_surface_filename.c_str());
    } else {
        load_surface_no_topology(s_surface_filename.c_str());
    }

    if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
        load_material_params(s_material_params_filename.c_str());
        load_stokes_params(s_stokes_filename.c_str(), scale);
    }

    if (params.calc_ssint == 1 || params.calc_steric == 1) {
        load_ssint(s_ssint_filename.c_str(), ssint_matrix->get_num_types(), params.ssint_type);
    }

    if (params.calc_preComp == 1) {
        if (strcmp(s_beads_filename.c_str(), "") == 0) {
            FFEA_CAUTION_MESSG("No beads file was assigned to Blob %d\n", blob_index);
        }
        load_beads(s_beads_filename.c_str(), scale);
    }

    if (params.calc_ctforces == 1) {
        load_ctforces(params.ctforces_fname);
    }

    // Kinetic binding sites are still a structural property
    load_binding_sites();

    // This can be defaulted by leaving it out of the script
    if (blob_state == FFEA_BLOB_IS_DYNAMIC && s_pin_filename != "") {
        load_pinned_nodes(s_pin_filename.c_str());
    } else {
        pinned_nodes_list.clear();
    }

    // Linearise all the elements
    for (auto &elem_i : elem) {
        elem_i.linearise_element();
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
        if((params.calc_stokes == 0 && params.calc_noise == 0) && (params.calc_springs == 1 || params.calc_ctforces == 1)) {
            if(linear_solver == FFEA_NOMASS_CG_SOLVER) {
                FFEA_CAUTION_MESSG("For Blob %d, Conformation %d:\n\tUsing Springs / Constant forces in conjuction with the CG_nomass solver without calc_stokes or calc_noise may not converge, specially if there are very few ctforces or springs.\n\tIf you encounter this problem, consider setting either calc_stokes or calc_noise (you may try a low temperature) re-run the system :)\n", blob_index, conformation_index)
            }
        }

    // Create the chosen linear equation Solver for this Blob
	// Only initialise for dynamic blobs. Static don not need to be mechanically solved
	if (blob_state == FFEA_BLOB_IS_DYNAMIC) {
		if (linear_solver == FFEA_DIRECT_SOLVER) {
		    solver = std::make_unique<SparseSubstitutionSolver>();
		    mass_in_blob = true;
		} else if (linear_solver == FFEA_ITERATIVE_SOLVER) {
		    solver = std::make_unique<ConjugateGradientSolver>();
		    mass_in_blob = true;
		} else if (linear_solver == FFEA_MASSLUMPED_SOLVER) {
		    solver = std::make_unique<MassLumpedSolver>();
		    mass_in_blob = true;
		} else if (linear_solver == FFEA_NOMASS_CG_SOLVER) {
		    solver = std::make_unique<NoMassCGSolver>();
		} else {
		    throw FFEAException("Error in Blob initialisation: linear_solver=%d is not a valid solver choice\n", linear_solver);
		}
		if (!solver) throw FFEAException("No solver to work with");
	}

        // Initialise the Solver (whatever it may be)
        printf("\t\tBuilding solver:\n");
        solver->init(node, elem, params, pinned_nodes_list, bsite_pinned_nodes_list);
    }


    // Allocate the force vector array for the whole Blob
    try {
        force = std::vector<arr3>(node.size(), {0,0,0});
    } catch (std::bad_alloc &) {
        throw FFEAException("Could not store the force vector array");
    }

    // Calculate how many faces each surface node is a part of
    try {
        num_contributing_faces = std::vector<int>(num_surface_nodes);
    } catch(std::bad_alloc &) {
        throw FFEAException("Failed to allocate num_contributing_faces");
    }
    
    for (int i = 0; i < num_surface_nodes; i++) {
        num_contributing_faces[i] = 0;
    }
    for (int i = 0; i < surface.size(); i++) {
        for (int j = 0; j < 3; j++) {
            num_contributing_faces[surface[i].n[j]->index]++;
        }
    }

    // Store stokes drag on nodes, for use in viscosity matrix
    if (params.calc_stokes == 1) {
        for (auto &node_i : node) {
            node_i.stokes_drag = 6.0 * ffea_const::pi * params.stokes_visc * node_i.stokes_radius;
        }
    }

    // Generate the Mass matrix for this Blob
    if (params.calc_es == 1) {
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
        for (int el = 0; el < elem.size(); el++) {
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
        poisson_solver = std::make_unique<CG_solver>();
        poisson_solver->init(num_interior_nodes, params.epsilon2, params.max_iterations_cg);

        // Create the vector containing all values of the potential at each interior node
        try {
            phi_Omega = std::vector<scalar>(num_interior_nodes, 0);
        } catch (std::bad_alloc&) {
            throw FFEAException("Could not allocate memory (for phi_Omega array).");
        }

        // Create the vector containing all values of the potential on each surface node
        try {
            phi_Gamma = std::vector<scalar>(num_surface_nodes, 0);
        } catch (std::bad_alloc&) {
            throw FFEAException("Could not allocate memory (for phi_Gamma array).");
        }

        /*
                        // Create the vector containing the smoothed out charge distribution across the elements
                        try {
                            q = std::vector<scalar>(node.size(), 0);
                        } catch (std::bad_alloc&) {
                            throw FFEAException("Could not allocate memory (for q array).\n");
                        }
                        // Create the vector containing all charges on each node of the Blob
                        try {
                            nodal_q = std::vector<scalar>(node.size(), 0);
                        } catch (std::bad_alloc&) {
                            throw FFEAException("Could not allocate memory (for nodal_q array).\n");
                        }

                        // Set the charges on the Blob
                        for(int i = 0; i < node.size(); i++) {
        //			nodal_q[i] = (node[i].pos[2] - (2.9e-8/2.0)) * 1e35; //1e26;
                                nodal_q[i] = 1e26;
        //			printf("%e %e\n", node[i].pos[2], nodal_q[i]);
                        }

                        // Apply mass matrix to obtain smooth extrapolated charge density across elements
                        M->apply(nodal_q, q);

                        for(int i = 0; i < node.size(); i++) {
                                node[i].rho = q[i];
                        }
         */

        // Create the vector containing the charge distribution across the elements
        try {
            q = std::vector<scalar>(node.size(), 0);
        } catch (std::bad_alloc&) {
            throw FFEAException("Could not allocate memory (for q array).");
        }

        // const scalar charge_density = 6.0e25;
        const scalar charge_density = 6.0e25 * mesoDimensions::volume / mesoDimensions::charge;
        for (int n = 0; n < elem.size(); n++) {
            for (int i = 0; i < 10; i++) {
                q[elem[n].n[i]->index] += charge_density * elem[n].vol_0;
            }
        }

        for (int i = 0; i < node.size(); i++) {
            node[i].rho = q[i];
        }

        // for(int i = 0; i < node.size(); i++) {
        //     printf("q[%d] = %e\n", i, q[i]);
        // }

        // Create the Right Hand Side vector for the Poisson solver
        try {
            poisson_rhs = std::vector<scalar>(num_interior_nodes, 0);
        } catch (std::bad_alloc&) {
            throw FFEAException("Could not allocate memory (for poisson_rhs array).");
        }

        // build_poisson_matrices_and_setup_for_solve();
        // poisson_solver->solve(poisson_interior_matrix, phi_Omega, poisson_rhs);

        printf("\t\tdone.\n");
    }

    //If calc_compress option set in input script (default off), compress by factor specified in script
    if (calc_compress == 1) {
        compress_blob(compress);
    }

    if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
        toBePrinted_nodes = std::vector<scalar>(10 * node.size());
    } else {
        if (params.calc_es == 0) {
            toBePrinted_nodes = std::vector<scalar>(3 * node.size());
        } else {
            toBePrinted_nodes = std::vector<scalar>(4 * node.size());
        }
    }
}

bool Blob::check_inversion() {
	int n;
	matrix3 J;

	vector<int> invEls;
	invEls.clear();

    for (n = 0; n < elem.size(); n++) {

        // calculate jacobian for this element
        elem[n].calculate_jacobian(J);

        // get the 12 derivatives of the shape functions (by inverting the jacobian)
        // and also get the element volume. The function returns an error in the
        // case of an element inverting itself (determinant changing sign since last step)
        if (elem[n].calc_shape_function_derivatives_and_volume(J)) {
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
		return true;
	}

	return false;
}

void Blob::update_internal_forces() {
    if (blob_state != FFEA_BLOB_IS_DYNAMIC) {
        return;
    }

    if(pinned_nodes_list.size() == node.size()) {
        return;
    }

    if (params.calc_ctforces) apply_ctforces();

    /* some "work" variables */
    matrix3 J; // Holds the Jacobian calculated for the *current* element being processed
    matrix3 stress; // Holds the current stress tensor (elastic stress, with thermal fluctuations)
    vector12 du; // Holds the force change for the current element
    int tid; // Holds the current thread id (in parallel regions)
    int num_inversions = 0; // Counts the number of elements that have inverted (if > 0 then simulation has failed)

    // Element loop
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel default(none) private(J, stress, du, tid) reduction(+:num_inversions)
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
        for (int n = 0; n < elem.size(); n++) {

            // calculate jacobian for this element
            elem[n].calculate_jacobian(J);

            // get the 12 derivatives of the shape functions (by inverting the jacobian)
            // and also get the element volume. The function returns an error in the
            // case of an element inverting itself (determinant changing sign since last step)
            if (elem[n].calc_shape_function_derivatives_and_volume(J)) {
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

            if (params.calc_noise == 1) {
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

            if (params.calc_es == 1) {
                elem[n].calculate_electrostatic_forces();
            }
        }
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    }
#endif

    // Check if any elements have inverted
    if (num_inversions != 0) {
        if (num_inversions == 1) {
            throw FFEAException("1 element has inverted since the last step. Aborting simulation.\n");
        } else {
            throw FFEAException("%d elements have inverted since the last step. Aborting simulation.\n", num_inversions);
        }
    }
}

void Blob::update_positions() {
    // Aggregate forces on nodes from all elements (if not static)
    if (get_motion_state() != FFEA_BLOB_IS_DYNAMIC) {
	    return;
    }

    aggregate_forces_and_solve();

    // Update node velocities and positions
    euler_integrate();

    // todo: remove this
    //for (int i=0; i<this->get_num_elements(); i++){
    //    for (int j=0; j<9; j++){
    //        printf("%f", this->get_element(i)->n[j]->pos[0]);
    //    }
    //}

    // Linearise the 2nd order elements
    for (int n = 0; n < elem.size(); n++) {
        elem[n].linearise_element();
    }
}

void Blob::reset_solver() {
    // Delete and rebuild (to make sure everything is overwritten)
    solver->init(node, elem, params, pinned_nodes_list, bsite_pinned_nodes_list);
}

void Blob::translate_linear(arr3 *vec) {
    // Get a mapping from all node indices to just linear node indices
    int num_linear_nodes = get_num_linear_nodes();
    vector<int> map(num_linear_nodes);
    int j = 0;
    for(int i = 0; i < node.size(); ++i) {
        if(node[i].am_I_linear()) {
            map[j++] = i;
        }
    }

    // Translate linear nodes
    for(int i = 0; i < num_linear_nodes; ++i) {
        node[map[i]].pos[0] += vec[i][0];
        node[map[i]].pos[1] += vec[i][1];
        node[map[i]].pos[2] += vec[i][2];
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
    arr3 com;
    scalar x, y, z;

    get_centroid(com);

    // Move all nodes to the origin:
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) shared(com)
#endif
    for (int i = 0; i < node.size(); i++) {
        node[i].pos[0] -= com[0];
        node[i].pos[1] -= com[1];
        node[i].pos[2] -= com[2];
    }

    // Do the actual rotation and bring the nodes back to its initial position:
    for (int i = 0; i < node.size(); i++) {
        x = node[i].pos[0];
        y = node[i].pos[1];
        z = node[i].pos[2];

        node[i].pos[0] = x * r11 + y * r12 + z * r13 + com[0];
        node[i].pos[1] = x * r21 + y * r22 + z * r23 + com[1];
        node[i].pos[2] = x * r31 + y * r32 + z * r33 + com[2];


    }


    if (beads) {
        // Move all beads to the origin:
        for (int i = 0; i < bead_position.size(); ++i) {
            bead_position[i][0] -= com[0];
            bead_position[i][1] -= com[1];
            bead_position[i][2] -= com[2];
        }
        // Do the actual rotation and bring the beads back to its initial position:
        for (int i = 0; i < bead_position.size(); ++i) {
            x = bead_position[i][0];
            y = bead_position[i][1];
            z = bead_position[i][2];

            bead_position[i][0] = x * r11 + y * r12 + z * r13 + com[0];
            bead_position[i][1] = x * r21 + y * r22 + z * r23 + com[1];
            bead_position[i][2] = x * r31 + y * r32 + z * r33 + com[2];
        }
    }
}


arr3 Blob::position(scalar x, scalar y, scalar z) {
    scalar centroid_x = 0.0, centroid_y = 0.0, centroid_z = 0.0;
    arr3 v;
    // scalar dx, dy, dz;

    // Calculate centroid of entire Blob mesh
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(shared) reduction(+:centroid_x,centroid_y,centroid_z)
#endif
    for (int i = 0; i < node.size(); ++i) {
        centroid_x += node[i].pos[0];
        centroid_y += node[i].pos[1];
        centroid_z += node[i].pos[2];
    }

    centroid_x *= (1.0 / node.size());
    centroid_y *= (1.0 / node.size());
    centroid_z *= (1.0 / node.size());

    // Calculate displacement vector required to move centroid to requested position
    v[0] = x - centroid_x;
    v[1] = y - centroid_y;
    v[2] = z - centroid_z;

    // Move all nodes in mesh by displacement vector
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) shared(v)
#endif
    for (int i = 0; i < node.size(); ++i) {
        node[i].pos[0] += v[0];
        node[i].pos[1] += v[1];
        node[i].pos[2] += v[2];
    }

    return v;

}

void Blob::position_beads(scalar x, scalar y, scalar z) {
    for (int i = 0; i < bead_position.size(); i ++) {
        bead_position[i][0] += x;
        bead_position[i][1] += y;
        bead_position[i][2] += z;
    }
}

void Blob::move(scalar dx, scalar dy, scalar dz) {
    for (auto& node_i : node) {
        node_i.pos[0] += dx;
        node_i.pos[1] += dy;
        node_i.pos[2] += dz;
    }
    for(int i = 0; i < get_num_faces(); ++i) {
        if(get_motion_state() != FFEA_BLOB_IS_DYNAMIC && surface[i].n[3] != nullptr) {
            surface[i].n[3]->pos[0] += dx;
            surface[i].n[3]->pos[1] += dy;
            surface[i].n[3]->pos[2] += dz;
            //fprintf(stderr, "Surface %d = %f %f %f\n", i, dx, dy, dz);
        }
    }
}

void Blob::get_CoM(arr3 &com) {
    com[0] = 0;
    com[1] = 0;
    com[2] = 0;

    for (int n = 0; n < elem.size(); n++) {
        com[0] += (elem[n].mass * elem[n].n[0]->pos[0] + elem[n].n[1]->pos[0] + elem[n].n[2]->pos[0] + elem[n].n[3]->pos[0]) / 4;
        com[1] += (elem[n].mass * elem[n].n[0]->pos[1] + elem[n].n[1]->pos[1] + elem[n].n[2]->pos[1] + elem[n].n[3]->pos[1]) / 4;
        com[2] += (elem[n].mass * elem[n].n[0]->pos[2] + elem[n].n[1]->pos[2] + elem[n].n[2]->pos[2] + elem[n].n[3]->pos[2]) / 4;
    }
    if (elem.empty()) {
        return;
    } else {
        com[0] /= elem.size();
        com[1] /= elem.size();
        com[2] /= elem.size();
    }
}

void Blob::get_centroid(arr3 &com) {
    com[0] = 0;
    com[1] = 0;
    com[2] = 0;
    for (auto& node_i : node) {
        com[0] += node_i.pos[0];
        com[1] += node_i.pos[1];
        com[2] += node_i.pos[2];
    }
    com[0] /= node.size();
    com[1] /= node.size();
    com[2] /= node.size();
}

void Blob::calc_and_store_centroid(arr3 &com) {
    get_centroid(com);

    CoG[0] = com[0];
    CoG[1] = com[1];
    CoG[2] = com[2];
}

// This one returns an array rather than arsing about with pointers
arr3 Blob::calc_centroid() {
    arr3 com;
    com[0] = 0.0;
    com[1] = 0.0;
    com[2] = 0.0;
    for (const auto &node_i:node) {
        com[0] += node_i.pos[0];
        com[1] += node_i.pos[1];
        com[2] += node_i.pos[2];
    }
    com[0] /= node.size();
    com[1] /= node.size();
    com[2] /= node.size();
    return com;
}

std::vector<arr3 *> &Blob::get_actual_node_positions() {
    return node_position;
}

void Blob::copy_node_positions(std::vector<arr3*> &nodes) {
    for(size_t i = 0; i < node.size(); ++i) {
        (*nodes[i])[0] = node[i].pos[0];
        (*nodes[i])[1] = node[i].pos[1];
        (*nodes[i])[2] = node[i].pos[2];
    }
}

void Blob::set_node_positions(const std::vector<arr3*> &node_pos) {
    for (size_t i = 0; i < node.size(); ++i) {
        node[i].pos[0] = (*node_pos[i])[0];
        node[i].pos[1] = (*node_pos[i])[1];
        node[i].pos[2] = (*node_pos[i])[2];
    }
}

void Blob::set_pos_0() {
    CoG_0[0] = 0.0;
    CoG_0[1] = 0.0;
    CoG_0[2] = 0.0;
    for (auto& node_i : node) {
        node_i.pos_0[0] = node_i.pos[0];
        node_i.pos_0[1] = node_i.pos[1];
        node_i.pos_0[2] = node_i.pos[2];
        CoG_0[0] += node_i.pos_0[0];
        CoG_0[1] += node_i.pos_0[1];
        CoG_0[2] += node_i.pos_0[2];
    }
    CoG_0[0] /= node.size();
    CoG_0[1] /= node.size();
    CoG_0[2] /= node.size();
}

void Blob::kinetically_set_faces(bool state) {
    for(int i = 0; i < surface.size(); ++i) {
        surface[i].set_kinetic_state(state);
    }
}

void Blob::linearise_elements() {
    for(auto &elem_i : elem) {
        elem_i.linearise_element();
    }
}

void Blob::linearise_force() {
    int nIdx[NUM_NODES_QUADRATIC_TET];
    for (int i=0; i< elem.size(); i++) {
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
    arr3 cog,cogaft;
    this->get_centroid(cog);
    //loop moves nodes in by
    for (auto& node_i : node) {
        node_i.pos[0] *= compress;
        node_i.pos[1] *= compress;
        node_i.pos[2] *= compress;
    }
    this->get_centroid(cogaft);

    if (!(fabs(cogaft[0] - cog[0]*compress) < 0.000001&&fabs(cogaft[1] - cog[1]*compress) < 0.000001&&fabs(cogaft[2] - cog[2]*compress) < 0.000001)) {
        printf("FRIENDLY WARNING: Centre of Geometry has moved during compression. Compression feature implemented for spherical objects so other shapes may experience issues.\n");
    }
}

void Blob::create_viewer_node_file(const char *node_filename, scalar scale) {
    FILE *out;
    char new_node_filename[] = "VIEWERNODE_";
    int i;

    // Name new node file
    strcat(new_node_filename, node_filename);


    //open the new node file
    if (!(out = fopen(new_node_filename, "w"))) {
        throw FFEAFileException(new_node_filename);
    }
    printf("\t\tWriting to viewer nodes file: %s\n", new_node_filename);

    fprintf(out, "ffea viewer node file\n");
    fprintf(out, "num_nodes %zu\n", node.size());
    fprintf(out, "num_surface_nodes %d\n", num_surface_nodes);
    fprintf(out, "num_interior_nodes %d\n", num_interior_nodes);

    // Write all the nodes to file
    fprintf(out, "surface nodes:\n");
    for (i = 0; i < num_surface_nodes; i++) {
        fprintf(out, "%le %le %le\n", node[i].pos[0] / scale, node[i].pos[1] / scale, node[i].pos[2] / scale);
    }

    fprintf(out, "interior nodes:\n");
    for (auto& node_i : node) {
        fprintf(out, "%le %le %le\n", node_i.pos[0] / scale, node_i.pos[1] / scale, node_i.pos[2] / scale);
    }

    fclose(out);
    printf("\t\t\tWrote %d nodes from %s\n", i, new_node_filename);
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
        for (size_t i = 0; i < node.size(); i++) {
            fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                    node[i].pos[0]*mesoDimensions::length, node[i].pos[1]*mesoDimensions::length, node[i].pos[2]*mesoDimensions::length,
                    node[i].vel[0]*mesoDimensions::velocity, node[i].vel[1]*mesoDimensions::velocity, node[i].vel[2]*mesoDimensions::velocity,
                    node[i].phi,
                    force[i][0]*mesoDimensions::force, force[i][1]*mesoDimensions::force, force[i][2]*mesoDimensions::force);
        }
    } else {
        if (params.calc_es == 0) {
            for (auto& node_i : node) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                    node_i.pos[0]*mesoDimensions::length, node_i.pos[1]*mesoDimensions::length, node_i.pos[2]*mesoDimensions::length,
                        0., 0., 0., 0., 0., 0., 0.);
            }
        } else {
            for (auto& node_i : node) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                    node_i.pos[0]*mesoDimensions::length, node_i.pos[1]*mesoDimensions::length, node_i.pos[2]*mesoDimensions::length,
                    0., 0., 0.,
                    node_i.phi, 0., 0., 0.);
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
        for (size_t i = 0; i < node.size(); i++) {
            fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                    toBePrinted_nodes[10*i], toBePrinted_nodes[10*i+1], toBePrinted_nodes[10*i+2],
                    toBePrinted_nodes[10*i+3], toBePrinted_nodes[10*i+4], toBePrinted_nodes[10*i+5],
                    toBePrinted_nodes[10*i+6], toBePrinted_nodes[10*i+7], toBePrinted_nodes[10*i+8],
                    toBePrinted_nodes[10*i+9]);
        }
    } else {
        if (params.calc_es == 0) {
            for (size_t i = 0; i < node.size(); i++) {
                fprintf(trajectory_out, "%e %e %e %e %e %e %e %e %e %e\n",
                        toBePrinted_nodes[3*i], toBePrinted_nodes[3*i+1], toBePrinted_nodes[3*i+2],
                        0., 0., 0., 0., 0., 0., 0.);
            }
        } else {
            for (size_t i = 0; i < node.size(); i++) {
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
        for (int i = 0; i < node.size(); i++) {
            toBePrinted_nodes[10*i   ] = node[i].pos[0]*mesoDimensions::length;
            toBePrinted_nodes[10*i +1] = node[i].pos[1]*mesoDimensions::length;
            toBePrinted_nodes[10*i +2] = node[i].pos[2]*mesoDimensions::length;
            toBePrinted_nodes[10*i +3] = node[i].vel[0]*mesoDimensions::velocity;
            toBePrinted_nodes[10*i +4] = node[i].vel[1]*mesoDimensions::velocity;
            toBePrinted_nodes[10*i +5] = node[i].vel[2]*mesoDimensions::velocity;
            toBePrinted_nodes[10*i +6] = node[i].phi;
            toBePrinted_nodes[10*i +7] = force[i][0]*mesoDimensions::force;
            toBePrinted_nodes[10*i +8] = force[i][1]*mesoDimensions::force;
            toBePrinted_nodes[10*i +9] = force[i][2]*mesoDimensions::force;
        }
    } else {
        if (params.calc_es == 0) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            #pragma omp parallel for default(none)
#endif
            for (int i = 0; i < node.size(); i++) {
                toBePrinted_nodes[3*i   ] = node[i].pos[0]*mesoDimensions::length;
                toBePrinted_nodes[3*i +1] = node[i].pos[1]*mesoDimensions::length;
                toBePrinted_nodes[3*i +2] = node[i].pos[2]*mesoDimensions::length;
            }
        } else {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            #pragma omp parallel for default(none)
#endif
            for (int i = 0; i < node.size(); i++) {
                toBePrinted_nodes[4*i   ] = node[i].pos[0]*mesoDimensions::length;
                toBePrinted_nodes[4*i +1] = node[i].pos[1]*mesoDimensions::length;
                toBePrinted_nodes[4*i +2] = node[i].pos[2]*mesoDimensions::length;
                toBePrinted_nodes[4*i +3] = node[i].phi;
            }
        }
    }
}

void Blob::read_nodes_from_file(FILE *trajectory_out) {
    char state_str[20];
    char *result = nullptr;

    // If blob is static, don't read any nodes. Simply read the word "STATIC"
    if (blob_state == FFEA_BLOB_IS_STATIC) {
        result = fgets(state_str, 20, trajectory_out);
        if (!result) {
            throw FFEAException("Problem when reading 'STATIC' (expected) line in trajectory file.");
        }
        if (strcmp(state_str, "STATIC\n") != 0) {
            throw FFEAException("When restarting from trajectory file, expected to read 'STATIC', but instead found '%s...'.", state_str);
        }
        return;
    } else {
        result = fgets(state_str, 20, trajectory_out);
        if (!result) {
            throw FFEAException("Problem when reading state line in trajectory file.");
        }
    }

    for (int i = 0; i < node.size(); i++) {
        if (fscanf(trajectory_out, "%le %le %le %le %le %le %le %le %le %le\n", &node[i].pos[0], &node[i].pos[1], &node[i].pos[2], &node[i].vel[0], &node[i].vel[1], &node[i].vel[2], &node[i].phi, &force[i][0], &force[i][1], &force[i][2]) != 10) {
            throw FFEAException("(When restarting) Error reading from trajectory file, for node %d", i);
        } else {
            node[i].pos[0] /= mesoDimensions::length;
            node[i].pos[1] /= mesoDimensions::length;
            node[i].pos[2] /= mesoDimensions::length;
            node[i].vel[0] /= mesoDimensions::velocity;
            node[i].vel[1] /= mesoDimensions::velocity;
            node[i].vel[2] /= mesoDimensions::velocity;
            force[i][0] /= mesoDimensions::force;
            force[i][1] /= mesoDimensions::force;
            force[i][2] /= mesoDimensions::force;
        }

    }
}

void Blob::calculate_deformation() {
    int num_inversions = 0;
    matrix3 J;
    for (int n = 0; n < elem.size(); n++) {

        // calculate jacobian for this element
        elem[n].calculate_jacobian(J);

        // get the 12 derivatives of the shape functions (by inverting the jacobian)
        // and also get the element volume. The function returns an error in the
        // case of an element inverting itself (determinant changing sign since last step)
        if (elem[n].calc_shape_function_derivatives_and_volume(J)) {
            FFEA_error_text();
            printf("Element %d has inverted during deformation calculation\n", n);
            num_inversions++;
        }

        // And F_ij
        elem[n].calc_deformation(J);
    }

    if(num_inversions > 0) {
        throw FFEAException("Deformation calculation found inversions!");  // Note, the previous error return was not checked
    }
}

scalar Blob::calc_volume() {
    scalar volume = 0.0;
    for(int i = 0; i < elem.size(); ++i) {
        volume += elem[i].calc_volume();
    }
    return volume;
}

scalar Blob::calculate_strain_energy() {
    int n;
    scalar C, detF, strain_energy = 0.0;
    calculate_deformation();
    for(n = 0; n < elem.size(); ++n) {
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
    arr3 total_ssint_xz_force;
    scalar total_ssint_xz_energy = 0.0, total_ssint_xz_area = 0.0;
    scalar temp1, temp2, temp3;
    scalar cogx = 0.0, cogy = 0.0, cogz = 0.0, lx = 0.0, ly = 0.0, lz = 0.0, brmsd = 0.0;
    scalar r[4][3];
    vector12 vec;

    kenergy = 0.0;
    senergy = 0.0;

    // OpenMP can't reduce members of classes :(
    //arr3_set_zero(L);
    //arr3_set_zero(CoG);
    arr3_set_zero(CoM);

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) reduction(+:kenergy, senergy) private(n, vec, temp1, temp2, temp3)
#endif
    for (n = 0; n < elem.size(); n++) {
        if (linear_solver != FFEA_NOMASS_CG_SOLVER) {
            /*
             * Kinetic energy contribution:
             */

            // Read the u vector for this element
            elem[n].get_element_velocity_vector(vec);

            // Apply the mass matrix
            elem[n].apply_element_mass_matrix(vec);

            // Dot u with M.u to get the contribution to the kinetic energy
            kenergy += elem[n].n[0]->vel[0] * vec[0] +
                       elem[n].n[1]->vel[0] * vec[1] +
                       elem[n].n[2]->vel[0] * vec[2] +
                       elem[n].n[3]->vel[0] * vec[3] +
                       elem[n].n[0]->vel[1] * vec[4] +
                       elem[n].n[1]->vel[1] * vec[5] +
                       elem[n].n[2]->vel[1] * vec[6] +
                       elem[n].n[3]->vel[1] * vec[7] +
                       elem[n].n[0]->vel[2] * vec[8] +
                       elem[n].n[1]->vel[2] * vec[9] +
                       elem[n].n[2]->vel[2] * vec[10] +
                       elem[n].n[3]->vel[2] * vec[11];

            /*
             * Centre of Mass contribution:
             */
            //com_x += elem[n].mass * (elem[n].n[0]->pos[0] + elem[n].n[1]->pos[0] + elem[n].n[2]->pos[0] + elem[n].n[3]->pos[0]) / 4;
            //com_y += elem[n].mass * (elem[n].n[0]->pos[1] + elem[n].n[1]->pos[1] + elem[n].n[2]->pos[1] + elem[n].n[3]->pos[1]) / 4;
            //com_z += elem[n].mass * (elem[n].n[0]->pos[2] + elem[n].n[1]->pos[2] + elem[n].n[2]->pos[2] + elem[n].n[3]->pos[2]) / 4;
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
    get_CoM(CoM);
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
            std::array{.1, .05, .05, .05},
            std::array{.05, .1, .05, .05},
            std::array{.05, .05, .1, .05},
            std::array{.05, .05, .05, .1}
        };

#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none) reduction(+:lx, ly, lz) shared(MM) private(n, r, temp1, temp2, temp3, i, j)
#endif
        for (n = 0; n < elem.size(); n++) {
            // Find the separation vectors for this element
            for (i = 0; i < 4; i++) {
                r[i][0] = elem[n].n[i]->pos[0] - CoM[0];
                r[i][1] = elem[n].n[i]->pos[1] - CoM[1];
                r[i][2] = elem[n].n[i]->pos[2] - CoM[2];
            }

            // Calculate contribution to angular momentum from this element
            temp1 = 0;
            temp2 = 0;
            temp3 = 0;
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    temp1 += MM[i][j] * (r[j][1] * elem[n].n[i]->vel[2] - r[j][2] * elem[n].n[i]->vel[1]);
                    temp2 += MM[i][j] * (r[j][2] * elem[n].n[i]->vel[0] - r[j][0] * elem[n].n[i]->vel[2]);
                    temp3 += MM[i][j] * (r[j][0] * elem[n].n[i]->vel[1] - r[j][1] * elem[n].n[i]->vel[0]);
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
    for (const auto &node_i : node) {
        /*
         * Center of geometry contribution (geometry better on nodes than elements)
         */
        cogx += node_i.pos[0];
        cogy += node_i.pos[1];
        cogz += node_i.pos[2];

    }
    CoG[0] = cogx / node.size();
    CoG[1] = cogy / node.size();
    CoG[2] = cogz / node.size();

    // Remove translation and rotation (maybe make optional in future)
    bool remTrans = false;
    // bool remRot = true;

    if(remTrans) {
        for(const auto& node_i : node) {
            temp1 = node_i.pos[0] - node_i.pos_0[0] + CoG_0[0] - CoG[0];
            temp2 = node_i.pos[1] - node_i.pos_0[1] + CoG_0[1] - CoG[1];
            temp3 = node_i.pos[2] - node_i.pos_0[2] + CoG_0[2] - CoG[2];
            brmsd += temp1 * temp1 + temp2 * temp2 + temp3*temp3;
        }
    } else {
        for(const auto& node_i : node) {
            temp1 = node_i.pos[0] - node_i.pos_0[0];
            temp2 = node_i.pos[1] - node_i.pos_0[1];
            temp3 = node_i.pos[2] - node_i.pos_0[2];
            brmsd += temp1 * temp1 + temp2 * temp2 + temp3*temp3;
        }
    }
    rmsd = sqrt(brmsd / node.size());
    L[0] = lx;
    L[1] = ly;
    L[2] = lz;
}

void Blob::write_measurements_to_file(FILE *fout) {
    // White space for blob index bit
    fprintf(fout, "     ");
    if(there_is_mass()) {
        fprintf(fout, "%-14.6e", kineticenergy * mesoDimensions::Energy);
    }
    fprintf(fout, "%-14.6e", strainenergy * mesoDimensions::Energy);
    fprintf(fout, "%-14.6e%-14.6e%-14.6e%-14.6e", CoG[0] * mesoDimensions::length, CoG[1] * mesoDimensions::length, CoG[2] * mesoDimensions::length, rmsd * mesoDimensions::length);
    fflush(fout);
}

void Blob::make_stress_measurements(FILE *stress_out, int blob_number) {
    int n;
    if (stress_out != nullptr) {
        fprintf(stress_out, "blob\t%d\n", blob_number);
        for (n = 0; n < elem.size(); n++) {
            fprintf(stress_out, "%e\n", elem[n].internal_stress_mag);
        }
        fprintf(stress_out, "\n");
    }
}

void Blob::calc_centroids_and_normals_of_all_faces() {
    int i;
    for (i = 0; i < surface.size(); i++)
        surface[i].calc_area_normal_centroid();
}

void Blob::calc_all_centroids() {
    int i;
    for (i = 0; i < surface.size(); i++) {
        surface[i].calc_area_normal_centroid();
    }
    for(i = 0; i < elem.size(); ++i) {
        elem[i].calc_centroid();
    }
}

int Blob::get_num_faces() {
    return surface.size();
}

int Blob::get_num_beads() {
    return bead_position.size();
}

int Blob::getNumBindingSites() {
    return binding_site.size();
}

bool Blob::is_using_beads() {
    return !bead_position.empty();
}

scalar Blob::get_rmsd() {
    return rmsd;
}

int Blob::get_linear_solver() {
    return linear_solver;
};

/*
 *
arr3 Blob::get_CoG() {
    return CoG;
}
 */
void Blob::get_stored_centroid(arr3 &cog){
    store(CoG, cog);
}

Face * Blob::get_face(int i) {
    if (surface[i].is_ssint_active() == true) {
        return &surface[i];
    } else {
        return nullptr;
    }
}

Face * Blob::absolutely_get_face(int i) {
    return &surface[i];
}

tetra_element_linear *Blob::get_element(int i) {
    return &elem[i];
}

/**
 * @brief returns the position of bead i.
 *
 * @ingroup FMM
 **/
void Blob::get_bead_position(int i, arr3 &v) {
    store(bead_position[i], v);
}

/**
 * @brief returns the bead_type pointer.
 *
 * @ingroup FMM
 **/
std::vector<int> &Blob::get_bead_types() {
    return bead_type;

}

/**
 * @brief returns the bead_type of bead i.
 *
 * @ingroup FMM
 **/
int Blob::get_bead_type(int i) const {
    return bead_type[i];
}

/**
 * @brief returns the list of nodes where bead i should be assigned to.
 *
 * @ingroup FMM
 **/
std::vector<int> &Blob::get_bead_assignment(int i) {
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

void Blob::build_linear_node_elasticity_matrix(Eigen::SparseMatrix<scalar> *A) {
    // matrix3 J, stress;
    vector12 elastic_force[2];
    vector<Eigen::Triplet<scalar> > components;

    // Firstly, get a mapping from all node indices to just linear node indices
    // int num_linear_nodes = get_num_linear_nodes();
    vector<int> map(node.size());
    {
        int j = 0;
        for (size_t i = 0; i < node.size(); ++i) {
            if (node[i].am_I_linear()) {
                map[i] = j;
                j++;
            } else {
                map[i] = -1;
            }
        }
    }

    // For each element
    for(int elem_index = 0; elem_index < elem.size(); ++elem_index) {
        // Calculate dx, how far each node should be moved for a linearisataion, as well as unstrained parameter
        scalar dx = cbrt(elem[elem_index].calc_volume()) / 1000.0;

        // For every node a in every direction i
        for(int a = 0; a < 4; ++a) {

            // Get global index for node a
            int global_a = elem[elem_index].n[a]->index;
            int global_a_lin = map[elem[elem_index].n[a]->index];

            for(int i = 0; i < 3; ++i) {
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
                for(int b = 0; b < 4; ++b) {
                    // Get global index for node b
                    // int global_b = elem[elem_index].n[b]->index;
                    int global_b_lin = map[elem[elem_index].n[b]->index];

                    for(int j = 0; j < 3; ++j) {
                        scalar val = (1.0 / (2 * dx)) * (elastic_force[1][4 * j + b] - elastic_force[0][4 * j + b]);

                        // Row is dE_p, column dx_q. Not that it should matter! Directions then nodes i.e. x0,y0,z0,x1,y1,z1..[0]n,yn,zn
                        //row = num_linear_nodes * i + global_a;
                        //column = num_linear_nodes * j + global_b;
                        int row = 3 * global_a_lin + i;
                        int column = 3 * global_b_lin + j;
                        components.push_back(Eigen::Triplet<scalar>(row, column, val));
                    }
                }
            }
        }
    }

    // Now build the matrix
    A->setFromTriplets(components.begin(), components.end());
}

void Blob::build_linear_node_viscosity_matrix(Eigen::SparseMatrix<scalar> *K) {
    vector<Eigen::Triplet<scalar> > components;

    // Firstly, get a mapping from all node indices to just linear node indices
    vector<int> map(node.size());
    int offset = 0;
    for(int i = 0; i < node.size(); ++i) {
        if(node[i].am_I_linear()) {
            map[i] = i - offset;
        } else {
            offset += 1;
            map[i] = -1;
        }
    }

    //int num_linear_nodes = get_num_linear_nodes();

    // For each element
    for(int elem_index = 0; elem_index < elem.size(); ++elem_index) {
        // Calculate a local viscosity matrix
        matrix3 J;
        elem[elem_index].calculate_jacobian(J);
        elem[elem_index].calc_shape_function_derivatives_and_volume(J);
        elem[elem_index].create_viscosity_matrix();

        // Add each component of the local matrix to the global matrix

        //For each node
        for(int a = 0; a < 4; ++a) {
            int global_a = map[elem[elem_index].n[a]->index];
            for(int b = 0; b < 4; ++b) {
                int global_b = map[elem[elem_index].n[b]->index];

                // And each direction
                for(int i = 0; i < 3; ++i) {
                    for(int j = 0; j < 3; ++j) {
                        scalar val = elem[elem_index].viscosity_matrix[4 * i + a][4 * j + b];
                        //row = num_linear_nodes * i + global_a;
                        //column = num_linear_nodes * j + global_b;
                        int row = 3 * global_a + i;
                        int column = 3 * global_b + j;
                        components.push_back(Eigen::Triplet<scalar>(row, column, val));
                    }
                }
            }
        }
    }

    // For each node, add stokes if necessary
    if(params.calc_stokes == 1) {
        for(int global_a = 0; global_a < node.size(); ++global_a) {
            if(map[global_a] == -1) {
                continue;
            } else {
                for(int i = 0; i < 3; ++i) {
                    scalar val = node[map[global_a]].stokes_drag;
                    //row = 4 * i + map[global_a];
                    int row = 3 * map[global_a] + i;
                    int column = row;
                    components.push_back(Eigen::Triplet<scalar>(row, column, val));
                }
            }
        }
    }

    // Now build the matrix and get the symmetric part (just in case)
    K->setFromTriplets(components.begin(), components.end());
    components.clear();
    for(int j = 0; j < K->outerSize(); ++j) {
        for(Eigen::SparseMatrix<scalar>::InnerIterator it(*K,j); it; ++it) {
            components.push_back(Eigen::Triplet<scalar>(it.row(), it.col(), 0.5 * it.value()));
            components.push_back(Eigen::Triplet<scalar>(it.col(), it.row(), 0.5 * it.value()));
        }
    }
    K->setFromTriplets(components.begin(), components.end());
}

void Blob::build_linear_node_mass_matrix(Eigen::SparseMatrix<scalar> *M) {
    //MassMatrixQuadratic M_el;
    MassMatrixLinear M_el;
    vector<Eigen::Triplet<scalar> > components;

    // Firstly, get a mapping from all node indices to just linear node indices

    vector<int> map(node.size());
    int offset = 0;
    for(int i = 0; i < node.size(); ++i) {

        if(node[i].am_I_linear()) {
            //if(true) {
            map[i] = i - offset;
        } else {
            offset += 1;
            map[i] = -1;
        }
    }

    // num_linear_nodes = get_num_linear_nodes();

    // For each element
    scalar sum = 0.0, sum2 = 0.0;
    for(int elem_index = 0; elem_index < elem.size(); ++elem_index) {

        // Build mass matrix
        elem[elem_index].construct_element_mass_matrix(M_el);

        // Add each component of the local matrix to the global matrix

        // For each node
        for(int a = 0; a < 4; ++a) {
            int global_a = map[elem[elem_index].n[a]->index];

            for(int b = 0; b < 4; ++b) {
                int global_b = map[elem[elem_index].n[b]->index];

                // From another function that get's actual memory location. Dereference for value
                // Val same for each direction
                scalar val = M_el.get_M_alpha_value(a, b);
                //sum2 += val;
                //cout << a << " " << b << " " << val * mesoDimensions::mass << endl;
                // And each direction
                for(int i = 0; i < 3; ++i) {
                    int row = 3 * global_a + i;
                    int column = 3 * global_b + i;
                    components.push_back(Eigen::Triplet<scalar>(row, column, val));
                }
            }
        }
    }

    // Now build the matrix and get the symmetric part (just in case)
    M->setFromTriplets(components.begin(), components.end());
    components.clear();
    for(int j = 0; j < M->outerSize(); ++j) {
        for(Eigen::SparseMatrix<scalar>::InnerIterator it(*M,j); it; ++it) {
            components.push_back(Eigen::Triplet<scalar>(it.row(), it.col(), 0.5 * it.value()));
            components.push_back(Eigen::Triplet<scalar>(it.col(), it.row(), 0.5 * it.value()));
            //sum += it.value();
        }
    }
    //cout << "Blob mass matrix sums: " << mesoDimensions::mass * sum2 << " " << mesoDimensions::mass * sum << endl;
    M->setFromTriplets(components.begin(), components.end());
}

void Blob::build_linear_node_rp_diffusion_matrix(Eigen_MatrixX *D) {

    int i, j, num_linear_nodes;
    scalar mod, mod2, a;
    Eigen_Matrix3 block, rr;
    Eigen_Vector3 sep;

    // Firstly, get a mapping from all node indices to just linear node indices
    vector<int> map(node.size());
    int offset = 0;
    for(size_t i = 0; i < node.size(); ++i) {
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
    for(size_t n = 0; n < node.size(); ++n) {

        // If secondary node, continue
        if(map[n] == -1) {
            continue;
        } else {
            i = map[n];
        }
        for(size_t m = n; m < node.size(); ++m) {

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
                block = Eigen::Matrix3d::Identity() *  params.kT / node[n].stokes_drag;
            } else {

                // Get distance between nodes and the outer product of the separation
                Eigen::Vector3d vecn(node[n].pos[0], node[n].pos[1], node[n].pos[2]);
                Eigen::Vector3d vecm(node[m].pos[0], node[m].pos[1], node[m].pos[2]);

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
                    block *= params.kT / (8 * ffea_const::pi * params.stokes_visc * mod2 * mod);

                } else {
                    block = Eigen::Matrix3d::Identity() * (1 - ((9 * mod) / (32.0 * a)));
                    block += rr * (3 / (32.0 * a * mod));
                    block *= params.kT / (6 * ffea_const::pi * params.stokes_visc * a);
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
}

void Blob::solve_poisson(scalar *phi_gamma_IN, scalar *J_Gamma_OUT) {
    if (num_interior_nodes > 0) {
        /* Convert the given potential on the surface faces into the potential on each node */
        for (int i = 0; i < num_surface_nodes; i++) {
            phi_Gamma[i] = 0;
        }

        for (int i = 0; i < surface.size(); i++) {
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

        for (int i = 0; i < surface.size(); i++) {
            J_Gamma_OUT[i] = surface[i].get_normal_flux();
        }
    }
}


void Blob::apply_ctforces() {
    // CTF 1 - Add the linear constant forces, stored in ctf_force.
    // If there are no forces, num_l_ctf will be zero.
    for (int i=0; i<num_l_ctf; i++) {
        force[ctf_l_nodes[i]][0] += ctf_l_forces[3*i  ];
        force[ctf_l_nodes[i]][1] += ctf_l_forces[3*i+1];
        force[ctf_l_nodes[i]][2] += ctf_l_forces[3*i+2];
    }
    // CTF 2 - Add the rotational forces:
    for (int i=0; i<num_r_ctf; i++) {
        // 2.1 - get the axis:
        arr3 r, f;
        scalar rsize;
        // 2.1.1 - if it is defined by points:
        if (ctf_r_type[2*i] == 'p') {
            // get r, and its length rsize
            intersectingPointToLine(node[ctf_r_nodes[i]].pos,// point out of the line
                                         arr_view<scalar, 3>(ctf_r_axis, 6 * i), //  point
                                    arr_view<scalar, 3>(ctf_r_axis, (6*i+3)), r); // vector, axis
            sub(r, node[ctf_r_nodes[i]].pos, r);
        } /* else if (ctf_l_type[2*i] == 'n') {
        // not working yet!
        that should read:
        axis[0] = blob[ctf_r_axis[6*i + 3]].conf[ctf_r_axis[6*i + 4]].node[ctf_r_axis[6*i + 5]]->pos[0] -
                  blob[ctf_r_axis[6*i + 0]].conf[ctf_r_axis[6*i + 1]].node[ctf_r_axis[6*i + 2]]->pos[0]
        axis[1] =  ...
        axis[2] =  ...
        arr3Normalize(axis);
        // We do not have this information in "Blob"!!
      } */
        arr3arr3VectorProduct(r, arr_view<scalar, 3>(ctf_r_axis, +(6*i+3)), f);
        // 2.2.a - apply a constant circular force:
        if ( ctf_r_type[2*i + 1] == 'f' ) {
            resize(ctf_r_forces[i], f);
            // 2.2.b - apply a constant torque:
        } else if (ctf_r_type[2*i + 1] == 't') {
            rsize = magnitude(r);
            resize(ctf_r_forces[i]/rsize, f);
        }
        force[ctf_r_nodes[i]][0] += f[0];
        force[ctf_r_nodes[i]][1] += f[1];
        force[ctf_r_nodes[i]][2] += f[2];
    }

    // CTF 3 - Add the linear surface forces:
    scalar totalArea;
    arr3 traction;
    int auxndx = 0;
    std::vector<scalar> faceAreas = std::vector<scalar>();
    for (int i=0; i<num_slsets_ctf; i++) {
        totalArea = 0;
        faceAreas.reserve(ctf_sl_surfsize[i]);
        for (int j=0; j<ctf_sl_surfsize[i]; j++) {
            faceAreas.push_back(surface[ctf_sl_faces[auxndx + j]].get_area());
            totalArea += faceAreas[j];
        }
        for (int j=0; j<ctf_sl_surfsize[i]; j++) {
            resize2(faceAreas[j]/(3*totalArea), arr_view<scalar,3>(ctf_sl_forces, 3*i), traction);
            surface[ctf_sl_faces[auxndx+j]].add_force_to_node(0, traction);
            surface[ctf_sl_faces[auxndx+j]].add_force_to_node(1, traction);
            surface[ctf_sl_faces[auxndx+j]].add_force_to_node(2, traction);
        }
        auxndx += ctf_sl_surfsize[i];
        faceAreas.clear();
    }
}

/*
 */
void Blob::zero_force() {
    for (int i = 0; i < elem.size(); i++) {
        elem[i].zero_force();
    }
    for (int i = 0; i < surface.size(); i++) {
        surface[i].zero_force();
    }
}

void Blob::set_forces_to_zero() {
    for (int i = 0; i < node.size(); ++i) {
        force[i][0] = 0;
        force[i][1] = 0;
        force[i][2] = 0;
    }
}

void Blob::get_node(int index, arr3 &v) {
    store(node[index].pos, v);
}

void Blob::get_node_0(int index, arr3 &v) {
    store(node[index].pos_0, v);
}

void Blob::add_force_to_node(const arr3& f, int index) {
    force[index][0] += f[0];
    force[index][1] += f[1];
    force[index][2] += f[2];
}

void Blob::velocity_all(scalar vel_x, scalar vel_y, scalar vel_z) {
    for (auto &node_i : node) {
        node_i.vel[0] = vel_x;
        node_i.vel[1] = vel_y;
        node_i.vel[2] = vel_z;
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
        for (int n = 0; n < elem.size(); n++) {
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
void Blob::enforce_box_boundaries(arr3 &box_dim) {
    if (params.wall_x_1 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos[0] < 0) {
                if (node[i].vel[0] < 0) {
                    node[i].vel[0] = 0;
                }
            }
        }
    }
    if (params.wall_x_2 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos[0] > box_dim[0]) {
                if (node[i].vel[0] > 0) {
                    node[i].vel[0] = 0;
                }
            }
        }
    }
    if (params.wall_y_1 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos[1] < 0) {
                if (node[i].vel[1] < 0) {
                    node[i].vel[1] = 0;
                }
            }
        }
    }
    if (params.wall_y_2 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos[1] > box_dim[1]) {
                if (node[i].vel[1] > 0) {
                    node[i].vel[1] = 0;
                }
            }
        }
    }
    if (params.wall_z_1 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos[2] < 0) {
                if (node[i].vel[2] < 0) {
                    node[i].vel[2] = 0;
                }
            }
        }
    }
    if (params.wall_z_2 == WALL_TYPE_HARD) {
        for (int i = 0; i < num_surface_nodes; i++) {
            if (node[i].pos[2] > box_dim[2]) {
                if (node[i].vel[2] > 0) {
                    node[i].vel[2] = 0;
                }
            }
        }
    }
}

int Blob::get_num_nodes() {
    return static_cast<int>(node.size());
}

int Blob::get_num_elements() {
    return static_cast<int>(elem.size());
}

int Blob::get_motion_state() {
    return blob_state;
}

scalar Blob::get_scale() {
    return scale;
}

scalar Blob::get_RandU01() {
    return (*rng)[0].RandU01();
}

int Blob::get_num_linear_nodes() {
    set<int> node_indices;
    for(int n = 0; n < elem.size(); ++n) {
        for(int i = 0; i < 4; ++i) {
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

void Blob::get_min_max(arr3 &blob_min, arr3 &blob_max) {
    blob_min[0] = INFINITY;
    blob_max[0] = -1 * INFINITY;
    blob_min[1] = INFINITY;
    blob_max[1] = -1 * INFINITY;
    blob_min[2] = INFINITY;
    blob_max[2] = -1 * INFINITY;

    for(const auto &node_i : node) {
        if(node_i.pos[0] > blob_max[0]) {
            blob_max[0] = node_i.pos[0];
        } else if (node_i.pos[0] < blob_min[0]) {
            blob_min[0] = node_i.pos[0];
        }

        if(node_i.pos[1] > blob_max[1]) {
            blob_max[1] = node_i.pos[1];
        } else if (node_i.pos[1] < blob_min[1]) {
            blob_min[1] = node_i.pos[1];
        }

        if(node_i.pos[2] > blob_max[2]) {
            blob_max[2] = node_i.pos[2];
        } else if (node_i.pos[2] < blob_min[2]) {
            blob_min[2] = node_i.pos[2];
        }
    }
}

/*
 */
void Blob::load_nodes(const char *node_filename, scalar scale) {
    FILE *in;
    double x, y, z;
    const int max_line_size = 50;
    char line[max_line_size];

    // open the node file
    if (!(in = fopen(node_filename, "r"))) {
        throw FFEAFileException(node_filename);
    }
    printf("\t\tReading in nodes file: %s\n", node_filename);

    // first line should be the file type "ffea node file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of node file.");
    }
    if (strcmp(line, "walrus node file\n") != 0 && strcmp(line, "ffea node file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea node file' (read '%s') \n", line);
    }

    // read in the number of nodes in the file
    int num_nodes;
    if (fscanf(in, "num_nodes %d\n", &num_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of nodes.");
    }
    printf("\t\t\tNumber of nodes = %d\n", num_nodes);

    // read in the number of surface nodes in the file
    if (fscanf(in, "num_surface_nodes %d\n", &num_surface_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of surface nodes.");
    }
    printf("\t\t\tNumber of surface nodes = %d\n", num_surface_nodes);

    // read in the number of interior nodes in the file
    if (fscanf(in, "num_interior_nodes %d\n", &num_interior_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of interior nodes.");
    }
    printf("\t\t\tNumber of interior nodes = %d\n", num_interior_nodes);

    // Allocate the memory for all these nodes
    try {
        node = std::vector<mesh_node>(num_nodes);
        node_position = std::vector<arr3*>(num_nodes, nullptr);
    } catch (std::bad_alloc&) {
        fclose(in);
        throw FFEAException("Unable to allocate memory for nodes array.");
    }

    // Check for "surface nodes:" line
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error when looking for 'surface nodes:' line.");
    }
    if (strcmp(line, "surface nodes:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'surface nodes:' line (found '%s' instead).", line);
    }

    // Read in all the surface nodes from file
    for (int i = 0; i < node.size(); i++) {

        if (i == num_surface_nodes) {
            // Check for "interior nodes:" line
            if (!fgets(line, max_line_size, in)) {
                fclose(in);
                throw FFEAException("Error when looking for 'interior nodes:' line.");
            }
            if (strcmp(line, "interior nodes:\n") != 0) {
                fclose(in);
                throw FFEAException("Could not find 'interior nodes:' line (found '%s' instead).", line);
            }
        }

        if (fscanf(in, "%le %le %le\n", &x, &y, &z) != 3) {
            fclose(in);
            throw FFEAException("Error reading from nodes file at node %d", i);
        } else {
            node[i].pos[0] = scale * (scalar) x;
            node[i].pos[1] = scale * (scalar) y;
            node[i].pos[2] = scale * (scalar) z;
            node_position[i] = &node[i].pos;
            node[i].vel[0] = 0;
            node[i].vel[1] = 0;
            node[i].vel[2] = 0;

            node[i].index = static_cast<int>(i);
        }
    }

    fclose(in);
    printf("\t\t\tRead %zu nodes from %s\n", node.size(), node_filename);
}

/*
 */
void Blob::load_topology(const char *topology_filename) {
    FILE *in;
    const int max_line_size = 50;
    char line[max_line_size];

    // Now open the topology file
    if (!(in = fopen(topology_filename, "r"))) {
        throw FFEAFileException(topology_filename);
    }
    printf("\t\tReading in topology file: %s\n", topology_filename);

    // first line should be the file type "ffea topology file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of topology file.");
    }
    if (strcmp(line, "walrus topology file\n") != 0 && strcmp(line, "ffea topology file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea topology file' (read '%s').", line);
    }

    // read in the total number of elements in the file
    int num_elements;
    if (fscanf(in, "num_elements %d\n", &num_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of elements.");
    }
    printf("\t\t\tNumber of elements = %d\n", num_elements);

    // read in the number of surface elements in the file
    if (fscanf(in, "num_surface_elements %d\n", &num_surface_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of surface elements.");
    }
    printf("\t\t\tNumber of surface elements = %d\n", num_surface_elements);

    // read in the number of interior elements in the file
    if (fscanf(in, "num_interior_elements %d\n", &num_interior_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of interior elements.");
    }
    printf("\t\t\tNumber of interior elements = %d\n", num_interior_elements);

    // Allocate the memory for all these elements
    try {
        elem = std::vector<tetra_element_linear>(num_elements);
    } catch (std::bad_alloc &) {
        fclose(in);
        throw FFEAException("Unable to allocate memory for element array.");
    }

    // Check for "surface elements:" line
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error when looking for 'surface elements:' line.");
    }
    if (strcmp(line, "surface elements:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'surface elements:' line (found '%s' instead).", line);
    }

    // Read in all the elements from file
    int i = 0;
    for (; i < num_elements; i++) {
        if (i == num_surface_elements) {
            // Check for "interior elements:" line
            if (!fgets(line, max_line_size, in)) {
                fclose(in);
                throw FFEAException("Error when looking for 'interior elements:' line.");
            }
            if (strcmp(line, "interior elements:\n") != 0) {
                fclose(in);
                throw FFEAException("Could not find 'interior elements:' line (found '%s' instead).", line);
            }
        }

        int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
        if (fscanf(in, "%d %d %d %d %d %d %d %d %d %d\n", &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &n9, &n10) != 10) {
            fclose(in);
            throw FFEAException("Error reading from elements file at element %d", i);
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= node.size() ||
                    n2 < 0 || n2 >= node.size() ||
                    n3 < 0 || n3 >= node.size() ||
                    n4 < 0 || n4 >= node.size() ||
                    n5 < 0 || n5 >= node.size() ||
                    n6 < 0 || n6 >= node.size() ||
                    n7 < 0 || n7 >= node.size() ||
                    n8 < 0 || n8 >= node.size() ||
                    n9 < 0 || n9 >= node.size() ||
                    n10 < 0 || n10 >= node.size()) {
                fclose(in);
                throw FFEAException("Error: Element %d references an out of bounds node index.", i);
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
}

/*
 */
void Blob::load_surface(const char *surface_filename) {
    FILE *in = nullptr;
    int i, n1, n2, n3, element;
    const int max_line_size = 50;
    char line[max_line_size];

    if (!(in = fopen(surface_filename, "r"))) {
        throw FFEAFileException(surface_filename);
    }
    printf("\t\tReading in surface file: %s\n", surface_filename);

    // first line should be the file type "ffea surface file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of surface file.");
    }
    if (strcmp(line, "walrus surface file\n") != 0 && strcmp(line, "ffea surface file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea surface file' (read '%s').", line);
    }

    // read in the number of faces in the file
    int num_surface_faces;
    if (fscanf(in, "num_surface_faces %d\n", &num_surface_faces) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of faces.");
    }
    printf("\t\t\tNumber of faces = %d\n", num_surface_faces);

    // Allocate the memory for all these faces
    try {
        surface = std::vector<Face>(num_surface_faces);
    } catch (std::bad_alloc &) {
        throw FFEAException("Failed to allocate memory for the faces.");
    }

    // Check for "faces:" line
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error when looking for 'faces:' line.");
    }
    if (strcmp(line, "faces:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'faces:' line (found '%s' instead).", line);
    }

    // Read in all the faces from file
    scalar smallest_A = INFINITY;
    for (i = 0; i < num_surface_faces; i++) {
        if (fscanf(in, "%d %d %d %d\n", &element, &n1, &n2, &n3) != 4) {
            fclose(in);
            throw FFEAException("Error reading from surface file at face %d. There should be 4 space separated integers.", i);
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= node.size() ||
                    n2 < 0 || n2 >= node.size() ||
                    n3 < 0 || n3 >= node.size()) {
                fclose(in);
                throw FFEAException("Error: Surface face %d references an out of bounds node index.", i);
            } else if (element < 0 || element >= elem.size()) {
                fclose(in);
                throw FFEAException("Error: Surface face %d references an out of bounds element index.", i);
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
                throw FFEAException("Error: Could not find the opposite node.");
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
}

/*
 */
void Blob::load_surface_no_topology(const char *surface_filename) {
    FILE *in = nullptr;
    int i, n1, n2, n3, element;
    const int max_line_size = 50;
    char line[max_line_size];

    if (!(in = fopen(surface_filename, "r"))) {
        throw FFEAFileException(surface_filename);
    }
    printf("\t\tReading in surface file: %s\n", surface_filename);

    // first line should be the file type "ffea surface file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of surface file.");
    }
    if (strcmp(line, "walrus surface file\n") != 0 && strcmp(line, "ffea surface file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea surface file' (read '%s').", line);
    }

    // read in the number of faces in the file
    int num_surface_faces;
    if (fscanf(in, "num_surface_faces %d\n", &num_surface_faces) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of faces.");
    }
    printf("\t\t\tNumber of faces = %d\n", num_surface_faces);

    // Allocate the memory for all these faces
    try {
        surface = std::vector<Face>(num_surface_faces);
    } catch (std::bad_alloc &) {
        throw FFEAException("Failed to allocate memory for the faces.");
    }

    // Check for "faces:" line
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error when looking for 'faces:' line.");
    }
    if (strcmp(line, "faces:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'faces:' line (found '%s' instead).", line);
    }

    // Read in all the faces from file (element will always be zero here, because no internal structure exists)
    scalar smallest_A = INFINITY;
    for (i = 0; i < surface.size(); i++) {
        if (fscanf(in, "%d %d %d %d\n", &element, &n1, &n2, &n3) != 4) {
            fclose(in);
            throw FFEAException("Error reading from surface file at face %d. There should be 4 space separated integers.", i);
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= node.size() ||
                    n2 < 0 || n2 >= node.size() ||
                    n3 < 0 || n3 >= node.size()) {
                fclose(in);
                throw FFEAException("Error: Surface face %d references an out of bounds node index.", i);
            }

            surface[i].init(i, &node[n1], &node[n2], &node[n3], nullptr, this, params);

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
}

/*
 */
void Blob::load_material_params(const char *material_params_filename) {
    FILE *in = nullptr;
    int num_material_elements;
    const int max_line_size = 50;
    char line[max_line_size];

    if (!(in = fopen(material_params_filename, "r"))) {
        throw FFEAFileException(material_params_filename);
    }
    printf("\t\tReading in material parameters file: %s\n", material_params_filename);

    // first line should be the file type "ffea material params file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of material params file.");
    }
    if (strcmp(line, "walrus material params file\n") != 0 && strcmp(line, "ffea material params file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea material params file' (read '%s').", line);
    }

    // read in the number of elements in the file
    if (fscanf(in, "num_elements %d\n", &num_material_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of elements.");
    }
    printf("\t\t\tNumber of elements in material params file = %d\n", num_material_elements);

    // Check that we have same number of elements in material params file as in topology file
    if (num_material_elements != elem.size()) {
        fclose(in);
        throw FFEAException("Number of elements in material params file (%d) does not match number of elements in topology file (%zu).", num_material_elements, elem.size());
    }

    // Set the material parameters for each element in the Blob
    scalar density = 0.0, shear_visc = 0.0, bulk_visc = 0.0, shear_mod = 0.0, bulk_mod = 0.0, dielectric = 0.0;
    int i;
    for (i = 0; i < elem.size(); i++) {
        if (fscanf(in, "%le %le %le %le %le %le\n", &density, &shear_visc, &bulk_visc, &shear_mod, &bulk_mod, &dielectric) != 6) {
            fclose(in);
            throw FFEAException("Error reading from material params file at element %d. There should be 6 space separated real values (density, shear_visc, bulk_visc, shear_mod, bulk_mod, dielectric).", i);
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
}

/*
 */
void Blob::load_stokes_params(const char *stokes_filename, scalar scale) {
    FILE *in = nullptr;
    int num_stokes_nodes;
    const int max_line_size = 50;
    char line[max_line_size];

    if (!(in = fopen(stokes_filename, "r"))) {
        throw FFEAFileException(stokes_filename);
    }
    printf("\t\tReading in material parameters file: %s\n", stokes_filename);

    // first line should be the file type "ffea stokes radii file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of stokes radii file.");
    }
    if (strcmp(line, "walrus stokes radii file\n") != 0 && strcmp(line, "ffea stokes radii file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea stokes radii file' (read '%s').", line);
    }

    // read in the number of nodes in the file
    if (fscanf(in, "num_nodes %d\n", &num_stokes_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of nodes.");
    }
    printf("\t\t\tNumber of nodes in stokes radii file = %d\n", num_stokes_nodes);

    // Check that we have same number of nodes in stokes radii file as in nodes file
    if (num_stokes_nodes != node.size()) {
        fclose(in);
        throw FFEAException("Number of nodes in stokes radii file (%d) does not match number of nodes in nodes file (%zu).", num_stokes_nodes, node.size());
    }

    // Set the stokes radius for each node in the Blob
    scalar stokes_radius = 0.0;
    int crap, i, check = 0;
    for (i = 0; i < node.size(); i++) {
        if (fscanf(in, "%le\n", &stokes_radius) != 1) {
            fclose(in);
            throw FFEAException("Error reading from stokes radii file at node %d. There should be 1 real value (stokes radius).", i);
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
                    throw FFEAException();
                } else {
                    printf("Please enter y or n:");
                }

            }
        }
    }

    fclose(in);

    printf("\t\t\tRead %d stokes radii from %s\n", i, stokes_filename);
}

/**
 * @brief Read the beads_filename, loading beads position and types.
 *
 * @ingroup FMM
 * @details
 */
void Blob::load_beads(const char *beads_filename, scalar scale) {

    ifstream fin;
    string line;
    vector<string> vec_line;
    // int typeBead;

    fin.open(beads_filename, std::ifstream::in);
    if (fin.fail()) {
        throw FFEAFileException(beads_filename);
    }
    printf("\t\tReading in Beads file: %s\n", beads_filename);
    printf("\t\tScaling beads positions using scale: %e\n", scale);


    std::vector<string> stypes;
    string type;
    scalar x, y, z;
    bead_position.clear();

    // a set of constant strings, and a number of temporary vectors
    //    to parse the lines in search of " < nodes = ... > ".
    const string nodesKeyword = "nodes";
    const char *openField = "<";
    const char *closeField = ">";
    const char *splitNodes = ",";
    const char *defineRange = "-";
    const string defineField = "=";
    std::vector<string> v1, v2, v3, v4;

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
        bead_position.push_back({ x, y, z });

        // 1.3 - look for node restrictions within "< nodes = ...  >"
        // 1.3.1 - split the line in a vector using "<":
        boost::split(v1, line, boost::is_any_of(openField));
        bead_assignment.push_back(std::vector<int>()); // add a row.
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
                        throw FFEAException("failed to parse this set of nodes: %s", v2[1].c_str());
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
    bead_position.shrink_to_fit();

    // 2.2 - bead types are integers starting from zero:
    vector<string>::iterator it;
    try {
        bead_type = std::vector<int>(stypes.size());
    } catch (std::bad_alloc&) {
        throw FFEAException("Failed to allocate memory for bead types.");
    }
    int index;
    for (unsigned int i=0; i<stypes.size(); i++) {
        it = std::find(pc_params.types.begin(), pc_params.types.end(), stypes[i]);
        if (it == pc_params.types.end()) { // type in beads file not matching the types in .ffea file!!
            throw FFEAException("Type '%s' read in beads file does not match any of the bead types specified in the .ffea file\n", stypes[i].c_str());
        }
        index = std::distance(pc_params.types.begin(), it);
        bead_type[i] = index;
    }

    // 2.3 - num_beads:
    beads_on_blob = true;
}

void Blob::load_ctforces(const string& ctforces_fname) {

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
        throw FFEAException("ctforces_fname '%s' is too short, something is wrong\n", ctforces_fname.c_str());
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
                throw FFEAException("Invalid number of ctforces read: %s\n", line_split[1].c_str());
            }
        }

        if (line_split[0] == "num_linear_forces") {
            try {
                n_ct_lforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                throw FFEAException("Invalid number of constant linear forces read: %s\n", line_split[1].c_str());
            }
        }

        if (line_split[0] == "num_rot_forces") {
            try {
                n_ct_rforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                throw FFEAException("Invalid number of circular forces read: %s\n", line_split[1].c_str());
            }
        }

        if (line_split[0] == "num_linear_surface_forces") {
            try {
                n_cts_lforces = boost::lexical_cast<int>(line_split[1]);
            } catch (boost::bad_lexical_cast) {
                throw FFEAException("Invalid number of n_cts_lforces read: %s\n", line_split[1].c_str());
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
        throw FFEAException("--- ABORTING. Something went wrong when parsing the header of the ctforces file\n"
            "--- Check it or turn calc_ctforces to 0 in the FFEA input file");
    }


    // 2 - READ AND STORE THE LINEAR PART:
    // 2.1 - First the header:
    // check that we have a line saying "linear forces:"
    if (ctforces_lines[now_reading].compare(0, 14, "linear forces:") != 0) {
        // if not, and needed: ABORT.
        if (n_ct_lforces > 0) {
            throw FFEAException("Error reading the ctforces file; it should announce the start of linear ctforces, but instead read: %s\n", ctforces_lines[now_reading].c_str());
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
            throw FFEAException("Invalid line in the ctforces file:\n \t %s\n", ctforces_lines[i].c_str());
        }
        int b_i = boost::lexical_cast<int>(line_split[4]); // read the blob index
        if (b_i != blob_index) continue;
        int c_i = boost::lexical_cast<int>(line_split[5]); // read the conformation index
        if (c_i != conformation_index) continue;
        if (line_split[6].compare("all") == 0) {
            num_l_ctf += node.size();
        } else {
            num_l_ctf += 1;
        }
        my_lines.push_back(ctforces_lines[i]);
    }

    if (num_l_ctf > 0) cout << num_l_ctf << " linear forces were loaded for blob " << blob_index << ":" << conformation_index << endl;

    // 2.3 - And finally store the stuff properly:
    try {
        ctf_l_nodes = std::vector<int>(num_l_ctf);       // allocate nodes
        ctf_l_forces = std::vector<scalar>(3 * num_l_ctf); // allocate forces
    } catch (std::bad_alloc &) {
        throw FFEAException("Failed to allocate ctf_l relevant arrays\n");
    }
    arr3 ctf_d; // direction of the force
    int cnt = 0;
    for (int i=0; i<my_lines.size(); i++) {
        boost::split(line_split, my_lines[i], boost::is_any_of(" \t"), boost::token_compress_on);
        scalar F = boost::lexical_cast<scalar>(line_split[0]);
        ctf_d[0] = boost::lexical_cast<scalar>(line_split[1]);
        ctf_d[1] = boost::lexical_cast<scalar>(line_split[2]);
        ctf_d[2] = boost::lexical_cast<scalar>(line_split[3]);
        normalize(ctf_d);  // normalise the direction of the force.
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
            for (int j=0; j< node.size(); j++) {
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
            throw FFEAException("Wrong number of lines in ctforces. ABORTING");
        } else return; // Oh, the work was over.
    }
    // check that we have a line saying "rotational forces:"
    if (ctforces_lines[now_reading].compare(0, 18, "rotational forces:") != 0) {
        // if not, and needed: ABORT.
        if (n_ct_rforces > 0) {
            throw FFEAException("Wrong header; it should announce the start of rotational ctforces, but instead read: %s\n", ctforces_lines[now_reading].c_str());
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
            throw FFEAException("Invalid line in the ctforces file:\n \t %s\n", ctforces_lines[i].c_str());
        }
        int b_i = boost::lexical_cast<int>(line_split[8]); // read the blob index
        if (b_i != blob_index) continue;
        int c_i = boost::lexical_cast<int>(line_split[9]); // read the conformation index
        if (c_i != conformation_index) continue;
        if (line_split[10].compare("all") == 0) {
            num_r_ctf += node.size();
        } else {
            num_r_ctf += 1;
        }
        my_lines.push_back(ctforces_lines[i]);
    }

    // 3.3 - And finally store the stuff properly:
    arr3 ctf_p; // temporary point in the axis.
    try {
        ctf_r_nodes = std::vector<int>(num_r_ctf); // allocate nodes
        ctf_r_forces = std::vector<scalar>(num_r_ctf); // allocate forces
        ctf_r_axis = std::vector<scalar>(6 * num_r_ctf); // allocate axis
        ctf_r_type = std::vector<char>(2 * num_r_ctf); // allocate type of rotational force
    } catch (std::bad_alloc &) {
        throw FFEAException("Failed to allocate ctf_r relevant arrays\n");
    }
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
            normalize(ctf_d);  //  and thus normalise.
            resize(mdlm1, ctf_p);  // and rescale CTFPENDING: check!!!
        } else if (type.compare(0,1,"n")) { // otherwise store as pairs of nodes, or complain.
            throw FFEAException("Invalid rotational force: %s, in line read: %s\n", type.substr(0).c_str(), my_lines[i].c_str());
        }
        if ((type.compare(1,1,"f")) && (type.compare(1,1,"t"))) { // check and store type force or torque
            throw FFEAException("Invalid rotational force type: %s, in line read: %s\n", type.substr(1).c_str(), my_lines[i].c_str());
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
            for (int j=0; j< node.size(); j++) {
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
            throw FFEAException("Wrong number of lines in ctforces. ABORTING");
        } else return; // Oh, the work was over.
    }
    // check that we have a line saying "linear surface forces:"
    if (ctforces_lines[now_reading].compare(0, 22, "linear surface forces:") != 0) {
        // if not, and needed: ABORT.
        if (n_cts_lforces > 0) {
            throw FFEAException("Wrong header; it should announce the start of the linear surface ctforces, but instead read: %s\n", ctforces_lines[now_reading].c_str());
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
            throw FFEAException("Invalid line in the ctforces file:\n \t %s\n", ctforces_lines[i].c_str());
        }
        int b_i = boost::lexical_cast<int>(line_split[4]); // read the blob index
        if (b_i != blob_index) continue;
        int c_i = boost::lexical_cast<int>(line_split[5]); // read the conformation index
        if (c_i != conformation_index) continue;
        if (line_split[6].compare("all") == 0) {
            num_sltotal_ctf += surface.size();
            throw FFEAException();
        } else {
            num_sltotal_ctf += line_split.size() - 6;  // there may be many faces in every 'surface'!
        }
        my_lines.push_back(ctforces_lines[i]);
    }

    // 4.3 - And finally store the stuff properly:
    //     - If a face were appearing different lines, force will be applied "serially",
    //            and not summed up. They may be part of different "surfaces".
    try {
        ctf_sl_faces = std::vector<int>(num_sltotal_ctf);       // allocate faces
        num_slsets_ctf = my_lines.size();
        ctf_sl_surfsize = std::vector<int>(num_slsets_ctf);
        ctf_sl_forces = std::vector<scalar>(3 * num_slsets_ctf); // allocate forces
    } catch(std::bad_alloc &) {
        throw FFEAException("Failed to allocate ctf_sl relevant arrays\n");
    }
    ctf_d; // direction of the force
    cnt = 0;
    for (int i=0; i<my_lines.size(); i++) {
        boost::split(line_split, my_lines[i], boost::is_any_of(" \t"), boost::token_compress_on);
        scalar F = boost::lexical_cast<scalar>(line_split[0]);
        ctf_d[0] = boost::lexical_cast<scalar>(line_split[1]);
        ctf_d[1] = boost::lexical_cast<scalar>(line_split[2]);
        ctf_d[2] = boost::lexical_cast<scalar>(line_split[3]);
        normalize(ctf_d);  // normalise the direction of the force.
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
            ctf_sl_surfsize[i] = surface.size();
            for (int j=0; j< surface.size(); j++) {
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
}


void Blob::add_steric_nodes() {
    for(int i = 0; i < surface.size(); ++i) {
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
void Blob::forget_beads() {
    bead_position.clear();
    bead_type.clear();
}

void Blob::print_node_positions() {
    for (const auto &node_i : node) {
        cout << "---n: " << node_i.pos[0] << " " << node_i.pos[1] << "  " << node_i.pos[2] << endl;
    }
}

void Blob::print_bead_positions() {
    for (const auto &bead : bead_position) {
        cout << "---b: " << bead[0] << " "
             << bead[1] << " "
             << bead[2] << endl;
    }
}

/*
 */
void Blob::load_ssint(const char *ssint_filename, int num_ssint_face_types, string ssint_method) {
    FILE *in = nullptr;
    const int max_line_size = 50;
    char line[max_line_size];

    if (!(in = fopen(ssint_filename, "r"))) {
        throw FFEAFileException(ssint_filename);
    }
    printf("\t\tReading in VDW file: %s\n", ssint_filename);

    // first line should be the file type "ffea vdw file" (.vdw) or "ffea ssint file" (.ssint)
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of VDW file.");
    }
    if (strcmp(line, "walrus vdw file\n") != 0 && strcmp(line, "ffea vdw file\n") != 0 && strcmp(line, "ffea ssint file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea vdw file' (read '%s').", line);
    }

    // read in the number of faces in the file
    int num_ssint_faces = 0;
    if (fscanf(in, "num_faces %d\n", &num_ssint_faces) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of faces.");
    }
    printf("\t\t\tNumber of faces = %d\n", num_ssint_faces);

    if (num_ssint_faces != surface.size()) {
        throw FFEAException("Number of faces specified in VDW file (%d) does not agree with number in surface file (%zu).", num_ssint_faces, surface.size());
    }

    // Check for "ssint params:" line
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error when looking for 'vdw params:' line.");
    }
    if (strcmp(line, "vdw params:\n") != 0 && strcmp(line, "ssint params:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'vdw params:' line (found '%s' instead).", line);
    }

    // Read in all the ssint parameters from the file, assigning them to the appropriate faces
    int ssint_type = 0;

    // If steric only, set all to type 0
    if (params.calc_ssint == 0) {
        for(int i = 0; i < surface.size(); ++i) {
            if (fscanf(in, "%d\n", &ssint_type) != 1) {
                fclose(in);
                throw FFEAException("Error reading from VDW file at face %d. There should be 1 integer denoting ssint face species (-1 - unreactive).", i);
            } else {
                if (ssint_type > num_ssint_face_types - 1) {
                    ssint_type = 0;
                }
                surface[i].set_ssint_interaction_type(ssint_type);
            }
        }
    } else {
        for (int i = 0; i < surface.size(); i++) {
            if (fscanf(in, "%d\n", &ssint_type) != 1) {
                fclose(in);
                throw FFEAException("Error reading from VDW file at face %d. There should be 1 integer denoting ssint face species (-1 - unreactive). \n", i);
            } else {
                if (ssint_type > num_ssint_face_types - 1) {
                    throw FFEAException("Error reading from VDW file at face %d. The given VDW face type (%d) is higher than that allowed by the ssint forcefield params file (%d). \n", i, ssint_type, num_ssint_face_types - 1);
                }
                surface[i].set_ssint_interaction_type(ssint_type);
            }
        }
    }

    // Set whether vdw is active on blob
    for(int i = 0; i < surface.size(); ++i) {
        if(surface[i].ssint_interaction_type != -1) {
            ssint_on_blob = true;
            break;
        }
    }
    fclose(in);

    printf("\t\t\tRead %zu VDW faces from %s\n", surface.size(), ssint_filename);
}

/*
 */
void Blob::load_binding_sites() {
    int num_binding_site_types = binding_matrix->get_num_interaction_types();

    // Return successful as params.calc_kinetics == 0 or no sites are required
    if (s_binding_filename.empty()) return;
    // Open file
    ifstream fin;
    fin.open(s_binding_filename.c_str());
    if(fin.fail()) {
        throw FFEAException("'binding_params_fname' %s not found.", s_binding_filename.c_str());
    }

    cout << "\t\tReading in Binding Sites file: " << s_binding_filename << endl;

    // Check if correct file
    const int MAX_BUF_SIZE = 255;
    char buf[MAX_BUF_SIZE];
    string buf_string;
    vector<string> string_vec;
    fin.getline(buf, MAX_BUF_SIZE);
    buf_string = string(buf);
    boost::trim(buf_string);

    if(buf_string != "ffea binding sites file") {
        throw FFEAException("This is not a 'ffea binding site file' (read '%s').", buf);
    }

    // read in the number of binding sites in the file
    int num_binding_sites;
    fin >> buf_string >> num_binding_sites;
    cout << "\t\t\tNumber of binding sites = " << num_binding_sites << endl;

    if (num_binding_sites > surface.size()) {
        throw FFEAException("Number of binding sites specified in binding sites file (%d) cannot exceed number of surface faces (%zu).", num_binding_sites, surface.size());
    }

    if (num_binding_sites == 0) {
        return;
    }

    // Create binding sites
    try {
        binding_site = std::vector<BindingSite>(num_binding_sites);
    } catch (std::bad_alloc&) {
        throw FFEAException("Failed to allocate array of binding sites.");
    }

    // Check for "binding sites:" line
    fin.getline(buf, MAX_BUF_SIZE);
    fin.getline(buf, MAX_BUF_SIZE);
    buf_string = string(buf);
    boost::trim(buf_string);
    if(buf_string != "binding sites:") {
        throw FFEAException("Could not find 'binding sites:' line (found '%s' instead).", buf);
    }

    // Get all binding sites
    unsigned int num_faces;
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
            throw FFEAException("Unable to read type %%d num_faces %%d line for binding site %d in %s.", i, s_binding_filename.c_str());
        }

        if(bind_type >= num_binding_site_types) {
            throw FFEAException("Binding site %d specifies site type %d, which is outside range of types allowed by the 'bsite_in_fname' matrix (%d types allowed).", i, bind_type, num_binding_site_types);
        }

        binding_site[i].set_type(bind_type);
        binding_site[i].set_num_faces(num_faces);

        // Now build list of faces
        fin.getline(buf, MAX_BUF_SIZE);
        buf_string = string(buf);
        boost::trim(buf_string);
        boost::split(string_vec, buf_string, boost::is_space(), boost::token_compress_on);
        if(string_vec.size() != num_faces + 1) {
            throw FFEAException("In %s, num_faces specified, %d, != num_faces in following line, %zd.", s_binding_filename.c_str(), num_faces, string_vec.size() - 1);
        }

        for(unsigned int j = 0; j < num_faces; ++j) {
            face_index = atoi(string_vec.at(j + 1).c_str());
            if(face_index >= surface.size()) {
                throw FFEAException("Face index %d specifies face outside range of surface faces defined in surface file (%zu).", face_index, surface.size());
            } else {
                binding_site[i].add_face(&surface[face_index]);
            }
        }

        // Properties continually change so no need to calculate stuff unless about to be used
    }

    fin.close();
}

/*
 */
void Blob::load_pinned_nodes(const char *pin_filename) {
    FILE *in;
    int i, pn_index;
    const int max_line_size = 50;
    char line[max_line_size];

    if (!(in = fopen(pin_filename, "r"))) {
        throw FFEAFileException(pin_filename);
    }
    printf("\t\tReading in pinned nodes file: %s\n", pin_filename);

    // first line should be the file type "ffea pinned nodes file"
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error reading first line of pin file.");
    }
    if (strcmp(line, "walrus pinned nodes file\n") != 0 && strcmp(line, "ffea pinned nodes file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea pinned nodes file' (read '%s').", line);
    }

    // read in the number of pinned node indices in the file
    int num_pinned_nodes;
    if (fscanf(in, "num_pinned_nodes %d\n", &num_pinned_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of pinned nodes.");
    }
    printf("\t\t\tNumber of pinned nodes = %d\n", num_pinned_nodes);

    // Allocate the memory for the list of pinned node indices
    try {
        pinned_nodes_list = std::vector<int>(num_pinned_nodes);
    } catch (std::bad_alloc &) {
        throw FFEAException("Failed to allocate pinned_nodes_list.");
    }

    // Check for "pinned nodes:" line
    if (!fgets(line, max_line_size, in)) {
        fclose(in);
        throw FFEAException("Error when looking for 'pinned nodes:' line.");
    }
    if (strcmp(line, "pinned nodes:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'pinned nodes:' line (found '%s' instead).", line);
    }

    // Read in all the pinned node indices from file
    for (i = 0; i < num_pinned_nodes; i++) {
        if (fscanf(in, "%d\n", &pn_index) != 1) {
            fclose(in);
            throw FFEAException("Error reading from pinned node file at face %d. There should be 1 integer per line.", i);
        } else {
            // check that this does not reference nodes outside of the node array
            if (pn_index < 0 || pn_index >= node.size()) {
                fclose(in);
                throw FFEAException("Error: Pinned node %d references an out of bounds node index.", i);
            }
            pinned_nodes_list[i] = pn_index;
        }
    }

    fclose(in);

    if (i == 1)
        printf("\t\t\tRead 1 pinned node index from %s\n", pin_filename);
    else
        printf("\t\t\tRead %d pinned node indices from %s\n", i, pin_filename);
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
    for (int i = 0; i < elem.size(); i++) {
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

    for (int i=0; i< surface.size(); i++) {
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
    for (int i = 0; i < elem.size(); i++) {
        mass += elem[i].mass;
    }
    printf("\t\tTotal mass of blob = %e\n", mass*mesoDimensions::mass);
}

/*
 */
void Blob::aggregate_forces_and_solve() {
    // Aggregate the forces on each node by summing the contributions from each element.
#ifdef FFEA_PARALLEL_WITHIN_BLOB
    #pragma omp parallel for default(none) schedule(guided)
#endif
    for (int n = 0; n < node.size(); ++n) {
        for (int m = 0; m < node[n].num_element_contributors; m++) {
            force[n][0] += (*node[n].force_contributions[m])[0];
            force[n][1] += (*node[n].force_contributions[m])[1];
            force[n][2] += (*node[n].force_contributions[m])[2];
        }
    }

    // Aggregate surface forces onto nodes
    for (int n = 0; n < surface.size(); ++n) {
        for (int i = 0; i < 4; i++) {
            int sni = surface[n].n[i]->index;
            force[sni][0] += surface[n].force[i][0];
            force[sni][1] += surface[n].force[i][1];
            force[sni][2] += surface[n].force[i][2];

            //			printf("force on %d from face %d = %e %e %e\n", sni, n, force[sni][0], force[sni][1], force[sni][2]);
        }
    }
    //	printf("----\n\n");
    if (params.calc_stokes == 1) {
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
                for (int i = 0; i < node.size(); ++i) {
                    force[i][0] -= node[i].vel[0] * node[i].stokes_drag;
                    force[i][1] -= node[i].vel[1] * node[i].stokes_drag;
                    force[i][2] -= node[i].vel[2] * node[i].stokes_drag;
                    if (params.calc_noise == 1) {
                        force[i][0] -= RAND(-.5, .5) * sqrt((24 * params.kT * node[i].stokes_drag) / (params.dt));
                        force[i][1] -= RAND(-.5, .5) * sqrt((24 * params.kT * node[i].stokes_drag) / (params.dt));
                        force[i][2] -= RAND(-.5, .5) * sqrt((24 * params.kT * node[i].stokes_drag) / (params.dt));
                    }
                }
#ifdef FFEA_PARALLEL_WITHIN_BLOB
            }
#endif
        } else {
            if (params.calc_noise == 1) {

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
                    for (int i = 0; i < node.size(); ++i) {
                        force[i][0] -= RAND(-.5, .5) * sqrt((24 * params.kT * node[i].stokes_drag) / (params.dt));
                        force[i][1] -= RAND(-.5, .5) * sqrt((24 * params.kT * node[i].stokes_drag) / (params.dt));
                        force[i][2] -= RAND(-.5, .5) * sqrt((24 * params.kT * node[i].stokes_drag) / (params.dt));
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
    for (int n = 0; n < pinned_nodes_list.size(); ++n) {
        int pn_index = pinned_nodes_list[n];
        force[pn_index][0] = 0;
        force[pn_index][1] = 0;
        force[pn_index][2] = 0;
    }

    // This should have a similar way of making the solver matrix become the identity matrix as the pinned nodes do
    for(set<int>::iterator it = bsite_pinned_nodes_list.begin(); it != bsite_pinned_nodes_list.end(); ++it) {
        force[*it][0] = 0;
        force[*it][1] = 0;
        force[*it][2] = 0;
    }

    // Use the linear solver to solve for Mx = f where M is the Blob's mass matrix,
    // or Kv = f where K is the viscosity matrix for the system
    // x/v is the (unknown) force solution and f is the force vector for the system.
    solver->solve(force);
}

/*
 *
 */
void Blob::euler_integrate() {
    // Update the velocities and positions of all the nodes
    if (linear_solver == FFEA_NOMASS_CG_SOLVER) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none) schedule(static)
#endif
        for (int i = 0; i < node.size(); ++i) {
            node[i].vel[0] = force[i][0];
            node[i].vel[1] = force[i][1];
            node[i].vel[2] = force[i][2];

            node[i].pos[0] += force[i][0] * params.dt; // really meaning v * dt
            node[i].pos[1] += force[i][1] * params.dt;
            node[i].pos[2] += force[i][2] * params.dt;
        }

    } else {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
        #pragma omp parallel for default(none) schedule(static)
#endif
        for (int i = 0; i < node.size(); ++i) {
            node[i].vel[0] += force[i][0] * params.dt;
            node[i].vel[1] += force[i][1] * params.dt;
            node[i].vel[2] += force[i][2] * params.dt;

            node[i].pos[0] += node[i].vel[0] * params.dt;
            node[i].pos[1] += node[i].vel[1] * params.dt;
            node[i].pos[2] += node[i].vel[2] * params.dt;
        }
    }
}

/*
 *
 */
void Blob::calculate_node_element_connectivity() {
    // initialise num_element_contributors to zero
    for(auto &node_i : node)
        node_i.num_element_contributors = 0;

    // count how many times each node is referenced in the list of elements
    for (int i = 0; i < elem.size(); ++i)
        for (int j = 0; j < NUM_NODES_QUADRATIC_TET; ++j)
            node[elem[i].n[j]->index].num_element_contributors++;

    // allocate the contributions array for each node to the length given by num_element_contributors
    for (int i = 0; i < node.size(); ++i) {
        try {
            node[i].force_contributions = std::vector<arr3*>(node[i].num_element_contributors, nullptr);
        } catch (std::bad_alloc &) {
            throw FFEAException("Failed to allocate memory for 'force_contributions' array (on node %d)\n", i);
        }
    }

    // create an array of counters keeping track of how full each contributions array is on each node
    std::vector node_counter(node.size(), 0);

    // go back through the elements array and fill the contributions arrays with pointers to the
    // appropriate force contributions in the elements
    int node_index;
    for (int i = 0; i < elem.size(); i++)
        for (int j = 0; j < NUM_NODES_QUADRATIC_TET; j++) {
            node_index = elem[i].n[j]->index;
            node[node_index].force_contributions[node_counter[node_index]] = &elem[i].node_force[j];
            node_counter[node_index]++;
        }
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

void Blob::create_pinned_nodes(set<int> list) {
    try {
        pinned_nodes_list.resize(list.size());
    } catch(std::bad_alloc &) {
        throw FFEAException("Could not allocate memory for pinned_nodes_list\n");
    }
    int i = 0;
    for(auto it = list.begin(); it != list.end(); ++it) {
        pinned_nodes_list[i++] = *it;
    }
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

void Blob::build_mass_matrix() {
    // Calculate the Sparsity Pattern for the Mass matrix
    printf("\t\tCalculating sparsity pattern for 2nd order Mass matrix\n");
    SparsityPattern sparsity_pattern_mass_matrix;
    sparsity_pattern_mass_matrix.init(node.size());
    
    std::vector<MassMatrixQuadratic> M_alpha;
    try {
        M_alpha = std::vector<MassMatrixQuadratic>(elem.size());
    } catch(std::bad_alloc &) {
        throw FFEAException("Failed to allocate memory for the mass matrix\n");
    }

    scalar *mem_loc;
    int ni_index, nj_index;
    for (int el = 0; el < elem.size(); el++) {
        elem[el].construct_element_mass_matrix(M_alpha[el]);
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
