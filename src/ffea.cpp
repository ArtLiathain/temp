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

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cctype>
#include <omp.h>

#include <boost/program_options.hpp>
#include "FFEA_return_codes.h"
#include "FFEA_user_info.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
//#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"
#include "ffea_test.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace b_po = boost::program_options;
namespace fs = std::filesystem;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
#ifdef USE_MPI
    MPI::Init();
#endif
    cout << "\n\n\n***************************************************\n\tFLUCTUATING FINITE ELEMENT ANALYSIS\n***************************************************\n\n" << endl;
    //cout << " Version:\t" << FFEA_VERSION << " [" << FFEA_MASCOT << "]" << endl;
    // cout << " Version:\tSuper Saiyan " << FFEA_MASCOT << " " << FFEA_VERSION << "(Version " << FFEA_VERSION << ")" << endl;
    cout << "\tTeam:   Ryan Cocking (bsrctb@leeds.ac.uk), Molly Gravett (bsmgr@leeds.ac.uk), Ben Hanson (B.S.Hanson@leeds.ac.uk), Oliver Harlen, Sarah Harris," << endl;
    cout << "\t        Joanna Leng, Robin Oliver, Daniel Read, Robin Richardson (pyrar@leeds.ac.uk), Tom Ridley (tomridley08@gmail.com)," << endl;
    cout << "\t        Jarvellis Rogers (J.F.Rogers1@leeds.ac.uk), Albert Solernou (a.solernou@leeds.ac.uk), Rob Welch (hello@robwel.ch)\n" << endl;     

    print_ffea_version(stdout);
    print_ffea_compilation_details(stdout);

    // Get some arguments using boost
    b_po::options_description desc("Allowed options");
    b_po::positional_options_description pdesc;
    b_po::variables_map var_map;

    string script_fname;
    int mode = 0;
    int frames_to_delete = 0;
    int verbose;
#ifdef USE_MPI
    double st,et,st1,et1;

    st1=MPI::Wtime();
#endif

    // Options for visible and non-visible cmd line params
    desc.add_options()
    ("help,h", "Print usage message")
    ("no-detail,d", "Stop measurements being outputted at a higher level of detail (to a .fdm file)")
    ("input-file,i", b_po::value<string>(&script_fname), "Input script filename")
    ("delete-frames,l", b_po::value<int>(&frames_to_delete)->default_value(0), "If restarting a simulation, this will delete the final 'arg' frames before restarting")
    ("mode,m", b_po::value<int>(&mode)->default_value(0), "Simulation Mode\n\t0 - FFEA\n\t1 - Elastic Network Model\n\t2 - Dynamic Mode Model\n\t3 - Timestep Calculator\n")
    ("verbose,v", b_po::value<int>(&verbose)->default_value(0), "Prints extra details to stdout on what FFEA is doing\n\t0 - Low\n\t1 - Medium\n\t2 - High\n\t3 - Manic")
    ;

    // 1 input file max! Option invisible (positional arg)
    pdesc.add("input-file", 1);

    // Parse
    try {
        b_po::store(b_po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), var_map);
        b_po::notify(var_map);

        // User related problems end here
    } catch(exception& e) {
        cerr << e.what() << endl;
    }

    // Boost makes sure options are valid. We check arguments are valid

    // Verbosity
    set_verbosity_level(verbose);

    // --mode
    if(var_map.count("mode")) {
        if(mode != 0 && mode != 1 && mode != 2 && mode != 3 && mode != 4) {
            FFEA_error_text();
            cout << "Unrecognised argument to --mode" << endl;
            cout << desc << endl;
            return FFEA_ERROR;
        }
    }

    // --delete_frames
    if(var_map.count("delete_frames")) {
        if(frames_to_delete < 0) {
            FFEA_error_text();
            cout << "Argument to --delete_frames cannot be less than zero. We can't magic up simulation frames for you! :D" << endl;
            cout << desc << endl;
            return FFEA_ERROR;
        }
    }

    /* Starting actual functionality */

    // Help text is built in to boost
    if (var_map.count("help")) {
        cout << desc << endl;
        return FFEA_OK;
    }

    // Check we have an input script
    if (! var_map.count("input-file")) {
        cout << "\n\nUsage: ffea [FFEA SCRIPT FILE (.ffea)] [OPTIONS]\n\n\n" << endl;
        cout << desc << endl;
        return FFEA_ERROR;
    }

    // set up a script_fname with the absolute path
    fs::path fs_script_fname = fs::absolute(script_fname);
    fs::path canonicalPath = fs::canonical(fs_script_fname.parent_path());
    fs_script_fname = canonicalPath / fs_script_fname.filename();
    script_fname = fs_script_fname.string();
    cout << "\tInput FFEA script - " << script_fname << "\n\n";

    // todo: check to see if the input file is a test script, if so, return the result from the test it specifies.
    std::string testscript_extension = "ffeatest";
    if((script_fname.size()>=testscript_extension.size()) && equal(testscript_extension.rbegin(), testscript_extension.rend(), script_fname.rbegin())){
        cout << "Test script detected. Entering test mode.\n";
        int test_return_value = ffea_test::do_ffea_test(script_fname);
        return test_return_value;
    }

    //The system of all proteins, electrostatics and water
    World *world;

    // Allocate the world
    world = new World();

    // Initialise the world, loading all blobs, parameters, electrostatics, kinetics etc.
    cout << "Initialising the world:\n" << endl;
    if(world->init(script_fname, frames_to_delete, mode, !var_map.count("no-detail")) == FFEA_ERROR) {
        FFEA_error_text();
        cout << "Errors during initialisation mean World cannot be constructed properly." << endl;

        // Delete the world (oh no!)
        cout << "Deleting world..." << endl;
        delete world;
        cout << "...done. World has been sucessfully destroyed." << endl;

        return FFEA_ERROR;
    }

    /* World is initialised. How shall we run FFEA? */
    int myreturn = FFEA_OK;
    if(mode == 0) {

        // Full FFEA
        cout << "FFEA simulation - Running world...\n" << endl;
        if(world->run() == FFEA_ERROR) {
            FFEA_error_text();
            cout << "Error occurred when running World\n" << endl;
            myreturn = FFEA_ERROR;
        } else {
            cout << "...done\n" << endl;
            myreturn = FFEA_OK;
        }

    } else if(mode == 1 || mode == 2 || mode == 4) {

        // Some form of network model
        if(mode == 1) {

            /* Linear Elastic Model */
            cout << "\n\n\n***************************************************\n\tFFEA - Linear Elastic Model\n***************************************************\n\n";
        } else if (mode == 2) {

            /* Dynamic Mode Model */
            cout << "\n\n\n***************************************************\n\tFFEA - Dynamic Mode Model\n***************************************************\n\n";
        } else if (mode == 4) {

            /* Rotne Prager Dynamic Mode Model */
            cout << "\n\n\n***************************************************\n\tFFEA - Rotne Prager Dynamic Mode Model\n***************************************************\n\n";
        }

        // Read some user input

        // Get a blob list
        cout << "\tFirstly, lets decide which blobs to analyse!" << endl << endl;
        set<int> blobs;
        int ablob;
        string buf;
        while(true) {
            try {
                cout << "\t\tEnter an index for the blob you would like an elastic network model for, or type 'q' to finish?:";

                // Get index as string
                cin >> buf;
                if(buf.compare("q") == 0 || buf.compare("Q") == 0 || buf.compare("") == 0) {
                    cout << endl << "\tThat's all the blobs!" << endl;
                    break;
                }

                // Convert to int
                ablob = atoi(buf.c_str());

                // Check int
                if(ablob < 0 || ablob >= world->get_num_blobs()) {
                    FFEA_error_text();
                    cout << "\tPlease enter a valid blob index (0 <= x <" << world->get_num_blobs() << ")" << endl;
                    continue;
                } else {
                    blobs.insert(ablob);
                }
            } catch (exception &e) {
                FFEA_error_text();
                cout << "Exception occured:" << endl;
                cout << e.what() << endl;
            }
        }

        if(blobs.size() == 0) {
            cout << "No blobs selected. No elastic network models will be found :(" << endl;
            return FFEA_OK;
        }

        // How many modes?
        int num_modes;
        while(true) {
            try {
                cout << "\n\tHow many modes would you like to visualise?:";

                // Get index as string
                cin >> buf;

                // Convert to int
                num_modes = atoi(buf.c_str());

                // Check int
                if(num_modes <= 0) {
                    FFEA_error_text();
                    cout << "\tYou must choose a whole number i > 0 and i < 3N" << endl;
                    continue;
                }
                break;

            } catch(exception &e) {
                FFEA_error_text();
                cout << "Exception occured:" << endl;
                cout << e.what() << endl;
                cout << "\tPlease enter a whole number less than 3 * the total number of nodes in your system and greater than zero." << endl;
            }
        }

        set<int>::iterator it;
        cout << endl << "\tYou've asked for " << num_modes << " modes to be calculated for blob(s):" << endl << "\t";
        for(it = blobs.begin(); it != blobs.end(); ++it) {
            cout << *it << " ";
        }

        if(mode == 1) {

            /* Linear Elastic Model */
            cout << endl << endl << "\tBeginning the calculation of the LEMs..." << endl;
            if(world->lem(blobs, num_modes) == FFEA_ERROR) {
                cout << endl;
                FFEA_error_text();
                cout << "Problem when calculating linear elastic models." << endl;
                myreturn = FFEA_ERROR;

            } else {
                cout << "...done! " << endl;
                myreturn = FFEA_OK;
            }

        } else if (mode == 2) {

            /* Dynamic Mode Model */
            cout << endl << endl << "\tBeginning the calculation of the DMMs..." << endl;
            if(world->dmm(blobs, num_modes) == FFEA_ERROR) {
                cout << endl;
                FFEA_error_text();
                cout << "Problem when calculating dynamic mode models." << endl;
                myreturn = FFEA_ERROR;
            } else {
                cout << "...done! " << endl;
                myreturn = FFEA_OK;
            }

        } else if (mode == 4) {

            /* Rotne-Prager Dynamic Mode Model */
            cout << endl << endl << "\tBeginning the calculation of the Rotne-Prager DMMs..." << endl;
            if(world->dmm_rp(blobs, num_modes) == FFEA_ERROR) {
                cout << endl;
                FFEA_error_text();
                cout << "Problem when calculating RP dynamic mode models." << endl;
                myreturn = FFEA_ERROR;
            } else {
                cout << "...done! " << endl;
                myreturn = FFEA_OK;
            }
        }

    } else if(mode == 3) {

        /* Calculate maximum allowed timestep for system */
        cout << "\n\n\n***************************************************\n\tFFEA - Timestep Calculator\n***************************************************\n\n";
        if(world->get_smallest_time_constants() == FFEA_ERROR) {
            FFEA_error_text();
            cout << "Error occurred in 'get_smallest_time_constants()' in World\n" << endl;
            myreturn = FFEA_ERROR;
        } else {
            myreturn = FFEA_OK;
        }
    } else {
        cout << "Error. Unrecognised argument to --mode \"" << mode << "\"" << endl << endl;
        cout << desc << endl;
    }

    /* Delete the world (oh no!) */
    cout << "Deleting world..." << endl;
    if (myreturn != FFEA_ERROR){
        std::cout << "Exiting gracefully...\n";
        //delete world;
    }
#ifdef USE_MPI
    st = MPI::Wtime();
    et = MPI::Wtime()-st;
    et1 = MPI::Wtime()-st1;
    cout<< "benchmarking--------Finalising time of ffea:"<<et<<"seconds"<<endl;
    cout<< "benchmarking--------total executing time:"<<et1<<"seconds"<<endl;
#endif
    cout << "...done. World has been sucessfully destroyed." << endl;

#ifdef USE_MPI
    MPI::Finalize();
#endif
    return myreturn;
}
