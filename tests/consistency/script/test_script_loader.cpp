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
#include <sstream>
#include "FFEA_input_reader.h"
#include "SimulationParams.h"
#include <filesystem>

int main (void)
{
	std::filesystem::path path( std::filesystem::current_path() );
	string file_to_open = path.string();
	file_to_open += string("/2r5u_8ang.ffea");
	// const int MAX_BUF_SIZE = 255;
	string buf_string;
	SimulationParams params;
	FFEA_input_reader *ffeareader;
	ffeareader = new FFEA_input_reader();

	// Copy entire script into string
	vector<string> script_vector;
	if(ffeareader->file_to_lines(file_to_open, &script_vector) == FFEA_ERROR) {
		return(1);
	}

	// Get params section
	cout << "Extracting Parameters..." << endl;
	params.FFEA_script_filename = file_to_open;  // includes absolute path.
	if(params.extract_params(script_vector) != 0) {
		FFEA_error_text();
		printf("Error parsing parameters in SimulationParams::extract_params()\n");
		return(1);
	}
	cout << "...done!" << endl;

	// Check for consistency
	cout << "\nVerifying Parameters..." << endl;
	if(params.validate(0) != 0) {
		FFEA_error_text();
		printf("Parameters found to be inconsistent in SimulationParams::validate()\n");
		return(1);
	}
	
	return(0);
}
