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

#ifndef FFEA_INPUT_READER
#define FFEA_INPUT_READER

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"

using namespace std;

/** Class to read FFEA files (pseudo-xml)*/
class FFEA_input_reader {
	
	public:
		FFEA_input_reader();
		~FFEA_input_reader();


		/** Get all lines from ffea, strip them of whitespace and return as a vector object */ 
		int file_to_lines(string script_fname, vector<string> *output);

		/** Extract any block from the current block */ 
		int extract_block(string block_title, int block_index, vector<string> input, vector<string> *output, bool mandatory=true); 

		/** Get rvalue from block */ 
		int parse_tag(string input, string *output);

		/** Specifically return map data */ 
		int parse_map_tag(string input, int *map_indices, string *map_fname);

		/** Split string around delim and return as strings */ 
		int split_string(string input, string *output, string delim, size_t output_length);

		/** Split string around delim and return as strings vector. */ 
		int split_string(string input, vector<string> &output, string delim);

		/** Split string around delim and return as ints */ 
		int split_string(string input, int *output, string delim, size_t output_length);

		/** Split string around delim and return as scalars */ 
		int split_string(string input, scalar *output, string delim, size_t output_length);
        
	private:

		string buf_string;
		int copying;
		vector<string>::iterator string_it;
		
};

#endif
