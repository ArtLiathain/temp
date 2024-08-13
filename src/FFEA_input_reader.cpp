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

#include "FFEA_input_reader.h"

#include <sstream>

FFEA_input_reader::FFEA_input_reader() {
	buf_string = "";
	copying = 0;
}

FFEA_input_reader::~FFEA_input_reader() {
	buf_string = "";
	copying = 0;
}

int FFEA_input_reader::file_to_lines(string script_fname, vector<string> *script_vector) {

	// Open script
	ifstream fin;
	fin.open(script_fname.c_str());
	if(fin.fail()) {
		FFEA_error_text();
		cout << script_fname << " does not exist.\n" << endl;
		return FFEA_ERROR;
	}

	// Copy entire script into string
	script_vector->clear();

	// The following variable declarations belong to 
	//   to the comment removal functionality
	size_t found, ini, end;
	int comment = 0;
	int count, count_0;
	string m_ini = "<!--";
	string m_end = "-->";
	string buf2_string;

	while(!fin.eof()) {
      getline(fin, buf_string);

                // start removing comments enclosed in <!-- --> 
                buf2_string.clear();
                found = 0;
                count = 0;
                count_0 = 0;
                ini = 0;
                end = buf_string.length(); 
                // essentially parse the line until no more 
                while (found != string::npos) {
                  if (comment == 0) {
                    found = buf_string.find(m_ini);
                    if (found!=string::npos) {
                       count += 1;
                       comment = 1;
                       ini = found;
                    }
                  }
                  if (comment == 1) {
                    found = buf_string.find(m_end);
                    if (found!=string::npos) {
                       count += 2;
                       comment = 0;
                       end = found + 3;
                    }
                  }
                  // the line end up without closing the comment: 
                  if (comment == 1) {
                    // cout << "!!! COMMENT IN\n";
                    buf_string = buf_string.substr(0,ini);
                    break;
                  // we're out of the comment.
                  } else if (comment == 0) {
                    if (count == count_0 + 3) {
                      buf2_string = buf_string.substr(0,ini);
                      buf2_string.append( buf_string.substr(end) );
                      buf_string = buf2_string;
                      count_0 = count;
                    } else if (count == count_0 + 2) {
                      buf_string = buf_string.substr(end);
                      count_0 = count;
                    }
                  }
                }
                // end removing comments. 


		boost::trim(buf_string);
		if(buf_string != "\n" && buf_string != "") {
			script_vector->push_back(buf_string);
		}
	}
	fin.close();
	return FFEA_OK;
}

int FFEA_input_reader::extract_block(string block_title, int block_index, vector<string> input, vector<string> *output, bool mandatory/*=true*/) {

	// Immediate error checking
	if(block_title != "param" && block_title != "system" && block_title != "blob" && block_title != "conformation" && block_title != "kinetics" && block_title != "maps" && block_title != "interactions" && block_title != "springs" && block_title != "precomp" && block_title != "ctforces" && block_title != "rod" && block_title != "coupling" ) {
		FFEA_error_text();
		cout << "Unrecognised block: " << block_title << ". Block structure is:" << endl;
		cout << "<param>\n</param>\n<system>\n\t<blob>\n\t\t<conformation>\n\t\t</conformation>\n\t\t\t.\n\t\t\t.\n\t\t\t.\n\t\t<kinetics>\n\t\t\t<maps>\n\t\t\t</maps>\n\t\t</kinetics>\n\t</blob>\n\t\t.\n\t\t.\n\t\t.\n\t<interactions>\n\t\t<springs>\n\t\t</springs>\n\t\t<precomp>\n\t\t</precomp>\n\t\t<ctforces>\n\t\t</ctforces>\n\t</interactions>\n</system>" << endl;
		return FFEA_ERROR;
	}

	int count = -1;
	output->clear();
	//cout << "New" << endl << endl;
	for(string_it = input.begin(); string_it != input.end(); ++string_it) {
		buf_string = boost::erase_last_copy(boost::erase_first_copy(*string_it, "<"), ">");
		boost::trim(buf_string);
        
        //Parse an enclosing blokc with attributes (e.g. <coupling type = rod>)
        vector<string> split_buf_string;
        boost::split(split_buf_string, buf_string, boost::is_any_of(" "));
        if (split_buf_string.size() > 1){
            buf_string = split_buf_string[0];
        }
        
		if(buf_string == block_title) {
			
			if(copying == 1) {
				FFEA_error_text();
				cout << "Shouldn't have found '" << buf_string << "' within " << block_title << " block." << endl;
				return FFEA_ERROR;
			}

			count++;
			if(count == block_index) {
				copying = 1;
				continue;
			}
		}

		if(copying == 1) {
			if(buf_string == "/" + block_title) {
				copying = 0;
				return FFEA_OK;
			} else {
				output->push_back(*string_it);
			}
		}
	}

	if(copying == 1) {
		FFEA_error_text();
		cout << "Never found closing tag '/" << block_title << "'." << endl;
		return FFEA_ERROR;
	} else {
		if (mandatory) {
			FFEA_error_text();
			cout << "Specified block_index " << block_index << " for block '" << block_title << "' not found." << endl;
		}
		return FFEA_ERROR;
	}
	
}
	
/** 
 * @brief parse an input ffea line
 * @param[in] string input e. g.,  string < blah = whatever>
 * @param[out] string[2] output; string[0] = blah, string[1] = whatever.
 */
int FFEA_input_reader::parse_tag(string input, string *output) {

	string tag;
	vector<string> lrvalvec;
   if (std::count(input.begin(), input.end(), '<') == 0)
     FFEA_ERROR_MESSG("Line '%s' is missing '<'\n", input.c_str())
   if (std::count(input.begin(), input.end(), '>') == 0)
     FFEA_ERROR_MESSG("Line '%s' is missing '>'\n", input.c_str())
     
	tag = boost::erase_last_copy(boost::erase_first_copy(input, "<"), ">");
	boost::trim(tag);

	// Split around "=", trim and return
	boost::split(lrvalvec, tag, boost::is_any_of("="));
     
	for(unsigned int i = 0; i < lrvalvec.size(); ++i) {
		output[i] = lrvalvec.at(i);
		boost::trim(output[i]);
	}

	return FFEA_OK;
}

int FFEA_input_reader::parse_map_tag(string input, int *map_indices, string *map_fname) {

	// Parse whole tag
	string lrvalue[2], indexlrvalue[2];
	parse_tag(input, lrvalue);

	// Check if map
	split_string(lrvalue[0], indexlrvalue, "(", 2);
	if(indexlrvalue[0] != "map") {
		cout << indexlrvalue[0] << endl;
		FFEA_error_text();
		cout << "Expected '<map (from,to) = fname>' but got " << input << endl;
		return FFEA_ERROR;
	}

	// Assign map!
	*map_fname = lrvalue[1];

	// Get indices
	indexlrvalue[1] = boost::erase_last_copy(indexlrvalue[1], ")");
	split_string(indexlrvalue[1], map_indices, ",", 2);
	return FFEA_OK;
}

int FFEA_input_reader::split_string(string input, string *output, string delim, size_t output_length) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		if (i >= output_length) {
			std::stringstream ss;
			ss << "String \"" << input << "\" contains more than the expected " << output_length << " items, it may be malformed.";
			throw std::out_of_range(ss.str());
		}
		output[i] = *it;
		boost::trim(output[i++]);
	}

	return lrvalvec.size();
}

int FFEA_input_reader::split_string(string input, vector<string> &output, string delim) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		output.push_back(*it);
		boost::trim(output[i++]);
	}

	return lrvalvec.size();
}

int FFEA_input_reader::split_string(string input, int *output, string delim, size_t output_length) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		if (i >= output_length) {
			std::stringstream ss;
			ss << "String \"" << input << "\" contains more than the expected " << output_length << " items, it may be malformed.";
			throw std::out_of_range(ss.str());
		}
		output[i++] = atoi((*it).c_str());
	}
	return lrvalvec.size();
}

int FFEA_input_reader::split_string(string input, scalar *output, string delim, size_t output_length) {

	vector<string> lrvalvec;
	vector<string>::iterator it;
	boost::split(lrvalvec, input, boost::is_any_of(delim));

	int i = 0;
	for(it = lrvalvec.begin(); it != lrvalvec.end(); it++) {
		if (i >= output_length) {
			std::stringstream ss;
			ss << "String \"" << input << "\" contains more than the expected " << output_length << " items, it may be malformed.";
			throw std::out_of_range(ss.str());
		}
		output[i++] = atof((*it).c_str());
	}
	return lrvalvec.size();
}
