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

#ifndef FFEA_RETURN_CODES_H_INCLUDED
#define FFEA_RETURN_CODES_H_INCLUDED

// #define FFEA_VERSION "2.0"
// #define FFEA_MASCOT "Mega Walrus"

#define FFEA_DIRECT_SOLVER		0
#define FFEA_ITERATIVE_SOLVER		1
#define FFEA_MASSLUMPED_SOLVER		2
#define FFEA_NOMASS_CG_SOLVER		3

#define FFEA_OK		0
#define FFEA_ERROR	-1
#define FFEA_CAUTION 	1

#define FFEA_FILE_ERROR_MESSG(F) {FFEA_error_text(); printf("Error opening file: %s\n", F); perror(NULL); return FFEA_ERROR;}
#define FFEA_ERROR_MESSG(...) {FFEA_error_text(); printf(__VA_ARGS__); return FFEA_ERROR;}
#define FFEA_CAUTION_MESSG(...) {FFEA_caution_text(); printf(__VA_ARGS__);}

#define FFEA_BLOB_IS_STATIC	0
#define FFEA_BLOB_IS_DYNAMIC	1
#define FFEA_BLOB_IS_FROZEN	2

#define FFEA_CONFORMATION_CHANGE	0
#define FFEA_BINDING_EVENT		1
#define FFEA_UNBINDING_EVENT		2
#define FFEA_IDENTITY_EVENT		3

#include <cstdio>

/** Prints "ERROR: " to stderr in red text */ 
void FFEA_error_text();
/** Prints "CAUTION: " to stderr in yellow text */
void FFEA_caution_text();
#endif
