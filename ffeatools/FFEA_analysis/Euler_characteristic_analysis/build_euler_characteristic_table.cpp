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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

//    3 - 2
//   /|  /|
//  0 - 1 |
//  | | | |
//  |7|_|_|6
//  |/__|/
//  4   5
//
// 01 12 23 30
// 04 15 26 37
// 45 56 67 74

int calc_euler_number_x_4(int flags[8])
{
	int edges[12][2] =	{
				{0, 1},
				{1, 2},
				{2, 3},
				{3, 0},
				{0, 4},
				{1, 5},
				{2, 6},
				{3, 7},
				{4, 5},
				{5, 6},
				{6, 7},
				{7, 4}
			};

	// V - E + F
	int V = 1, half_E = 0, quarter_F = 0, num_occupied = 0, duplicate_faces = 0;

	for(int i = 0; i < 8; i++) {
		if(flags[i] == 1) {
			num_occupied++;
		}
	}

	if(num_occupied == 0 || num_occupied == 8) {
		return 0;
	}

	// get number of links between 0 and 1
	for(int i = 0; i < 12; i++) {
		int s1 = flags[edges[i][0]], s2 = flags[edges[i][1]];
		int link_0_1 = s1 ^ s2;
		int link_1_1 = s1 & s2;
		half_E += link_0_1;
		duplicate_faces += link_1_1;
	}

	quarter_F = num_occupied * 3 - duplicate_faces * 2;

	// chi = V - E/2 + F/4
	// 4chi = 4V - 2E + F
	return 4 * V - 2*half_E + quarter_F;
}

int main(int argc, char **argv)
{
	if(argc < 2) {
		printf("Usage: ./build_euler_characteristic_table [OUTPUT C HEADER FILE]\n");
		return 0;
	}

	int num_possibilities = pow(2, 8) + 1;
    std::vector<int> euler_characteristic_lookup(num_possibilities);
	for(int i = 0; i < num_possibilities; i++) {
		int flags[8] = {(i & 1)>0, (i & 2)>0, (i & 4)>0, (i & 8)>0, (i & 16)>0, (i & 32)>0, (i & 64)>0, (i & 128)>0};
		euler_characteristic_lookup[i] = calc_euler_number_x_4(flags);
	}

	FILE *h_out = NULL;
	if((h_out = fopen(argv[1], "w")) == NULL) {
		printf("Couldn't open %s for writing\n", argv[1]);
		return 0;
	}

	fprintf(h_out, "#ifndef EULER_CHARACTERISTIC_TABLE\n#define EULER_CHARACTERISTIC_TABLE\n");
	fprintf(h_out, "\t\tint get_euler_characteristic_x_4(int flags[8])\n\t\t{\n");
	fprintf(h_out, "\t\t\tconst int euler_characteristic_x_4_lookup_table[%d] = {\n", num_possibilities);
	for(int i = 0; i < num_possibilities; i++) {
		fprintf(h_out, "\t\t\t\t%d,\n", euler_characteristic_lookup[i]);
	}
	fprintf(h_out, "\t\t\t};\n\n");
	fprintf(h_out, "\t\t\tint lookup_index = flags[0] + flags[1] * 2 + flags[2] * 4 + flags[3] * 8 + flags[4] * 16 + flags[5] * 32 + flags[6] * 64 + flags[7] * 128;\n\t\t\treturn euler_characteristic_x_4_lookup_table[lookup_index];\n\t\t}\n");
	fprintf(h_out, "#endif\n");

	fclose(h_out);
	return 0;
}
