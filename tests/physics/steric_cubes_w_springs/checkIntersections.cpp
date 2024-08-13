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

#include <cstring>
#include <iostream>
#include "mat_vec_types.h"
#include "mat_vec_fns_II.h"
#include "FFEA_return_codes.h"
#include "BlobLite.h"
#include "VolumeIntersection.h"


int skipnlines(FILE *iFile, int n);
FILE *set_up_trajfile(const char *traj_filename);
int print_tet(scalar (&T)[4][3]);

FILE *set_up_trajfile(const char *traj_filename){
  
   FILE *trj; 
   char c = 0;
   if ((trj = fopen(traj_filename, "r")) == NULL) {
      cout << "Failed to open: " << traj_filename << endl; 
      return NULL; 
   }
 
   while(c != '*') {
     c = fgetc(trj);
   }  

   skipnlines(trj, 2); 

   return trj;

}

int skipnlines(FILE *iFile, int n) {

   int i=0;
   char *ignore;
   constexpr int len_crap = 256;
   char crap[len_crap];
   for (i=0; i<n; i++){
     ignore = fgets(crap, len_crap, iFile);
     //printf("crap: %s", crap); 
   }
   return 0;

}

int print_tet(scalar (&T)[4][3]){
 
   for (int i=0; i<4; i++) {
     for (int j=0; j<3; j++) {
       cout << T[i][j] << " ";
     }
     cout << endl;
   }
   // cout << endl;
   return 0; 

}

int main(int argc, char** argv) { 
 
   // usage:
   // ./find_intersections scale1 nodes_file1 topology_file1 scale2 nodes_file2 topology_file2 traj_file num_steps
   cout << "argc = " << argc << endl;
   for(int i = 0; i < argc; i++)
      cout << "argv[" << i << "] = " << argv[i] << endl;

   int num_steps = atoi(argv[8]);
 
   // initialise the blobs:
   BlobLite b1, b2;
   b1.load_nodes(argv[2], atof(argv[1])); 
   b1.load_topology(argv[3]); 
   b2.load_nodes(argv[5], atof(argv[4]));
   b2.load_topology(argv[6]); 

   // initialise the trajectory file:
   FILE *trj = set_up_trajfile(argv[7]);
   if (trj == NULL) return 1;

   int checks = 0;
   int i_vol = 0;
   int i_cms = 0;
   arr3 cm1, cm2, aux; 
   scalar d0 = 1000; 
   scalar d = 10*d0;
   // for every time step:
   for (int t=0; t<num_steps; t++) {
     // load the coordinates of this time step:
     if (b1.read_nodes_from_file(trj)) break; 
     skipnlines(trj,1);
     if (b2.read_nodes_from_file(trj)) break; 
     skipnlines(trj,6);
     if (t < 3) continue;

     b1.center_of_coord(cm1);
     b2.center_of_coord(cm2);
     d = arr3arr3Distance<scalar,arr3>(cm2, cm1);
     if ( (d - d0) > 1e-4) {
       cout << " CMs separating... something went wrong! t: " << t << endl; 
       cout << " cm1: " << cm1[0] << endl;
       cout << " cm2: " << cm2[0] << endl; 
       cout << " d: " << d << " d0: " << d0 << endl; 
       return 1;
     } 
     d0 = d;

     // calculate if there is any tetrahedra intersecting:
     checks = 0;
     i_vol = 0; 
     scalar tet1[4][3], tet2[4][3];
     for (int i=0; i<b1.num_elements; i++){ 
       // set up tetrahedron 1:
       for (int nn = 0; nn<4; nn++) {
         for (int nc = 0; nc<3; nc++) {
           tet1[nn][nc] = b1.coord[ b1.elem[ i*NUM_NODES_QUADRATIC_TET + nn ]*3 + nc ];
         }
       }
       for (int j=0; j<b2.num_elements; j++){
         // set up tetrahedron 2:
         for (int nn = 0; nn<4; nn++) {
           for (int nc = 0; nc<3; nc++) {
             tet2[nn][nc] = b2.coord[ b2.elem[ j*NUM_NODES_QUADRATIC_TET + nn ]*3 + nc ];
           }
         }
         checks += 1;
         if (volumeIntersection<scalar,arr3>(tet1, tet2, false, aux)){
           // cout << " ivol " << endl;
           i_vol += 1;
         }
       }   // elements in b2 
     }     // elements in b1
   }       // num_steps
   fclose(trj);

   cout << "There were vol:" << i_vol << " in " << checks << " checkings" << endl; 
   if (( ((float) i_vol) / checks ) > 0.00500 ) {
     cout << " There were too many intersections: " << (( (float) i_vol) / checks ) << endl; 
     return 1;
   } 

   return 0; 

} 


 
     // scalar *tet3[4], tet4[4][3];
         // tet3[nn] = &b1.coord[ b1.elem[ i*NUM_NODES_QUADRATIC_TET + nn ]*3];
       /*for (int nn = 0; nn<4; nn++) {
         cout << "tet1: n" << nn << ": " << tet1[nn][0] << " " << tet1[nn][1] << " " << tet1[nn][2] << endl; 
         cout << "tet3: n" << nn << ": " << tet3[nn][0] << " " << tet3[nn][1] << " " << tet3[nn][2] << endl; 
       }*/
