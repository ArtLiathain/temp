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

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <iterator>
#include <cstring>


using namespace std;

string GetFileExtension(const string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

class vector3
{
	public:
		
		vector3() 
		{
			this->x = 0.0;
			this->y = 0.0;
			this->z = 0.0;
		}
		
		vector3(double x, double y, double z) 
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}
		
		// Vector functions
		void normalise() {
			double l = sqrt(x * x + y * y + z * z);
			if(l > 0.0) {
				x /= l;
				y /= l;
				z /= l;
			}
		}
	
		vector3 get_normal() {
			double l = sqrt(x * x + y * y + z * z);
			if(l > 0.0) {
				return vector3(this->x/l, this->y/l, this->z/l);
	
			} else {
				return vector3(this->x, this->y, this->z);
			}
		}
		double get_magnitude() {
			return sqrt(x * x + y * y + z * z);
		}
		
		void translate(vector3 a) {
			this->x += a.x;
			this->y += a.y;
			this->z += a.z;
			
		}
		
		void set_midpoint(vector3 a, vector3 b) {
			 this->x = (a.x + b.x) / 2.0;
			 this->y = (a.y + b.y) / 2.0;
			 this->z = (a.z + b.z) / 2.0;
		}
		
		void set_pos(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
			
		}
		
		
		// Overloaded operators
		vector3 operator+(vector3 b) {
			return vector3(x + b.x, y + b.y, z + b.z);
		}
		
		vector3 operator-(vector3 b) {
			return vector3(x - b.x, y - b.y, z - b.z);
		}
		
		vector3 operator*(double b) {
			return vector3(x * b, y * b, z * b);
		}
		void operator+=(vector3 b) {
			x += b.x;
			y += b.y;
			z += b.z;
		}
		
		void operator*=(double b) {
			x *= b;
			y *= b;
			z *= b;
		}
		
		static vector3 cross_product(vector3 a, vector3 b) {
			vector3 c;
			c.x = a.y * b.z - a.z * b.y;
			c.y = b.x * a.z - a.x * b.z;
			c.z = a.x * b.y - a.y * b.x;
			return c;
		}
		
		static double dot_product(vector3 a, vector3 b) {
			return a.x * b.x + a.y * b.y + a.z * b.z;
		}
		
		// Vector components
		double x;
		double y;
		double z;
};

class Edge
{
	public:
	
		Edge()
		{
			this->n_index[0] = 0;
			this->n_index[1] = 0;
			this->length = 0.0;	
		}
		
		Edge(int n0, int n1, vector3 *node)
		{
			this->n_index[0] = n0;
			this->n_index[1] = n1;
			this->length = calc_length(node);
		}
		
		// Edge functions
		int init(int n0, int n1, vector3 *node) {
			this->n_index[0] = n0;
			this->n_index[1] = n1;
			this->length = calc_length(node);
			return 0;
		}

		double calc_length(vector3 *node) {
			vector3 v(node[this->n_index[1]].x - node[this->n_index[0]].x, node[this->n_index[1]].y - node[this->n_index[0]].y, node[this->n_index[1]].z - node[this->n_index[0]].z);
			this->length = v.get_magnitude();
			return v.get_magnitude();
		}
		
		vector3 calc_centroid(vector3 *node) {
			vector3 centroid;
			centroid.x = (node[n_index[0]].x + node[n_index[1]].x) / 2.0;
			centroid.y = (node[n_index[0]].y + node[n_index[1]].y) / 2.0;
			centroid.z = (node[n_index[0]].z + node[n_index[1]].z) / 2.0;
			
			return centroid;
		}
	
		bool within_limits(double *limits, vector3 *node) {
			vector3 centroid = calc_centroid(node);
			if(centroid.x >= limits[0] && centroid.x <= limits[1] && centroid.y >= limits[2] && centroid.y <= limits[3] && centroid.z >= limits[4] && centroid.z <= limits[5]) {
				return true;
			} else {
				return false;
			}
		}

		// Edge components
		int n_index[2];
		double length;
};
class Face
{
	public:
		
		Face(int n0, int n1, int n2, vector3 *node)
		{
			this->n[0] = n0;
			this->n[1] = n1;
			this->n[2] = n2;
			this->remaining_node = n0;
			edge = new Edge[3];
			edge[0].init(n0, n1, node);
			edge[1].init(n1, n2, node);
			edge[2].init(n2, n0, node); 
			calc_shortest_edge();
		}
		
		// Face functions
		void calc_shortest_edge() {
			shortest_edge = INFINITY;
			for(int i = 0; i < 3; ++i) {
				if(edge[i].length < shortest_edge) {
					shortest_edge = edge[i].length;
					shortest_edge_index = i;
				}
			}
		}
		
		double get_edge_length(int n_a, int n_b) {
			for(int i = 0; i < 3; ++i) {
				if((edge[i].n_index[0] == n_a && edge[i].n_index[1] == n_b) || (edge[i].n_index[0] == n_b && edge[i].n_index[1] == n_a)) {
					return edge[i].length;
				}
			}
		}
		
		bool same_as(Face *f) {
			int i, j, count = 0;
			for(i = 0; i < 3; ++i) {
				for(j = 0; j < 3; ++j) {
					if(f->n[i] == n[j]) {
						count++;
						break;
					}
				}
			}
			if(count == 3) {
				return true;
			} else {
				return false;
			}
		}

		bool within_limits(double *limits, vector3 *node) {
			
			if(edge[shortest_edge_index].within_limits(limits, node)) {
				return true;
			} else {
				return false;
			}
		}

		// Face components
		int n[3];
		int remaining_node;
		Edge *edge;
		double shortest_edge;
		int shortest_edge_index;
};

class Tetrahedron
{
	public:
		
		Tetrahedron(vector3 *n0, vector3 *n1, vector3 *n2, vector3 *n3) {
			this->n[0] = n0;
			this->n[1] = n1;
			this->n[2] = n2;
			this->n[3] = n3;
			volume = 0.0;
		}
		
		~Tetrahedron() {
			for(int i = 0; i < 4; ++i) {
				n[i] = NULL;
			}
			volume = 0.0;
		}
		
		double calc_volume() {
			
			// Calc new edges
			edge[0] = *n[1] - *n[0];
			edge[1] = *n[2] - *n[0];
			edge[2] = *n[3] - *n[0];
			vector3 work_vec;
			work_vec = vector3::cross_product(edge[1], edge[0]);
			volume = vector3::dot_product(edge[2], work_vec);
			volume /= 6.0;
			return volume;
		}
	
		vector3 *n[4];
		// Only require 3 edges for volume
		vector3 edge[3];
		double volume;
};

// Not a proper volume, just the bare minimum for this implementation i.e. we don't care about it's connection consistency
class Volume
{
	public:
		
		Volume() 
		{
			element.clear();
			num_elements = 0; 
			volume = 0.0;
		}
		
		~Volume() 
		{
			element.clear();
			num_elements = 0; 
			volume = 0.0;
		}
		
		void reset() {
			element.clear();
			num_elements = 0;
			volume = 0.0;
		}
		
		void build_from_partial_surf(list <Face*> partial_surf, vector3 *new_node, vector3 *node) {
			list<Face*>::iterator face_iterator;
			
			for(face_iterator = partial_surf.begin(); face_iterator != partial_surf.end(); face_iterator++) {
				element.push_back(new Tetrahedron(&node[(*face_iterator)->n[0]], &node[(*face_iterator)->n[1]], &node[(*face_iterator)->n[2]], new_node));
			}
		}
		
		double calc_volume() {
			volume = 0.0;
			list<Tetrahedron*>::iterator elem_iter;
			for(elem_iter = element.begin(); elem_iter != element.end(); elem_iter++) {
				volume += (*elem_iter)->calc_volume();
			}
			
			return volume;
		}
		
	private:
			
		list<Tetrahedron*> element;
		int num_elements;
		double volume;	
};


// Main class defining surface
class Surface
{
    /**
     * Replacement for unistd.h getline()
     *
     * reads an entire line from stream, storing the address
     * of the buffer containing the text into *lineptr.  The buffer is
     * null-terminated and includes the newline character, if one was
     * found.
     */
    size_t getline(char** lineptr, size_t* n, FILE* stream) {
        // First pass, we should do a false run first to count
		long i = 0;
		char c;
		do {
			c = getc(stream);
			if (c == EOF)
				break;
			++i;
		} while (c != '\n');
		// if lineptr needs malloc/realloc perform that
		if (!*lineptr || *n < i + 1) {
			*n = i;
			*lineptr = static_cast<char*>(realloc(*lineptr, i + 1));
			memset(*lineptr, 0, i + 1);
		}
		// Reset the stream ptr	
		if(fseek(stream, ftell(stream) - i, SEEK_CUR)) {
			printf("getline() failed\n");
		}
		// Second pass, copy the data
		size_t j = 0;
		do {
			c = getc(stream);
			if (c == EOF)
				break;
			*lineptr[j++] = c;
		} while (c != '\n');
		// Return the length of the string, without terminating char
		return j;
	}
	public:
		
		Surface() 
		{
			node = NULL;
			num_nodes = 0;
			num_faces = 0;
			num_nodes_initial = 0;
			num_faces_initial = 0;
		}
		
		~Surface()
		{
			delete[] node;
			node = NULL;
			num_nodes = 0;
			num_faces = 0;
			num_nodes_initial = 0;
			num_faces_initial = 0;
		}
		
		// Surface functions
		int init(char *surf_fname) {
		
			int i = 0;
			string fname(surf_fname), ext;
			ext = GetFileExtension(fname);

			// Open file
			FILE *surf_file;
			surf_file = fopen(surf_fname, "r");

			if(surf_file == NULL) {
				printf("Error. Filename %s does not exist.\n", surf_fname);
				return -1;
			}

			// Read differently depending on extension
			if(ext.compare("surf") == 0) {
			
				// Check file type is correct
				char line[255];
				if(fscanf(surf_file, "%s\n", line) != 1) {
					printf("Error. 'surfacemesh' line does not exist.\n");
					return -1;
				}
				if(strcmp(line, "surfacemesh") != 0) {
					printf("Error. File %s is not a NETGEN surface file\n", surf_fname);
					return -1;
				}
			
				// Get num_nodes
				if(fscanf(surf_file, "%d\n", &num_nodes) != 1) {
					printf("Error. 'num_nodes' line does not exist.\n");
					return -1;
				};
				num_nodes_initial = num_nodes;
			
				// Create nodes
				node = new vector3[num_nodes];
			
				// Read nodes
				for(i = 0; i < num_nodes; ++i) {
					if(fscanf(surf_file, "%lf %lf %lf\n", &node[i].x, &node[i].y, &node[i].z) != 3) {
						printf("Error. Expected %d nodes, found only %d.\n", num_nodes, i);
					}
				}
			
				// Get num_faces
				if(fscanf(surf_file, "%d\n", &num_faces) != 1) {
					printf("Error. 'num_faces' line does not exist.\n");
					return -1;
				};
				num_faces_initial = num_faces;
			
				// Read faces
				int n0, n1, n2;
				for(i = 0; i < num_faces; ++i) {
					if(fscanf(surf_file, "%d %d %d\n", &n0, &n1, &n2) != 3) {
						printf("Error. Expected %d faces, found only %d.\n", num_faces, i);
					}
				
					// Add to list
					face.push_back(new Face(n0 - 1, n1 - 1, n2 - 1, node));
				}
			
			} else if (ext.compare("obj") == 0) {

				// Check file type is correct
				size_t len;
				char line_t, *line2 = new char[255];
				list<vector3> verts;
				list<int> fs;
				vector3 v;
				int n[3];
				long int loc;
				int ignore;

				line_t = getc (surf_file);
				while(line_t != EOF) {
					
					if(line_t == 'v' || line_t == 'V') {

						// vn is a normal. Not needed
						loc = ftell(surf_file);
						line_t = getc (surf_file);
						if(line_t == 'n' || line_t == 'N') {
							ignore = getline(&line2, &len, surf_file);
						} else {
							fseek(surf_file, loc, SEEK_SET);
							ignore = fscanf(surf_file, "%lf %lf %lf\n", &v.x, &v.y, &v.z);
							verts.push_back(v);
							num_nodes++;
						}


					} else if (line_t == 'f' || line_t == 'F') {

						// Sometimes there is weird formatting :/
						loc = ftell(surf_file);
						if(fscanf(surf_file, "%d %d %d\n", &n[0], &n[1], &n[2]) != 3) {
							fseek(surf_file, loc, SEEK_SET);
							ignore = fscanf(surf_file, " %d//%*d %d//%*d %d//%*d\n", &n[0], &n[1], &n[2]);
						}

						// Add to list
						fs.push_back(n[0]);
						fs.push_back(n[1]);
						fs.push_back(n[2]);
						num_faces++;
					} else {
						ignore = getline(&line2, &len, surf_file);
					}
					line_t = getc (surf_file);
				}
				//delete[] line2;

				// Create nodes and faces and add them from the vectors

				// Build node list
				num_nodes_initial = num_nodes;
				node = new vector3[num_nodes];
				i = -1;
				for(list<vector3>::iterator it = verts.begin(); it != verts.end(); ++it) {
					i += 1;
					node[i].x = (*it).x;
					node[i].y = (*it).y;
					node[i].z = (*it).z;
				}
				
				// Build face list
				for(i = 0; i < num_faces; ++i) {
					n[0] = fs.front();
					fs.pop_front();
					n[1] = fs.front();
					fs.pop_front();
					n[2] = fs.front();
					fs.pop_front();
					face.push_back(new Face(n[0] - 1, n[1] - 1, n[2] - 1, node));
				}

			} else if (ext.compare("stl") == 0) {
				cout << "Error. Currently do not read stl files. Please use an alternate format.\nIf you are making intermediate structures, only use .stl as your last output. Use .surf for intermediaries." << endl;
				cout << flush;				
				return 1;
				/*// Connectivity not present in stl files. Will need to sort it out :(
				list<vector3> tempnodes;
				size_t len;
				char *line = new char[255];
				
				getline(&line, &len, surf_file);*/
			}

			// Close file and return 
			fclose(surf_file);
			return 0;
		}
		
		int coarsen(double limit, char *conserve_vol, char *find_smallest, double *limits) {
			
			// Iterate through faces while any edge is less than limit
			int i, j, k, runs = 0, breaks = 0, completed_check, num_shared_faces, num_shared_nodes;
			int node_to_delete[2], remaining_node, do_not_delete, twod_section_exists, face_sharing_remaining_node[2];
			double length_deleted, original_volume, new_volume, tolerance;
			vector3 midpoint, connecting_node, distance_node, normal_vector, translate_vector, upper_limit_pos, lower_limit_pos, nodal_distance, nodal_distance_2;
			list<Face*>::iterator face_iterator, face_to_delete[4];
			for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
				//cout << (*face_iterator)->n[0] << " " << (*face_iterator)->n[1] << " " << (*face_iterator)->n[2] << endl;
			}
			while(true) {
				
				// Output completion data
				// Output Info
				if((num_faces_initial - num_faces) % 10 == 0) {
					printf("\rFaces Remaining = %d\tNodes Remaining = %d\tLast Length Deleted = %f", num_faces, num_nodes, length_deleted);
					//printf("Nodes Remaining = %d\n", num_nodes);
					//printf("Last Length Deleted = %f\n\n", length_deleted);					
				}

				runs++;
				completed_check = 1;

				// Find an edge requiring deletion, flag face and two nodes
				if(strcmp(find_smallest, "y") == 0 || strcmp(find_smallest, "Y") == 0) {

					// Find smallest edge in system first and face containing it
					length_deleted = INFINITY;
					for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
						if((*face_iterator)->shortest_edge < length_deleted) {

							// Check it's within limits
							if((*face_iterator)->within_limits(limits, node)) {
								length_deleted = (*face_iterator)->shortest_edge;
								face_to_delete[0] = face_iterator;
							}
						}
					}
					//printf("%d %d %d\n", (*face_to_delete[0])->n[0], (*face_to_delete[0])->n[1], (*face_to_delete[0])->n[2]);
					//printf("%e %e %d\n", length_deleted, limit, length_deleted < limit);
					// Not finished if length < limit
					if(length_deleted < limit) {
						completed_check = 0;
					}
		
					// Get nodes to delete
					for(i = 0; i < 2; ++i) {
						node_to_delete[i] = (*face_to_delete[0])->edge[(*face_to_delete[0])->shortest_edge_index].n_index[i];
					}
					
					// Get remaining node in that face
					for(i = 0; i < 3; ++i) {
						if((*face_to_delete[0])->n[i] != node_to_delete[0] && (*face_to_delete[0])->n[i] != node_to_delete[1]) {
							(*face_to_delete[0])->remaining_node = (*face_to_delete[0])->n[i];
						}
					}
				} else {

					// Find any edge less than limit and face contatining it
					for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
						if((*face_iterator)->shortest_edge < limit) {
							length_deleted = (*face_iterator)->shortest_edge;
							face_to_delete[0] = face_iterator;
							
							// Get nodes to delete
							for(i = 0; i < 2; ++i) {
								node_to_delete[i] = (*face_iterator)->edge[(*face_iterator)->shortest_edge_index].n_index[i];
							}
							for(i = 0; i < 3; ++i) {
								if((*face_iterator)->n[i] != node_to_delete[0] && (*face_iterator)->n[i] != node_to_delete[1]) {
									(*face_iterator)->remaining_node = (*face_iterator)->n[i];
								}
							}
							completed_check = 0;
							break;
						}
					}
				}

				// Break loop if completed for remapping nodes to create a consistent surface
				if(completed_check == 1) {
					printf("\nCoarsening Completed.\n");
					break;
				}
				// Find the other face and the neighbourhood about the nodes
				// Empty previous neighbourhood
				neighbourhood.clear();

				// Find faces and neighbourhood based on how many of the nodes_to_delete they are connected to
				num_shared_faces = 0;
				for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
					
					num_shared_nodes = 0;
					for(i = 0; i < 3; ++i) {
						for(j = 0; j < 2; ++j) {
							if((*face_iterator)->n[i] == node_to_delete[j]) {
								num_shared_nodes++;
								break;
							}
						}
					}

					// If only shares one node, counts as a neighbourhood face
					if(num_shared_nodes == 1) {
						neighbourhood.push_back(*face_iterator);
					}

					// If shares two nodes, shares the whole edge. These faces need dealing with (may be more than two)
					// Will find original face, then overwrite with the second
					else if(num_shared_nodes == 2) {
					
						// Remember how many faces are shared (2 is easy, 4 is not. Others are weird :/)
						num_shared_faces++;

						// Remember remaining node
						for(i = 0; i < 3; ++i) {
							if((*face_iterator)->n[i] != node_to_delete[0] && (*face_iterator)->n[i] != node_to_delete[1]) {
								(*face_iterator)->remaining_node = (*face_iterator)->n[i];
							}
						}
						face_to_delete[num_shared_faces - 1] = face_iterator;
					//} else {
						//printf("Num shared node = %d\n", num_shared_nodes);
					}
				}

				// Coarsen in different ways depending on how many faces are sharing the edge
	
				if(num_shared_faces == 2) {
					
					// To conserve volume, make a local volume mesh before moving the midpoint and calculate the volume
					if(strcmp(conserve_vol, "y") == 0 || strcmp(conserve_vol, "Y") == 0) {
						local_volume_surf.clear();
						for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
							local_volume_surf.push_back(*face_iterator);
						}
						local_volume_surf.push_back(*face_to_delete[0]);
						local_volume_surf.push_back(*face_to_delete[1]);
					
						// Also, calculate the vector along which node will be translated to conserve volume and shape
						for(i = 0; i < 2; ++i) {
							normal_vector = vector3::cross_product(node[(*face_to_delete[i])->n[1]] - node[(*face_to_delete[i])->n[0]], node[(*face_to_delete[i])->n[2]] - node[(*face_to_delete[i])->n[1]]);
							normal_vector.normalise();
							translate_vector += normal_vector;
							
						}
						translate_vector.normalise(); 
						
						// Create the original volume mesh
						connecting_node.set_midpoint(node[node_to_delete[0]], node[node_to_delete[1]]);
						connecting_node.translate(translate_vector * -100);
						distance_node.set_midpoint(node[node_to_delete[0]], node[node_to_delete[1]]);
						distance_node.translate(translate_vector * 100);
						local_volume.reset();
						local_volume.build_from_partial_surf(local_volume_surf, &connecting_node, node);
						original_volume = local_volume.calc_volume();
					}
					
					// For later
					nodal_distance = node[node_to_delete[1]] - node[node_to_delete[0]];
					
					// Move lowest node to midpoint
					midpoint.set_midpoint(node[node_to_delete[0]], node[node_to_delete[1]]);
					node[node_to_delete[0]].set_pos(midpoint.x, midpoint.y, midpoint.z);

					// Move highest node to infinity
					node[node_to_delete[1]].set_pos(INFINITY, INFINITY, INFINITY);
					num_nodes--;

					// Reconnect all neighbourhood faces to other node and recalculate shortest edge
					for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
					
						// Node indices
						for(i = 0; i < 3 ; ++i) {
							if((*face_iterator)->n[i] == node_to_delete[1]) {
								(*face_iterator)->n[i] = node_to_delete[0];
							}
						}	
					
						// Edge lengths
						for(i = 0; i < 3; ++i) {
							for(j = 0; j < 2; ++j) {
								if((*face_iterator)->edge[i].n_index[j] == node_to_delete[1]) {
									(*face_iterator)->edge[i].n_index[j] = node_to_delete[0];
									break;	
								}
							}
							(*face_iterator)->edge[i].calc_length(node);
						}
					
						// Shortest length
						(*face_iterator)->calc_shortest_edge();
					}
					
					// Delete faces
					face.erase(face_to_delete[1]);
					face.erase(face_to_delete[0]);
					num_faces -= 2;

					// To conserve volume, make a new volume mesh after moving the midpoint and iterate new node position until volumes are equal
					if(strcmp(conserve_vol, "y") == 0 || strcmp(conserve_vol, "Y") == 0) {
						
						// Create new volume mesh
						local_volume_surf.clear();
						for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
							local_volume_surf.push_back(*face_iterator);	
						}
						local_volume.reset();
						local_volume.build_from_partial_surf(local_volume_surf, &connecting_node, node);
						
						// Iterate node position until volume is equal within tolerance
						// Firstly, calculate lower and upper limit volumes
						tolerance = INFINITY;
						node[node_to_delete[0]].set_pos(distance_node.x, distance_node.y, distance_node.z);
						upper_limit_pos.set_pos(distance_node.x, distance_node.y, distance_node.z);
						node[node_to_delete[0]].set_pos(connecting_node.x, connecting_node.y, connecting_node.z);
						lower_limit_pos.set_pos(connecting_node.x, connecting_node.y, connecting_node.z);
						int while_counts = 0;
						while(tolerance > 0.01) {
							while_counts++;
							//printf("%e %e %e\n", connecting_node.x, connecting_node.y, connecting_node.z);
							if(while_counts > 100) {
								node[node_to_delete[0]].set_pos(midpoint.x, midpoint.y, midpoint.z);
								break;
							}
							node[node_to_delete[0]].set_midpoint(lower_limit_pos, upper_limit_pos);
							new_volume = local_volume.calc_volume();
							// If this is true, volume is too small for double precision. Makes both limit 0.0 :(
							if(new_volume < 1e-14) {
								break;
							}
							
							if(new_volume < original_volume) {
								lower_limit_pos.set_pos(node[node_to_delete[0]].x, node[node_to_delete[0]].y, node[node_to_delete[0]].z);
							} else {
								upper_limit_pos.set_pos(node[node_to_delete[0]].x, node[node_to_delete[0]].y, node[node_to_delete[0]].z);
							}
							tolerance = fabs((new_volume - original_volume) / original_volume);
						}
						
						// For weird surfaces, volume conservation can require a big spike to form. Don't let it!
						for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
							for(int i = 0; i < 3; ++i) {
								if((*face_iterator)->n[i] != node_to_delete[0]) {
									nodal_distance_2 = node[node_to_delete[0]] - node[(*face_iterator)->n[i]];
									if(nodal_distance_2.get_magnitude() / nodal_distance.get_magnitude() > 3) {
										node[node_to_delete[0]].set_pos(midpoint.x, midpoint.y, midpoint.z);
										face_iterator = neighbourhood.end();
										face_iterator--;
										breaks++;
										break;
									} 
								}
							}	
						}
					}
					
				} else if(num_shared_faces == 4) {

					// Firstly check for a tetrahedral collapse (leaves two triangles sticking out of an otherwise consistent surface)
					// Delete the two evil triangles
					twod_section_exists = 0;
					for(i = 0; i < num_shared_faces; ++i) {
						for(j = i + 1; j < num_shared_faces; ++j) {
							if((*face_to_delete[i])->remaining_node == (*face_to_delete[j])->remaining_node) {
								twod_section_exists = 1;
							//	printf("Face %d : %d %d %d\n", i, (*face_to_delete[i])->n[0], (*face_to_delete[i])->n[1], (*face_to_delete[i])->n[2]);
							//	printf("Face %d : %d %d %d\n", j, (*face_to_delete[j])->n[0], (*face_to_delete[j])->n[1], (*face_to_delete[j])->n[2]);
								face_sharing_remaining_node[0] = i;
								face_sharing_remaining_node[1] = j;
								remaining_node = (*face_to_delete[i])->remaining_node;
								
								// Exit these loops
								j = num_shared_faces;
								i = num_shared_faces;
								break;
							}
						}
					}
					
					if(twod_section_exists == 1) {
								
						// Erase the triangles
						face.erase(face_to_delete[face_sharing_remaining_node[0]]);
						face.erase(face_to_delete[face_sharing_remaining_node[1]]);
						num_faces -= 2;

						// Check that this node isn't contained anywhere else
						// If it is, we have created two separate surfaces, but they're at least consistent!
						do_not_delete = 0;
						for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
							for(k = 0; k < 3; ++k) {
								if((*face_iterator)->n[k] == remaining_node) {
									do_not_delete = 1;
					
									// Exit both loops
									face_iterator = face.end();
									face_iterator--;
									break;
								}				
							}
						}

						// If it isn't, delete it
						if(do_not_delete == 0) {
							node[remaining_node].set_pos(INFINITY, INFINITY, INFINITY);
							num_nodes--;
						}
						
					} else {

						// If not, then we have a very complicated structure
						// 4-way coarsening may work. If there is a seg fault, or a bug, it is likely here
						// Volume conservation may be possible via more accurate definition of the surface structure, 
						// but currently not implemented
						
						// Move lowest node to midpoint
						midpoint.set_midpoint(node[node_to_delete[0]], node[node_to_delete[1]]);
						node[node_to_delete[0]].set_pos(midpoint.x, midpoint.y, midpoint.z);
	
						// Move highest node to infinity
						node[node_to_delete[1]].set_pos(INFINITY, INFINITY, INFINITY);
						num_nodes--;
	
						// Reconnect all neighbourhood faces to other node and recalculate shortest edge
						for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
						
							// Node indices
							for(i = 0; i < 3 ; ++i) {
								if((*face_iterator)->n[i] == node_to_delete[1]) {
									(*face_iterator)->n[i] = node_to_delete[0];
								}
							}	
						
							// Edge lengths
							for(i = 0; i < 3; ++i) {
								for(j = 0; j < 2; ++j) {
									if((*face_iterator)->edge[i].n_index[j] == node_to_delete[1]) {
										(*face_iterator)->edge[i].n_index[j] = node_to_delete[0];
										break;	
									}
								}
								(*face_iterator)->edge[i].calc_length(node);
							}
						
							// Shortest length
							(*face_iterator)->calc_shortest_edge();
						}
	
						// Delete faces
						face.erase(face_to_delete[3]);
						face.erase(face_to_delete[2]);
						face.erase(face_to_delete[1]);
						face.erase(face_to_delete[0]);
						num_faces = num_faces - 4;
				
						printf("4-way coarsening attempted. Good luck!\n");
					}
				} else if (num_shared_faces == 1) {

					// In a consistent manifold, a face with 1 edge cannot exist. So, get rid of it! Possibly get rid of some nodes too
					printf("Deleting a single face (which shouldn't even exist anyway!)\n");
					int deletenode[3] = {1,1,1};
					for(i = 0; i < 3; ++i) {
						for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
							for(j = 0; j < 3; ++j) {
								if((*face_to_delete[0])->n[i] == (*face_iterator)->n[j]) {
									deletenode[i] = 0;
								}
							}
						}
					}
					printf("Deleting nodes");
					for(i = 0; i < 3; ++i) {
						if(deletenode[i] == 1) {
							printf("%d ", (*face_to_delete[0])->n[i]);
							node[(*face_to_delete[0])->n[i]].set_pos(INFINITY, INFINITY, INFINITY);
							num_nodes--;
						}
					}
					printf("\n");
					face.erase(face_to_delete[0]);
					num_faces--;
				} else {
					printf("Something has gone wrong. %d faces are sharing an edge\n", num_shared_faces);
					printf("A face:%d %d %d\n", (*face_to_delete[0])->n[0], (*face_to_delete[0])->n[1], (*face_to_delete[0])->n[2]);
					break;
				}
				
				// Check neighbourhood consistency
				for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
					
					// Checking that faces haven't collapsed
					for(i = 0; i < 3; ++i) {
						for(j = i + 1; j < 3; ++j) {
							if((*face_iterator)->n[i] == (*face_iterator)->n[j]) {
								printf("In neighbourhood checker\n");
								exit(0);
								face.erase(face_iterator);
								num_faces--;
							}
						}
					}
				}
			}
			
			// Remap all nodes
			node_map = new int[num_nodes_initial];
			int new_node_index = 0;
			for(i = 0; i < num_nodes_initial; ++i) {
				if(node[i].x != INFINITY) {
					node_map[i] = new_node_index;
					node[new_node_index] = node[i];
					new_node_index++;
				} else {
					node_map[i] = -1;
				}
			}
			
			// Remap nodes in faces
			for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
				for(i = 0; i < 3; ++i) {
					(*face_iterator)->n[i] = node_map[(*face_iterator)->n[i]];
				}
			}
			return 0;
		}
		
		int clean() {
			printf("Cleaning...\n");
			list<Face*>::iterator face_iterator, neighbour_iterator, temp_iterator;
			//list<Face*> to_erase;
			list<int> to_erase;
			list<int>::iterator it, it2;

			int index[2];
			// Clear all non-manifold faces
			index[0] = -1;
			for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
				index[0]++;
				index[1] = index[0];
				temp_iterator = face_iterator;
				advance(temp_iterator, 1);
				for(neighbour_iterator = temp_iterator; neighbour_iterator != face.end(); ++neighbour_iterator) {
					index[1]++;
					// Check if same face
					if((*neighbour_iterator)->same_as(*face_iterator)) {
						//to_erase.push_back(*neighbour_iterator);
						//to_erase.push_back(*face_iterator);
						to_erase.push_back(index[0]);
						to_erase.push_back(index[1]);
					}
				}
				
			}
			to_erase.sort();
			int count = 0;
			for(it = to_erase.begin(); it != to_erase.end(); ++it) {
				face_iterator = face.begin();
				advance(face_iterator, (*it) - count);
				face.erase(face_iterator);
				count++;
			}

			return 0;
		}

		int write_to_file(char *surf_fname) {

			int i;
			list<Face*>::iterator face_iterator;	
			string fname(surf_fname), ext;
			ext = GetFileExtension(fname);

			// Open file
			FILE *surf_out;
			surf_out = fopen(surf_fname, "w");
			if(ext.compare("surf") == 0) {
				fprintf(surf_out, "surfacemesh\n");
				fprintf(surf_out, "%d\n", num_nodes);
				for(i = 0; i < num_nodes; ++i) {
					fprintf(surf_out, "%6.6lf %6.6lf %6.6lf\n", node[i].x, node[i].y, node[i].z);
				}
				fprintf(surf_out, "%d\n", num_faces);
				for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
					fprintf(surf_out, "%d %d %d\n", (*face_iterator)->n[0] + 1, (*face_iterator)->n[1] + 1, (*face_iterator)->n[2] + 1);
				}
			} else if (ext.compare("obj") == 0) {
				for(i = 0; i < num_nodes; ++i) {
					fprintf(surf_out, "v %6.6lf %6.6lf %6.6lf\n", node[i].x, node[i].y, node[i].z);
				}
				for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
					fprintf(surf_out, "f %d %d %d\n", (*face_iterator)->n[0] + 1, (*face_iterator)->n[1] + 1, (*face_iterator)->n[2] + 1);
				}
			} else if (ext.compare("stl") == 0) {
				
				// Header
				fprintf(surf_out, "solid ffeatoolsSTL\n");
				for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {

					// Calculate normal
					vector3 e0, e1, n;
					e0 = node[(*face_iterator)->n[1]] - node[(*face_iterator)->n[0]];
					e1 = node[(*face_iterator)->n[2]] - node[(*face_iterator)->n[0]];

					n = vector3::cross_product(e0, e1);
					n.normalise();

					// Now write everything
					fprintf(surf_out, "facet normal %3.2f %3.2f %3.2f\n", n.x, n.y, n.z);
					fprintf(surf_out, "outer loop\n");
					for(i = 0; i < 3; ++i) {
						fprintf(surf_out, "vertex %10.6f %10.6f %10.6f\n", node[(*face_iterator)->n[i]].x, node[(*face_iterator)->n[i]].y, node[(*face_iterator)->n[i]].z);
					}
					fprintf(surf_out, "endloop\n");
					fprintf(surf_out, "endfacet\n");
				}
				fprintf(surf_out, "endsolid ffeatoolsSTL\n");
			}
			fclose(surf_out);
			return 0;
		}
		
		void check_edges() {
			int i, j, index0, index1;
			list<Face*>::iterator face_iterator, face_iterator2;
			for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator){
				for(i = 0; i < 3; ++i) {
					index0 = (*face_iterator)->edge[i].n_index[0];
					index1 = (*face_iterator)->edge[i].n_index[1];
					face_iterator2 = face_iterator;
					face_iterator2++;
					while(face_iterator2 != face.end()) {
						for(j = 0; j < 3; ++j) {
							if((*face_iterator2)->edge[j].n_index[0] == index0 && (*face_iterator2)->edge[j].n_index[1] == index1) {
								printf("Edge %d %d has length %lf\n", (*face_iterator)->edge[i].n_index[0], (*face_iterator)->edge[i].n_index[1], (*face_iterator)->edge[i].length);
								printf("Edge %d %d has length %lf\n\n", (*face_iterator2)->edge[j].n_index[0], (*face_iterator2)->edge[j].n_index[1], (*face_iterator2)->edge[j].length);
							} else if ((*face_iterator2)->edge[j].n_index[0] == index1 && (*face_iterator2)->edge[j].n_index[1] == index0) {
								printf("Edge %d %d has length %lf\n", (*face_iterator)->edge[i].n_index[0], (*face_iterator)->edge[i].n_index[1], (*face_iterator)->edge[i].length);
								printf("Edge %d %d has length %lf\n\n", (*face_iterator2)->edge[j].n_index[0], (*face_iterator2)->edge[j].n_index[1], (*face_iterator2)->edge[j].length);
							}
						}
						face_iterator2++;
					}
				}
			}
		}
		
	private:
		
		// Surface components
		int num_nodes;
		int num_faces;
		int num_nodes_initial;
		int num_faces_initial;
		int *node_map;
		vector3 *node; 
		list<Face*> face, neighbourhood, local_volume_surf;
		Volume local_volume;
};


// Main program begins
int main(int argc, char **argv) {
   int ignore; 

	if(argc != 6 && argc != 12) {
		cout << "Input Error" << endl << "USAGE: " << argv[0] << " [Input .surf fname] [Output .surf fname] [Coarseness level] [Volume conserve (y/n)] [Find smallest edge? (y/n)] OPTIONAL[Coarsening range (xmin, xmax, ymin,ymax, zmin,zmax)]" << endl;
		return -1;
	}
	
	char achar[1];
	if(strcmp(argv[5], "N") == 0 || strcmp(argv[5], "n") == 0) {
		printf("\nYou have selected not to find the smallest edge at each pass.\n");
		printf("Although slightly quicker, this can result in a 'blockier' final mesh.\n");
		if(strcmp(argv[5], "Y") == 0 || strcmp(argv[5], "y") == 0) {
			printf("This method can also result in instabilities when conserving volume\n");
		}
		printf("Would you like to continue (y/n)?:");
		ignore = fscanf(stdin, "%s", achar);
		if(strcmp(achar, "n") == 0 || strcmp(achar, "N") == 0) {
			printf("Bye!");
			return 0;
		}
	}

	// Set limits
	double limits[6] = {-1 * INFINITY, INFINITY, -1 * INFINITY, INFINITY, -1 * INFINITY, INFINITY};
	if(argc == 12) {
		for(int i = 0; i < 6; ++i) {
			limits[i] = atof(argv[i + 6]);
		}
	}

	// Create a surface
	Surface *surf = new Surface();
	if(surf->init(argv[1]) != 0) {
		cout << "Error initialising surface\n.";
		cout << flush;
		return -1;
	};

	// Coarsen surface
	surf->coarsen(atof(argv[3]), argv[4], argv[5], limits);
	
	// Get rid of surface errors
	//surf->clean();

	// Output surface filebuf
	surf->write_to_file(argv[2]);
	
	// Delete surface
	delete surf;
	
	return 0;
}
