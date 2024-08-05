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

#include "mat_vec_fns_II.h"
#include <iostream>
using std::cout;
using std::endl;

///////////////// SECTION 0 ////////////////////
////////  Constants and scalar functions ///////
////////////////////////////////////////////////
//
///  check whether two scalars have the same sign 
template <class t_scalar> bool sameSign(t_scalar a, t_scalar b){
 // return a <= ffea_const::zero == b <= ffea_const::zero;
 if (a*b >= 0) return true;
 return false;
}

/** Given 3 integers n0, n1, n2
 *    return the index missing in the list [0,1,2,3] 
 */
int getMissingNode(int n0, int n1, int n2) {

  if ((n0 == 0) || (n1 == 0) || (n2 ==0)) {
    if ((n0 == 1) || (n1 == 1) || (n2 ==1)) {
      if ((n0 == 2) || (n1 == 2) || (n2 ==2)) {
        return 3;
      } else {
        return 2;
      }
    } else {
      return 1;
    }
  } else {
    return 0;
  }

}

/** Given 1 integers iN,
 *    return the index missing indices of the list [0,1,2,3] 
 */
void getRestOfNodes(int iN, int &iO0, int &iO1, int &iO2) {

  if (iN == 0) {
    iO0 = 1; 
    iO1 = 2; 
    iO2 = 3; 
  } else if (iN == 1) {
    iO0 = 0; 
    iO1 = 2; 
    iO2 = 3; 
  } else if (iN == 2) {
    iO0 = 0; 
    iO1 = 1; 
    iO2 = 3; 
  } else {
    iO0 = 0; 
    iO1 = 1; 
    iO2 = 2; 
  }
}

/** Given 1 integers iN,
 *    return the index missing indices of the list [0,1,2,3] 
 */
void getMissingPair(int in0, int in1, int &on0, int &on1) {

  for (int i=0; i<4; i++) {
    if ((i != in0) && (i != in1)) {
      on0 = i;
      on1 = getMissingNode(in0, in1, on0);
      break; 
    }
  }
}


////////////////// SECTION 1 ///////////////////////
///  Basic operations for arr3, i. e., scalar v[3]// 
////////////////////////////////////////////////////
/** Add vectors vecA and vecB into res. */
/** res can also be vecA or vecB */
template <class t_scalar, class brr3> void arr3arr3Add(arr3_view<t_scalar,brr3> vecA, arr3_view<t_scalar,brr3> vecB, arr3_view<t_scalar,brr3> res){

   for (int i=0; i<3; i++) {
     res[i] = vecA[i] + vecB[i];
    }

}

/** res = vecA - vecB */
/** res can be either vecA or vecB */
// void arr3arr3Substract(arr3 &vecA, arr3 &vecB, arr3 &res){
// template <class brr3> void arr3arr3Substract(brr3 &vecA, brr3 &vecB, brr3 &res){
template <class t_scalar, class brr3> void arr3arr3Substract(arr3_view<t_scalar,brr3> vecA, arr3_view<t_scalar,brr3> vecB, arr3_view<t_scalar,brr3> res){

   for (int i=0; i<3; i++) {
     res[i] = vecA[i] - vecB[i];
    }

}
/** w = u x v */
/**  (w != u) && (w != v) */
// void arr3arr3VectorProduct(arr3 (&u), arr3 (&v), arr3 (&w)){
template <class t_scalar, class brr3> void arr3arr3VectorProduct(arr3_view<t_scalar,brr3> u, arr3_view<t_scalar,brr3> v, arr3_view<t_scalar,brr3> w){

    w[0] = u[1]*v[2] - v[1]*u[2];
    w[1] = -u[0]*v[2] + v[0]*u[2];
    w[2] = u[0]*v[1] - v[0]*u[1];


}

/** return the dot product for arrays vecA and vecB */
// scalar arr3arr3DotProduct(arr3 &vecA, arr3 &vecB) {
template <class t_scalar, class brr3> t_scalar arr3arr3DotProduct(arr3_view<t_scalar,brr3> vecA, arr3_view<t_scalar,brr3> vecB) {

  t_scalar result = 0.0;
  for (int i=0; i<3; i++) {
     // cout << "vecA[" << i << "]: " << vecA[i] << " vecB[" << i << "]: " << vecB[i] << endl; 
     result += vecA[i] * vecB[i];
  }
  return result;
}

/** Normalise vector arr3 e */
template <class t_scalar, class brr3> void arr3Normalise(arr3_view<t_scalar,brr3> e){

   t_scalar norm = 0.0;
   for (int i=0; i<3; i++) {
     norm += e[i]*e[i];
   }
   norm = sqrt(norm);
   for (int i=0; i<3; i++) {
     e[i] /= norm;
   }

}

/** get the normalised vector of arr3 e into arr3 n */
// template <class t_scalar, class brr3> void arr3Normalise2(brr3 &e, brr3 &n){
template <class t_scalar, class brr3> void arr3Normalise2(arr3_view<t_scalar,brr3> e, arr3_view<t_scalar,brr3> n){

   t_scalar norm = 0.0;
   for (int i=0; i<3; i++) {
     norm += e[i]*e[i];
   }
   if (norm == 0.0) throw -1; 
   norm = sqrt(norm);
   for (int i=0; i<3; i++) {
     n[i] = e[i]/norm;
   }

}

/** resize vector u, given scalar f */
template <class t_scalar, class brr3> void arr3Resize(t_scalar f, arr3_view<t_scalar,brr3> u){

   for (int i=0; i<3; i++) {
      u[i] = f*u[i];
    }

}

/** resize vector u into vector v, given scalar f */
template <class t_scalar, class brr3> void arr3Resize2(t_scalar f, arr3_view<t_scalar,brr3> u, arr3_view<t_scalar,brr3> v){

   #pragma omp simd
   for (int i=0; i<3; i++) {
      v[i] = f*u[i];
    }

}

/** Given a scalar f, change v so that v += f*u */
template <class t_scalar, class brr3> void arr3Resize3(t_scalar f, arr3_view<t_scalar,brr3> u, arr3_view<t_scalar,brr3> v){

   #pragma omp simd
   for (int i=0; i<3; i++) {
      v[i] += f*u[i]; 
   }
}


/** cp arr3 u into arr3 v */ 
template <class t_scalar, class brr3> void arr3Store(arr3_view<t_scalar,brr3> u, arr3_view<t_scalar,brr3> v){

   #pragma omp simd
   for (int i=0; i<3; i++) {
     v[i] = u[i];
    }

}

/** return the distance from vecA to vecB */
template <class t_scalar, class brr3> t_scalar arr3arr3Distance(arr3_view<t_scalar,brr3> vecA, arr3_view<t_scalar,brr3> vecB){
 
    t_scalar d=0.0;
    for (int i=0; i<3; i++){ 
      d  += (vecA[i] - vecB[i])*(vecA[i] - vecB[i]); 
    } 
    return sqrt(d); 
}


/** Return the length of a vector v */
template <class t_scalar, class brr3> t_scalar mag(arr3_view<t_scalar,brr3> v) {

   t_scalar s=0.0;
   #pragma omp simd reduction(+:s)
   for (int i=0; i<3; i++) {
      s += v[i] * v[i];
   }
   return sqrt(s);

}


/** Return the squared length of a vector v */
template <class t_scalar, class brr3> t_scalar mag2(arr3_view<t_scalar,brr3> v) {

   t_scalar s=0.0;
   for (int i=0; i<3; i++) {
      s += v[i] * v[i];
   }
   return s;

}


template <class brr3> void arr3Initialise(brr3 &v){

    for (int i=0; i<3; i++) {
      v[i] = ffea_const::zero;
    }

}

template <class t_scalar, class brr3> t_scalar detByRows(arr3_view<t_scalar,brr3> a, arr3_view<t_scalar,brr3> b, arr3_view<t_scalar,brr3> c){

  t_scalar det = 0;
  det  = a[0] * (b[1] * c[2] - b[2] * c[1]);
  det += a[1] * (b[2] * c[0] - b[0] * c[2]);
  det += a[2] * (b[0] * c[1] - b[1] * c[0]);
  return det; 

}

template <class t_scalar, class brr3> t_scalar detByCols(arr3_view<t_scalar,brr3> a, arr3_view<t_scalar,brr3> b, arr3_view<t_scalar,brr3> c){

  t_scalar det = 0;
  det  = a[0] * (b[1] * c[2] - b[2] * c[1]);
  det += b[0] * (c[1] * a[2] - c[2] * a[1]);
  det += c[0] * (a[1] * b[2] - a[2] * b[1]);
  return det; 

}




////////////////////////////////////////////////
////////////// END OF SECTION 1 ////////////////
////////////////////////////////////////////////


///////////////// SECTION 2 ////////////////////
////// Geometric functions for arr3 types //////
////////////////////////////////////////////////
/** t = unit(vecA - vecB) */
/** t can be either vecA or vecB */
// void tangent(arr3 &vecA, arr3 &vecB, arr3 &t){
template <class t_scalar, class brr3>  void tangent(brr3 &vecA, brr3 &vecB, brr3 &t){

   t_scalar w=0;
   for (int i=0; i<3; i++) {
     t[i] = vecA[i] - vecB[i];
     w += t[i] * t[i];
    }
    w = sqrt(w);
    for (int i=0; i<3; i++) {
      t[i] /= w;
    }

}

/** w = unit(u x v) */
/**  (w != u) && (w != v) */
template <class t_scalar, class brr3> void getUnitNormal(brr3 &u, brr3 &v, brr3 &w){

   w[0] = u[1]*v[2] - v[1]*u[2];
   w[1] = -u[0]*v[2] + v[0]*u[2];
   w[2] = u[0]*v[1] - v[0]*u[1];
   t_scalar l;
   l = sqrt( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] );
   for (int i=0; i<3; i++) {
     w[i] /= l;
   }
}


/** calculate the normal vector n to the plane defined by the three points */
template <class t_scalar, class brr3> void getNormal(brr3 &v1, brr3 &v2, brr3 &v3, brr3 &n){

   brr3 pl1, pl2;
   arr3arr3Substract<t_scalar,brr3>(v2, v1, pl1);
   arr3arr3Substract<t_scalar,brr3>(v3, v1, pl2);
   getUnitNormal<t_scalar,brr3>(pl1, pl2, n);

}


/* Given the face formed by tetA[0]:tetA[1]:tetA[2] 
 * get n, the normal to a face pointing inwards.
 */
template <class t_scalar, class brr3> void getNormalInwards(brr3 (&tetA)[4], int n0, int n1, int n2, brr3 (&n)){

   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   brr3 pl1, pl2, aux;
   arr3arr3Substract<t_scalar,brr3>(tetA[n1], tetA[n0], pl1);
   arr3arr3Substract<t_scalar,brr3>(tetA[n2], tetA[n0], pl2);
   // n is a unit vector normal to the face:
   arr3arr3VectorProduct<t_scalar,brr3>(pl1, pl2, aux);
   arr3Normalise2<t_scalar,brr3>(aux, n);
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   int n3 = getMissingNode(n0, n1, n2);
   arr3arr3Add<t_scalar,brr3>(aux, tetA[n0], aux);
   t_scalar d = - arr3arr3DotProduct<t_scalar,brr3>(n, tetA[n0]);
   t_scalar t1 = arr3arr3DotProduct<t_scalar,brr3>(tetA[n3], n) + d;
   t_scalar t2 = arr3arr3DotProduct<t_scalar,brr3>(aux, n) + d;
   if (!sameSign(t1,t2)) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n);

}

/** Given the face formed by f0, f1, and f2,
 *     and knowing the remaining p3 for a tetrahedron,
 * get n, the normal to a face pointing inwards.
 */
template <class t_scalar, class brr3> void getNormalInwards(brr3 &f0, brr3 &f1, brr3 &f2, brr3 &p3, brr3 (&n)){

   // pl1 and pl2 are the unit vectors (tetA[n1] - tetA[n0]) and (tetA[n2] - tetA[n0]), 
   brr3 pl1, pl2, aux;
   arr3arr3Substract<t_scalar,brr3>(f1, f0, pl1);
   arr3arr3Substract<t_scalar,brr3>(f2, f0, pl2);
   // n is a unit vector normal to the face:
   arr3arr3VectorProduct<t_scalar,brr3>(pl1, pl2, aux);
   arr3Normalise2<t_scalar,brr3>(aux, n);
   // but it must be inwards, i. e., on the same side of the plane than n3. 
   arr3arr3Add<t_scalar,brr3>(aux, f0, aux);
   t_scalar d = - arr3arr3DotProduct<t_scalar,brr3>(n, f0);
   t_scalar t1 = arr3arr3DotProduct<t_scalar,brr3>(p3, n) + d;
   t_scalar t2 = arr3arr3DotProduct<t_scalar,brr3>(aux, n) + d;
   if (!sameSign(t1,t2)) arr3Resize<t_scalar,brr3>(ffea_const::mOne, n);

}


/** check if points vec and test are at the same side
 *  of the plane formed by p1, p2 and p3 
 */
template <class t_scalar, class brr3> bool sameSidePlane(brr3 &vec, brr3 &test, brr3 &p1, brr3 &p2, brr3 &p3){

   brr3 pl1, pl2, n;
   arr3arr3Substract<t_scalar,brr3>(p2, p1, pl1);
   arr3arr3Substract<t_scalar,brr3>(p3, p1, pl2);
   arr3arr3VectorProduct<t_scalar,brr3>(pl1, pl2, n);
   arr3Normalise<t_scalar,brr3>(n);
   t_scalar d = - arr3arr3DotProduct<t_scalar,brr3>(n, p1);
   t_scalar t1 = arr3arr3DotProduct<t_scalar,brr3>(vec, n) + d;
   t_scalar t2 = arr3arr3DotProduct<t_scalar,brr3>(test, n) + d;
   // if ((t1 * t2) >= 0 ) return true;
   // cout << "t1: " << t1 << " t2: " << t2 << endl; 
   if (sameSign(t1,t2)) return true;
   else return false;

}

/**  Given 4 co-planar points, check if ip and p1 lay on the same side 
 *     of the of the line formed by p2 and p3. 
 *   More specifically we check whether pl23 x pl21 and pl23 x pl2e 
 *     are parallel or antiparallel. 
 */
template <class t_scalar, class brr3> bool sameSideLine(brr3 &e, brr3 &p1, brr3 &p2, brr3 &p3) {

   // if (!samePlane(e, p1, p2, p3)) cout << "alarm" << endl;
   brr3 pl23, pl21, pl2e, v1, ve;
   arr3arr3Substract<t_scalar,brr3>(p3, p2, pl23);
   arr3arr3Substract<t_scalar,brr3>(p1, p2, pl21);
   arr3arr3Substract<t_scalar,brr3>(e, p2, pl2e);
   arr3arr3VectorProduct<t_scalar,brr3>(pl21, pl23, v1);
   arr3arr3VectorProduct<t_scalar,brr3>(pl2e, pl23, ve);
   if (arr3arr3DotProduct<t_scalar,brr3>(v1,ve) >= 0) return true;
   return false;

}

/**  check whether vector vec is in tetrahedron B. */
/**  more specifically, it will be there if 
 *     for each plane of the tetrahedron, 
 *     the point is on the same side as the remaining vertex */
template <class t_scalar,class brr3> bool nodeInTet(brr3 &vec, brr3 (tet)[4]){

   if (!sameSidePlane<t_scalar,brr3>(vec, tet[0], tet[1], tet[2], tet[3])) return false;
   if (!sameSidePlane<t_scalar,brr3>(vec, tet[1], tet[2], tet[3], tet[0])) return false;
   if (!sameSidePlane<t_scalar,brr3>(vec, tet[2], tet[3], tet[0], tet[1])) return false;
   if (!sameSidePlane<t_scalar,brr3>(vec, tet[3], tet[0], tet[1], tet[2])) return false;

   return true;

}

template <class t_scalar,class brr3> bool nodeInTet(brr3 &vec, brr3 &tet0, brr3 &tet1, brr3 &tet2, brr3 &tet3){

   if (!sameSidePlane<t_scalar,brr3>(vec, tet0, tet1, tet2, tet3)) return false;
   if (!sameSidePlane<t_scalar,brr3>(vec, tet1, tet2, tet3, tet0)) return false;
   if (!sameSidePlane<t_scalar,brr3>(vec, tet2, tet3, tet0, tet1)) return false;
   if (!sameSidePlane<t_scalar,brr3>(vec, tet3, tet0, tet1, tet2)) return false;

   return true;

}

/** find the intersection point of the line that passes through the points e1 and e2, 
 *   and the plane defined by points p1, p2 and p3.
 *  \warning {this function should be called ONLY in the case 
 *            that intersection is known to occur.} 
 */
template <class t_scalar, class brr3> void linePlaneIntersectionPoint(brr3 &ip, brr3 &e1, brr3 &e2, brr3 &p1, brr3 &p2, brr3 &p3) {

   // v is the vector of the line L(t) = e1 + v*t
   brr3 v;
   arr3arr3Substract<t_scalar,brr3>(e2, e1, v);

   // now we need the unit vector that defines the plane, pn:
   brr3 pl1, pl2, pn;
   arr3arr3Substract<t_scalar,brr3>(p2, p1, pl1);
   arr3arr3Substract<t_scalar,brr3>(p3, p2, pl2);
   getUnitNormal<t_scalar,brr3>(pl1, pl2, pn);

   // the plane is defined through: ax + by + cz + d = 0; 
   //   (a,b,c) = pn
   // so to find d we simply:
   t_scalar d = - arr3arr3DotProduct<t_scalar,brr3>(pn, p1);

   // now find t and the point:
   t_scalar t = - (arr3arr3DotProduct<t_scalar,brr3>(pn, e1) + d) / arr3arr3DotProduct<t_scalar,brr3>(pn, v);
   arr3Resize<t_scalar,brr3>(t, v);
   arr3arr3Add<t_scalar,brr3>(e1, v, ip);

}

/** Return true and 
 *         the intersection point of the line that passes through the points e1 and e2
 *           and the plane defined by points p1, p2, p3 if this intersection actually occurs,
 *         and false otherwise. 
 *         \warning p1, p2, and p3 are assumed to be non-colinear.
 */
template <class t_scalar, class brr3> bool safeLinePlaneIntersectionPoint(brr3 &ip, brr3 &e1, brr3 &e2, brr3 &p1, brr3 &p2, brr3 &p3) {
   // v is the vector of the line L(t) = e1 + v*t
   brr3 v;
   arr3arr3Substract<t_scalar,brr3>(e2, e1, v);

   // now we need the unit vector that defines the plane, pn:
   brr3 pl1, pl2, pn;
   arr3arr3Substract<t_scalar,brr3>(p2, p1, pl1);
   arr3arr3Substract<t_scalar,brr3>(p3, p2, pl2);
   getUnitNormal<t_scalar,brr3>(pl1, pl2, pn);

   // CHECK! 
   scalar pnv = arr3arr3DotProduct<t_scalar,brr3>(pn, v);
   if ( abs(pnv) < ffea_const::threeErr ) return false; // the line won't intersect the plane!

   // the plane is defined through: ax + by + cz + d = 0; 
   //   (a,b,c) = pn
   // so to find d we simply:
   t_scalar d = - arr3arr3DotProduct<t_scalar,brr3>(pn, p1);

   // now find t and the point:
   t_scalar t = - (arr3arr3DotProduct<t_scalar,brr3>(pn, e1) + d) / arr3arr3DotProduct<t_scalar,brr3>(pn, v);
   arr3Resize<t_scalar,brr3>(t, v);
   arr3arr3Add<t_scalar,brr3>(e1, v, ip);

   return true; 

}

/** Return true and 
 *         the intersection point ip of the line that passes through the points e1 and e2
 *           and face defined by points p1, p2, p3 if this intersection actually occurs,
 *         and false otherwise. 
 *         \warning p1, p2, and p3 are assumed to be non-colinear.
 */
template <class t_scalar,class brr3> bool lineFaceIntersectionPoint(brr3 (&ip), brr3 (&e1), brr3 (&e2), brr3 (&p1), brr3 (&p2), brr3 (&p3)){

  // look for the intersection point... if it exists 
  if ( !safeLinePlaneIntersectionPoint<t_scalar,brr3>(ip, e1, e2, p1, p2, p3)) return false;

  // and finally check whether this point ip belongs to the triangular face:
  if ( (isPointInFace<t_scalar,brr3>(ip, p1, p2, p3)) ) return true;

  return false;

}


/** Check whether point ip is  
  *    inside of the three half-planes formed by the triangle's edges p1, p2, p3.
  */
template <class t_scalar, class brr3> bool isPointInFace(brr3 &ip, brr3 &p1, brr3 &p2, brr3 &p3) {

  if (! sameSideLine<t_scalar,brr3>(ip, p1, p2, p3) ) return false;
  if (! sameSideLine<t_scalar,brr3>(ip, p3, p1, p2) ) return false;
  if (! sameSideLine<t_scalar,brr3>(ip, p2, p3, p1) ) return false;

  return true;
}


/** Check whether an edge and a plane intersect, 
 *    and return the intersection point ip and true if found, false otherwise.
 * more specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (tet[f2] - tet[f1]) and (tet[f3] - tet[f1]).
 *    - the intersection of a line is a point in the plane 
 */
// bool intersectionPoint(arr3 &(ip), arr3 (&e1), arr3 (&e2), arr3 (&tet)[4], int f1, int f2, int f3){
template <class t_scalar,class brr3> bool intersectionPoint(brr3 &(ip), brr3 (&e1), brr3 (&e2), brr3 (&tet)[4], int f1, int f2, int f3){

  // check it e1 and e2 are on the same side of the plane:
  if ( sameSidePlane<t_scalar,brr3>(e1, e2, tet[f1], tet[f2], tet[f3]) ) return false;

  // given that they are on different sides of the plane look for the intersection point.
  linePlaneIntersectionPoint<t_scalar,brr3>(ip, e1, e2, tet[f1], tet[f2], tet[f3]);

  // and finally check whether this point ip belongs to the triangular face:
  if ( (isPointInFace<t_scalar,brr3>(ip, tet[f1], tet[f2], tet[f3])) ) return true;


  return false;

}


/** Check whether an edge and a plane intersect, 
 *    and return the intersection point ip and true if found, false otherwise.
 * more specifically check that both:
 *    - both ends of the edge (e1 and e2) are on different sides
 *           of the plane defined by the vectors (f2 - f1) and (f3 - f1).
 *    - the intersection of a line is a point in the plane 
 */
template <class t_scalar,class brr3> bool intersectionPoint(brr3 &ip, brr3 &e1, brr3 &e2, brr3 &f1, brr3 &f2, brr3 &f3){

  // check it e1 and e2 are on the same side of the plane:
  if ( sameSidePlane<t_scalar,brr3>(e1, e2, f1, f2, f3) ) return false;

  // given that they are on different sides of the plane look for the intersection point.
  linePlaneIntersectionPoint<t_scalar,brr3>(ip, e1, e2, f1, f2, f3);

  // and finally check whether this point ip belongs to the triangular face:
  if ( (isPointInFace<t_scalar,brr3>(ip, f1, f2, f3)) ) return true;


  return false;

}


template <class t_scalar, class brr3> void intersectingPointToLine(arr3_view<t_scalar,brr3> p0, arr3_view<t_scalar,brr3> p1, arr3_view<t_scalar,brr3> p2p1, arr3_view<t_scalar,brr3> p3) {

   brr3 p0p1, cp3p0, v2, v1v2;
   t_scalar c;
   
   arr3arr3Substract<t_scalar,brr3>(p0, p1, p0p1);
   arr3arr3VectorProduct<t_scalar,brr3>(p0p1, p2p1, v2);
   // CHECK! if the cross product is null, then distance is zero, and p3 is p0 itself.
   if (mag<t_scalar,brr3>(v2) < ffea_const::threeErr) {
     arr3Store(p0, p3); 
     return; 
   } 
   // get cp3p0, the vector to (or from) to the line p2p1:
   arr3arr3VectorProduct<t_scalar,brr3>(p2p1, v2, cp3p0); 

   // now calculate c, the amount of cp3p0 from p0 to the intersection point p3.
   arr3arr3VectorProduct<t_scalar,brr3>(p2p1, cp3p0, v1v2);
   c = detByRows<t_scalar, brr3>(p0p1, p2p1, v1v2); 
   c /= ( v1v2[0]*v1v2[0] + v1v2[1]*v1v2[1] + v1v2[2]*v1v2[2]);

   // and finally calculate the intersection point: 
   arr3Resize<t_scalar,brr3>(c, cp3p0); 
   arr3arr3Add<t_scalar,brr3>(p0, cp3p0, p3); 

   /*
   // CHECKS START HERE! // 
   brr3 p3p0; 
   arr3arr3Substract<t_scalar,brr3>(p3, p0, p3p0);
   c = mag<t_scalar,brr3>(p3p0);

   // check the distance:
   t_scalar d = distanceFromPointToLine<t_scalar, brr3>(p0, p1, p2); 
   if (abs (d - abs(c)) > ffea_const::threeErr) {
     cout << "something is wrong: d = " << d << " differs from c = " << c << endl; 
   } 
   */

   return;

}


template <class t_scalar, class brr3> void intersectingPointToLine(vector3 &p0, arr3_view<t_scalar,brr3> p1, arr3_view<t_scalar,brr3> p2p1, arr3_view<t_scalar,brr3> p3) {

   brr3 p0p1, cp3p0, v2, v1v2;
   t_scalar c;
   
   // arr3arr3Substract<t_scalar,brr3>(p2, p1, p2p1);
   // arr3arr3Substract<t_scalar,brr3>(p0, p1, p0p1);
   p0p1[0] = p0.x - p1[0]; 
   p0p1[1] = p0.y - p1[1]; 
   p0p1[2] = p0.z - p1[2]; 
   arr3arr3VectorProduct<t_scalar,brr3>(p0p1, p2p1, v2);
   // CHECK! if the cross product is null, then distance is zero, and p3 is p0 itself.
   if (mag<t_scalar,brr3>(v2) < ffea_const::threeErr) {
     // arr3Store(p0, p3); 
     p3[0] = p0.x;
     p3[1] = p0.y;
     p3[2] = p0.z;
     return;
   } 
   // get cp3p0, the vector to (or from) to the line p2p1:
   arr3arr3VectorProduct<t_scalar,brr3>(p2p1, v2, cp3p0); 

   // now calculate c, the amount of cp3p0 from p0 to the intersection point p3.
   arr3arr3VectorProduct<t_scalar,brr3>(p2p1, cp3p0, v1v2);
   c = detByRows<t_scalar, brr3>(p0p1, p2p1, v1v2); 
   c /= ( v1v2[0]*v1v2[0] + v1v2[1]*v1v2[1] + v1v2[2]*v1v2[2]);

   // and finally calculate the intersection point: 
   arr3Resize<t_scalar,brr3>(c, cp3p0); 
   p3[0] = p0.x + cp3p0[0]; 
   p3[1] = p0.y + cp3p0[1]; 
   p3[2] = p0.z + cp3p0[2]; 

   /*
   // CHECKS START HERE! // 
   brr3 p3p0; 
   arr3arr3Substract<t_scalar,brr3>(p3, p0, p3p0);
   c = mag<t_scalar,brr3>(p3p0);

   // check the distance:
   t_scalar d = distanceFromPointToLine<t_scalar, brr3>(p0, p1, p2); 
   if (abs (d - c) > ffea_const::threeErr) {
     cout << "something is wrong: d = " << d << " differs from c = " << c << endl; 
   } 
   */

   return; 

}


template <class t_scalar, class brr3> t_scalar distanceFromPointToLine(arr3_view<t_scalar,brr3> p0, arr3_view<t_scalar,brr3> p1, arr3_view<t_scalar,brr3> p2) {

   brr3 p2p1, p0p1, p0p2, tmp1;
   arr3arr3Substract<t_scalar,brr3>(p2, p1, p2p1);
   arr3arr3Substract<t_scalar,brr3>(p0, p2, p0p2);
   arr3arr3Substract<t_scalar,brr3>(p0, p1, p0p1);
   arr3arr3VectorProduct<t_scalar,brr3>(p0p1, p0p2, tmp1);
   scalar d = mag<t_scalar,brr3>(tmp1) / mag<t_scalar,brr3>(p2p1);
   return d; 

}


/* template <class t_scalar, class brr3> t_scalar getTetrahedraVolume(arr3_view<t_scalar,brr3> p0, arr3_view<t_scalar,brr3> p1, arr3_view<t_scalar,brr3> p2, arr3_view<t_scalar,brr3> p3){

   scalar v = 0;   
   v += detByCols<scalar,arr3>(p1, p2, p3);
   v -= detByCols<scalar,arr3>(p0, p2, p3);
   v += detByCols<scalar,arr3>(p0, p1, p3);
   v -= detByCols<scalar,arr3>(p0, p1, p2);
   return v/6.; 
  
}*/ 


template <class t_scalar, class brr3> t_scalar getTetrahedraVolume(arr3_view<t_scalar,brr3> p0, arr3_view<t_scalar,brr3> p1, arr3_view<t_scalar,brr3> p2, arr3_view<t_scalar,brr3> p3){

  return  (p1[0] - p0[0])*( (p1[1] - p2[1])*(p2[2] - p3[2]) - (p2[1] - p3[1])*(p1[2] - p2[2]) ) +
          (p2[0] - p1[0])*( (p2[1] - p3[1])*(p0[2] - p1[2]) - (p0[1] - p1[1])*(p2[2] - p3[2]) ) +
          (p3[0] - p2[0])*( (p0[1] - p1[1])*(p1[2] - p2[2]) - (p1[1] - p2[1])*(p0[2] - p1[2]) );
  
}


template <class brr3> void getTetrahedraCM(brr3 &p1, brr3 &p2, brr3 &p3, brr3 &p4, brr3 &c){

   for (int i=0; i<3; i++){
      c[i] = 0.25 * (p1[i] + p2[i] + p3[i] + p4[i]);
   }

}

template <class t_scalar, class brr3, class brr4> void getLocalCoordinatesForLinTet(arr3_view<t_scalar,brr3> t0, arr3_view<t_scalar,brr3> t1, arr3_view<t_scalar,brr3> t2, arr3_view<t_scalar,brr3> t3, arr3_view<t_scalar,brr3> p, brr4 &phi){

   phi[0] = getTetrahedraVolume(p, t1, t2, t3);
   phi[1] = - getTetrahedraVolume(p, t0, t2, t3);
   phi[2] = getTetrahedraVolume(p, t0, t1, t3);
   phi[3] = - getTetrahedraVolume(p, t0, t1, t2);
   t_scalar v = getTetrahedraVolume(t0, t1, t2, t3);
   for (int i=0; i<4; i++) {
     phi[i] /= v; 
   } 

   /*  // CHECK!! for the time: 
   arr3 tmp, pc;
   arr3Initialise<arr3>(pc); // WT*?
   arr3Resize2<scalar,arr3>(phi[0], t0, tmp);
   arr3arr3Add<scalar,arr3>(pc, tmp, pc);
   arr3Resize2<scalar,arr3>(phi[1], t1, tmp);
   arr3arr3Add<scalar,arr3>(pc, tmp, pc);
   arr3Resize2<scalar,arr3>(phi[2], t2, tmp);
   arr3arr3Add<scalar,arr3>(pc, tmp, pc);
   arr3Resize2<scalar,arr3>(phi[3], t3, tmp);
   arr3arr3Add<scalar,arr3>(pc, tmp, pc);
   arr3arr3Substract<scalar,arr3>(p, pc, tmp);
   if (mag<scalar,arr3>(tmp) > 1e-6) {
       cout << "local coordinates were not correctly calculated!" << endl; 
       cout << "p: " << p[0] << ", " << p[1] << ", " << p[2] << endl;
       cout << "pc: " << pc[0] << ", " << pc[1] << ", " << pc[2] << endl;
       cout << "diff: " << mag<scalar,arr3>(tmp) << endl; 
   } */


   
}




////////////////////////////////////////////////
////////////// END OF SECTION 2 ////////////////
////////////////////////////////////////////////


///////////////// SECTION 3 ////////////////////
/// Transition functions from vector3 to arr3 // 
////////////////////////////////////////////////
template <class brr3> void vec3Vec3SubsToArr3(vector3 &u, vector3 &v, brr3 (&w)){

    w[0] = u.x - v.x;
    w[1] = u.y - v.y;
    w[2] = u.z - v.z;

}

void vec3Arr3SubsToArr3(vector3 &u, arr3 (&v), arr3 (&w)){

    w[0] = u.x - v[0];
    w[1] = u.y - v[1];
    w[2] = u.z - v[2];

}

void arr3Vec3SubsToArr3(arr3 (&u), vector3 &v, arr3 (&w)){

    w[0] = u[0] - v.x;
    w[1] = u[1] - v.y;
    w[2] = u[2] - v.z;

}

void vec3Arr3AddToArr3(vector3 &u, arr3 (&v), arr3 (&w)){

    w[0] = u.x + v[0];
    w[1] = u.y + v[1];
    w[2] = u.z + v[2];

}

void vec3ResizeToArr3(scalar f, vector3 &u, arr3 (&v)){

    v[0] = f*u.x;
    v[1] = f*u.y;
    v[2] = f*u.z;

}

scalar vec3Arr3DotProduct(vector3 &u, arr3 &v) {

     scalar s;
     s = u.x * v[0];
     s += u.y * v[1];
     s += u.z * v[2];
     return s;

}


////////////////////////////////////////////////
////////////// END OF SECTION 3 ////////////////
////////////////////////////////////////////////


///////////////// SECTION 4 ////////////////////
// // // // Instantiate templates // // // // //
////////////////////////////////////////////////
template bool sameSign<scalar>(scalar a, scalar b);

template void arr3arr3Add<scalar, arr3>(arr3_view<scalar,arr3> vecA, arr3_view<scalar,arr3> vecB, arr3_view<scalar,arr3> res);

template void arr3arr3Substract<scalar,arr3>(arr3_view<scalar,arr3> vecA, arr3_view<scalar,arr3> vecB, arr3_view<scalar,arr3> res);

template void arr3arr3VectorProduct<scalar,arr3>(arr3_view<scalar,arr3> u, arr3_view<scalar,arr3> v, arr3_view<scalar,arr3> w);

template scalar arr3arr3DotProduct<scalar,arr3>(arr3_view<scalar,arr3> vecA, arr3_view<scalar,arr3> vecB);

template void arr3Normalise<scalar,arr3>(arr3_view<scalar,arr3> e);

template void arr3Normalise2<scalar,arr3>(arr3_view<scalar,arr3> e, arr3_view<scalar,arr3> n);

template void arr3Resize<scalar,arr3>(scalar f, arr3_view<scalar,arr3> u);

template void arr3Resize2<scalar,arr3> (scalar f, arr3_view<scalar,arr3> u, arr3_view<scalar,arr3> v);
template void arr3Resize3<scalar,arr3> (scalar f, arr3_view<scalar,arr3> u, arr3_view<scalar,arr3> v);

template void arr3Store<scalar,arr3>(arr3_view<scalar,arr3> u, arr3_view<scalar,arr3> v);

template scalar arr3arr3Distance<scalar,arr3>(arr3_view<scalar,arr3> vecA, arr3_view<scalar,arr3> vecB); 

template scalar mag<scalar,arr3>(arr3_view<scalar,arr3> v);

template scalar mag2<scalar,arr3>(arr3_view<scalar,arr3> v);

template void arr3Initialise<arr3>(arr3 &v);

template scalar detByRows<scalar,arr3>(arr3_view<scalar,arr3> a, arr3_view<scalar,arr3> b, arr3_view<scalar,arr3> c);

template scalar detByCols<scalar,arr3>(arr3_view<scalar,arr3> a, arr3_view<scalar,arr3> b, arr3_view<scalar,arr3> c);


#ifndef USE_DOUBLE
template bool sameSign<geoscalar>(geoscalar a, geoscalar b);
template void arr3arr3Add<geoscalar, grr3>(arr3_view<geoscalar,grr3> vecA, arr3_view<geoscalar,grr3> vecB, arr3_view<geoscalar,grr3> res);
template void arr3arr3Substract<geoscalar,grr3>(arr3_view<geoscalar,grr3> vecA, arr3_view<geoscalar,grr3> vecB, arr3_view<geoscalar,grr3> res);
template void arr3arr3VectorProduct<geoscalar,grr3>(arr3_view<geoscalar,grr3> u, arr3_view<geoscalar,grr3> v, arr3_view<geoscalar,grr3> w);
template geoscalar arr3arr3DotProduct<geoscalar,grr3>(arr3_view<geoscalar,grr3> vecA, arr3_view<geoscalar,grr3> vecB);
template void arr3Normalise<geoscalar,grr3>(arr3_view<geoscalar,grr3> e);
template void arr3Normalise2<geoscalar,grr3>(arr3_view<geoscalar,grr3> e, arr3_view<geoscalar,grr3> n);
template void arr3Resize<geoscalar,grr3>(geoscalar f, arr3_view<geoscalar,grr3> u);
template void arr3Resize2<geoscalar,grr3> (geoscalar f, arr3_view<geoscalar,grr3> u, arr3_view<geoscalar,grr3> v);
template void arr3Resize3<geoscalar,grr3> (geoscalar f, arr3_view<geoscalar,grr3> u, arr3_view<geoscalar,grr3> v);
template void arr3Store<geoscalar,grr3>(arr3_view<geoscalar,grr3> u, arr3_view<geoscalar,grr3> v);
template geoscalar arr3arr3Distance<geoscalar,grr3>(arr3_view<geoscalar,grr3> vecA, arr3_view<geoscalar,grr3> vecB); 
template geoscalar mag<geoscalar,grr3>(arr3_view<geoscalar,grr3> v);
template geoscalar mag2<geoscalar,grr3>(arr3_view<geoscalar,grr3> v);
template void arr3Initialise<grr3>(grr3 &v);
template geoscalar detByRows<geoscalar,grr3>(arr3_view<geoscalar,grr3> a, arr3_view<geoscalar,grr3> b, arr3_view<geoscalar,grr3> c);
template geoscalar detByCols<geoscalar,grr3>(arr3_view<geoscalar,grr3> a, arr3_view<geoscalar,grr3> b, arr3_view<geoscalar,grr3> c);
#endif


template void tangent<scalar,arr3>(arr3 &vecA, arr3 &vecB, arr3 &t);

template void getUnitNormal<scalar,arr3>(arr3 &u, arr3 &v, arr3 &w);

template void getNormal<scalar,arr3>(arr3 &v1, arr3 &v2, arr3 &v3, arr3 &n);

template void getNormalInwards<scalar,arr3>(arr3 (&tetA)[4], int n0, int n1, int n2, arr3 (&n));
template void getNormalInwards<scalar,arr3>(arr3 &f0, arr3 &f1, arr3 &f2, arr3 &p3, arr3 (&n));

template bool sameSidePlane<scalar,arr3>(arr3 &vec, arr3 &test, arr3 &p1, arr3 &p2, arr3 &p3);

template bool sameSideLine<scalar,arr3>(arr3 &e, arr3 &p1, arr3 &p2, arr3 &p3);

template bool nodeInTet<scalar,arr3>(arr3 &vec, arr3 (tet)[4]);
template bool nodeInTet<scalar,grr3>(arr3 &vec, arr3 &tet0, arr3 &tet1, arr3 &tet2, arr3 &tet3);

template void linePlaneIntersectionPoint<scalar,arr3>(arr3 &ip, arr3 &e1, arr3 &e2, arr3 &p1, arr3 &p2, arr3 &p3);

template bool safeLinePlaneIntersectionPoint<scalar,arr3>(arr3 &ip, arr3 &e1, arr3 &e2, arr3 &p1, arr3 &p2, arr3 &p3);

template bool lineFaceIntersectionPoint<scalar,arr3>(arr3 (&ip), arr3 (&e1), arr3 (&e2), arr3 (&p1), arr3 (&p2), arr3 (&p3));

template bool isPointInFace<scalar,arr3>(arr3 &ip, arr3 &p1, arr3 &p2, arr3 &p3);

template bool intersectionPoint<scalar,arr3>(arr3 &(ip), arr3 (&e1), arr3 (&e2), arr3 (&tet)[4], int f1, int f2, int f3);
template bool intersectionPoint<scalar,arr3>(arr3 &ip, arr3 &e1, arr3 &e2, arr3 &f1, arr3 &f2, arr3 &f3);


// template void intersectingPointToLine<scalar, arr3>(arr3 &p0, arr3 &p1, arr3 &p2, arr3 &p3);
template void intersectingPointToLine<scalar, arr3>(arr3_view<scalar, arr3> p0, arr3_view<scalar,arr3> p1, arr3_view<scalar,arr3> p2p1, arr3_view<scalar,arr3> p3);

template void intersectingPointToLine<scalar, arr3>(vector3 &p0, arr3_view<scalar,arr3> p1, arr3_view<scalar,arr3> p2p1, arr3_view<scalar,arr3> p3);


template scalar distanceFromPointToLine<scalar,arr3>(arr3_view<scalar, arr3> p0, arr3_view<scalar, arr3> p1, arr3_view<scalar, arr3> p2);

template scalar getTetrahedraVolume<scalar,arr3>(arr3_view<scalar,arr3> p0, arr3_view<scalar,arr3> p1, arr3_view<scalar,arr3> p2, arr3_view<scalar,arr3> p3);

template void getTetrahedraCM<arr3>(arr3 &p1, arr3 &p2, arr3 &p3, arr3 &p4, arr3 &c);

template void getLocalCoordinatesForLinTet<scalar,arr3,arr4>(arr3_view<scalar,arr3> t0, arr3_view<scalar,arr3> t1, arr3_view<scalar,arr3> t2, arr3_view<scalar,arr3> t3, arr3_view<scalar,arr3> p, arr4 &phi);

   //////////////     
template void vec3Vec3SubsToArr3<arr3>(vector3 &u, vector3 &v, arr3 (&w));



#ifndef USE_DOUBLE
template void tangent<geoscalar,grr3>(grr3 &vecA, grr3 &vecB, grr3 &t);
template void getUnitNormal<geoscalar,grr3>(grr3 &u, grr3 &v, grr3 &w);
template void getNormal<geoscalar,grr3>(grr3 &v1, grr3 &v2, grr3 &v3, grr3 &n);
template void getNormalInwards<geoscalar,grr3>(grr3 (&tetA)[4], int n0, int n1, int n2, grr3 (&n));
template void getNormalInwards<geoscalar,grr3>(grr3 &f0, grr3 &f1, grr3 &f2, grr3 &p3, grr3 (&n));
template bool sameSidePlane<geoscalar,grr3>(grr3 &vec, grr3 &test, grr3 &p1, grr3 &p2, grr3 &p3);
template bool sameSideLine<geoscalar,grr3>(grr3 &e, grr3 &p1, grr3 &p2, grr3 &p3);
template bool nodeInTet<geoscalar,grr3>(grr3 &vec, grr3 (tet)[4]);
template bool nodeInTet<geoscalar,grr3>(grr3 &vec, grr3 &tet0, grr3 &tet1, grr3 &tet2, grr3 &tet3);
template void linePlaneIntersectionPoint<geoscalar,grr3>(grr3 &ip, grr3 &e1, grr3 &e2, grr3 &p1, grr3 &p2, grr3 &p3);
template bool safeLinePlaneIntersectionPoint<geoscalar,grr3>(grr3 &ip, grr3 &e1, grr3 &e2, grr3 &p1, grr3 &p2, grr3 &p3);
template bool lineFaceIntersectionPoint<geoscalar,grr3>(grr3 (&ip), grr3 (&e1), grr3 (&e2), grr3 (&p1), grr3 (&p2), grr3 (&p3));
template bool isPointInFace<geoscalar,grr3>(grr3 &ip, grr3 &p1, grr3 &p2, grr3 &p3);
template bool intersectionPoint<geoscalar,grr3>(grr3 &(ip), grr3 (&e1), grr3 (&e2), grr3 (&tet)[4], int f1, int f2, int f3);
template bool intersectionPoint<geoscalar,grr3>(grr3 &ip, grr3 &e1, grr3 &e2, grr3 &f1, grr3 &f2, grr3 &f3);
template void intersectingPointToLine<geoscalar, grr3>(arr3_view<geoscalar, grr3> p0, arr3_view<geoscalar,grr3> p1, arr3_view<geoscalar,grr3> p2p1, arr3_view<geoscalar,grr3> p3);
template void intersectingPointToLine<geoscalar, grr3>(vector3 &p0, arr3_view<geoscalar,grr3> p1, arr3_view<geoscalar,grr3> p2p1, arr3_view<geoscalar,grr3> p3);
template geoscalar distanceFromPointToLine<geoscalar,grr3>(arr3_view<geoscalar, grr3> p0, arr3_view<geoscalar, grr3> p1, arr3_view<geoscalar, grr3> p2);
template geoscalar getTetrahedraVolume<geoscalar,grr3>(arr3_view<geoscalar,grr3> p0, arr3_view<geoscalar,grr3> p1, arr3_view<geoscalar,grr3> p2, arr3_view<geoscalar,grr3> p3);
template void getTetrahedraCM<grr3>(grr3 &p1, grr3 &p2, grr3 &p3, grr3 &p4, grr3 &c);
template void getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(arr3_view<geoscalar,grr3> t0, arr3_view<geoscalar,grr3> t1, arr3_view<geoscalar,grr3> t2, arr3_view<geoscalar,grr3> t3, arr3_view<geoscalar,grr3> p, grr4 &phi);
template void vec3Vec3SubsToArr3<grr3>(vector3 &u, vector3 &v, grr3 (&w));
#endif 

