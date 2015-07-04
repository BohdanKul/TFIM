#include <stdlib.h>     // std::exit
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector
#include <array>
#include <cmath>
using namespace std;

//*****************************************************************************
// Determine the smallest constant necessarily to 
// offset the most negative diagonal weight to zero.  
//*****************************************************************************
float getDiagonalOffset(const float& J12, const float& h1, const float& h2){
    float offset = 0;
    float maxh   = max(abs(h1),abs(h2));
    float minh   = min(abs(h1),abs(h2));
    // For a ferromagnetic bond 
    if (J12 < 0)
        if (h1*h2>0) offset = maxh + abs(J12) + minh;   // alligned fields
        else         offset = maxh + abs(J12  - minh);  // anti-alligned
    // For an antiferromagnetic bond
    else
        if (h1*h2<0) offset = maxh + abs(J12) + minh;   // anti-alligned fields
        else         offset = maxh + abs(J12  - minh);  // alligned

    return offset;
}

//*****************************************************************************
// Compute the weight of a diagonal vertex 
//*****************************************************************************
float getDiagonalWeight(const float& J12, const float& h1, const float& h2, int Vertexid){
    // States of vertex's bottom legs
    int l0 = getLeg(Vertexid, 0); 
    int l1 = getLeg(Vertexid, 1);

    return J12*l0*l1 + h1*l0 + h2*l1;
}

//*****************************************************************************
// Determine vertex's type:
//     0 - diagonal
//     1 - illegal
//     2 - off-diagonal
//*****************************************************************************
int getVertexType(int Vertexid){
    // States of vertex's legs
    int l0; int l1; int l2; int l3;
    getAllLegs(Vertexid, l0, l1, l2, l3); 
    
    if       ((l0 == l3) and (l1 == l2)) return 0;
    else if  ((l0 != l3) and (l1 != l2)) return 1;
    else                                 return 2;
}

//*****************************************************************************
// Return the state of a given leg on a vertex 
//*****************************************************************************
int getLeg(int Vertexid, int leg){
    return  2*((Vertexid >> (3-leg))&1) - 1;
}

//*****************************************************************************
// Return the state of a all legs on a vertex 
//*****************************************************************************
void getAllLegs(int Vertexid, int& l0, int& l1, int& l2, int& l3){
    l0 = getLeg(Vertexid,0);
    l1 = getLeg(Vertexid,1);
    l2 = getLeg(Vertexid,2);
    l3 = getLeg(Vertexid,3);
}

//*****************************************************************************
// Determine whether a particular off-diagonal vertex has a given leg flipped
//*****************************************************************************
int isFlippedLeg(int Vertexid, int leg){
    // States of vertex's legs
    int l0; int l1; int l2; int l3;
    getAllLegs(Vertexid, l0, l1, l2, l3); 
    
    if (leg == 0) return int(l0 != l3);
    if (leg == 1) return int(l1 != l2); 
    
   cout << "Error: isFlipped() encountered a diagonal vertex" << endl;
   exit(0);    
}

//*****************************************************************************
// Return which one of an off-diagonal vertex' legs is flipped 
//*****************************************************************************
int getFlippedLeg(int Vertexid){
    // States of vertex's legs
    int l0; int l1; int l2; int l3;
    getAllLegs(Vertexid, l0, l1, l2, l3); 
    
    if (l0 != l3) return 0;
    if (l1 != l2) return 1;
    
   cout << "Error: getFlippedLeg() encountered a diagonal vertex" << endl;
   exit(0);    
}


//*****************************************************************************
// Get weights of all possible vertices 
//*****************************************************************************
void getVertices(const float& _J12, const float& _h1, const float& _h2, const float& _d1, const float _d2 , float& diagshift, array<float,16>& VWeights, array<float,4>& diagVWeights){
    diagshift = getDiagonalOffset(_J12, _h1, _h2);
    //array<bool, 16> VConform;
    int vtype;
    int diagi = 0;
    for (auto i = 0; i != 16; i++){
        vtype = getVertexType(i);
        if (vtype == 0){ 
            VWeights[i] = getDiagonalWeight(_J12, _h1, _h2, i) + diagshift;
            diagVWeights[diagi] = VWeights[i];
            diagi += 1;
        }
        if (vtype == 1) VWeights[i] = 0.0;
        if (vtype == 2) VWeights[i] = _d1*isFlippedLeg(i,0) + _d2*isFlippedLeg(i,1); 

        //if (VWeights == 0.0) VConform[i] = false;
        //else                 VConform[i] = true; 
    }
}


