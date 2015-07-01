// lower_bound/upper_bound example
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector
#include <array>
#include <cmath>

using namespace std;

// Determine the smallest constant necessarily to 
// offset the most negative diagonal weight to zero.  
float DiagonalOffset(const float& J12, const float& h1, const float& h2){
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

// Compute the weight 
float DiagonalWeight(const float& J12, const float& h1, const float& h2, int Vertexid){
    // States of vertex's bottom legs
    int l0 = 2*((Vertexid >> 3)&1) - 1;
    int l1 = 2*((Vertexid >> 2)&1) - 1;

    return J12*l0*l1 + h1*l0 + h2*l1;
}

// Determine vertex's type:
//     0 - diagonal
//     1 - illegal
//     2 - off-diagonal
int VertexType(int Vertexid){
    // States of vertex's legs
    int l0 = 2*((Vertexid >> 3)&1) - 1;
    int l1 = 2*((Vertexid >> 2)&1) - 1;
    int l2 = 2*((Vertexid >> 1)&1) - 1;
    int l3 = 2*((Vertexid >> 0)&1) - 1;
    //cout << l0 << " " << l1 << " " << l2 << " " << l3 << endl;
    
    if       ((l0 == l3) and (l1 == l2)) return 0;
    else if  ((l0 != l3) and (l1 != l2)) return 1;
    else                                 return 2;
}

//int FlipVertex(int Vertexid, int el){
//    Vertexid

int FlipBit(const int& BitStr, const int& index){
    return BitStr ^ (1 << index);
}

//vector<vector<vector<vector<int>>>> 
//for (int exflip=0: {false,true}){
//    for (auto i=0; i!=2*16-1; i++){
//        
//int VertexMoves[16][;
//for (auto i=0; i!=16; i++){
//    for (int exflip=0: {false,true}){
//        
//int VertexId = 0;
//
//
//array<float,16> VWeights;
//int vtype;
//float diagshift = DiagonalOffset(J12, h1, h2);

//float J12 = 1.0; float h1 = -1.0; float h2 = 1.0; float delta = 1.0;
//for (auto i = 0; i != 16; i++){
//    vtype = VertexType(i);
//    if (vtype == 0) VWeights[i] = DiagonalWeight(J12, h1, h2, i) + diagshift;
//    if (vtype == 1) VWeights[i] = 0;
//    if (vtype == 2) VWeights[i] = delta; 
//}

// For a given vertex and leg, construct a vector of possible 
// vertices that can result in the off-diagonal update. 
//void getPossibleVertices(const int& l, const int& vtype, vector<int>& vertices){
//    int ntype = vtype;
//    vertices.push_back(ntype); // the original vertex is possible upon bounce
//    ntype = FlipBit(ntype,l);  // all other vertices have the entrance leg flipped 
//    vertices.push_back(ntype); // resulting vertex when no exit leg flip occurs
//
//    // Resulting vertices when the exit leg is flipped
//    int ttype;
//    for (auto ileg = 0; ileg != 4; ileg++){
//        if (ileg != l){
//            ttype = FlipBit(ntype, ileg); 
//            if (VWeights[ttype] != 0) vertices.push_back(ttype);
//        }
//    }
//}

//array
//float maxW = max_element(

int main () {
  
    array<array<array<int, 8>, 2>, 16> myarray;     
    cout << myarray[15][1][100] << endl;
    cout << myarray.at(2).at(1).at(7) << endl;
    int vt;
    for (int i=0; i!=16; i++){
        cout << " Vertex: " << i << " Legs: ";
        vt = VertexType(i);
        if (vt==0) cout << "   diagonal" << endl; 
        if (vt==1) cout << "   illegal" << endl; 
        if (vt==2) cout << "   off-diagonal" << endl; 
    }
    return 0;
}
