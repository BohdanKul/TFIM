// lower_bound/upper_bound example
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector
#include <array>
#include <cmath>
#include "helper.cpp"   // FlipBit() and sort_indexes()

using namespace std;


//*****************************************************************************
// Compute cumulative probability of exiting on one of 8 possible legs 
// on leg-switch move during the off-diagonal update 
//*****************************************************************************
void getSwitchLegProb(const int& enleg, const int& vtype, vector<int>& Ws, array<float,8> prob){
    
    float totalW = 0;  // sum of all vertex weigths 

    // Construct all possible resulting vertices 
    
    // Start with vertices produced by flipping only the entrance leg
    int ntype;
    ntype = FlipBit(ntype, enleg);  // flip the entrance leg 
    if (VConform[ntype])
        for (auto ileg = 0; ileg != 4; ileg++){
            Ws[4+ileg] = VWeights[ntype];
            totalW    += VWeights[ntype];
        }

    // Add vertices produced by flipping the entraince and exit legs
    for (auto ileg = 0; ileg != 4; ileg++){
        ntype = FlipBit(vtype, ileg); 
        Ws[ileg] = VWeights[ntype];
        totalW  += VWeights[ntype];
    }


    // Convert weights and legs to the format required by the directed loop solutions
    array<int,   8> reIndex = sort_indexes(Ws); // reindex legs such based on their weights
    array<float, 8> SortedWs; // sorted weights
    int enLegSort = -1;       // leg index in the sorted array that corresponds the entrance leg
    for (auto leg = 0; leg != 8; leg++){
        if (reIndex[leg] == enLeg) enLegSort = leg;
        SortedWs[leg] = Ws[reIndex[leg]];
    }
    
   
    for (auto i=0; i!=8; i++) prob[i] = 0;

    // If possible, adopt the bounce-free solution B 
    if !(totalW < 2.0*SortedWs[0]) prob = noBounceSolutionB(SortedWs, enLegSort);
    else                           prob = BounceSolution(SortedWs, enLegSort); 
    
    // Convert probability density to a cumulative distribution 
    for (auto i=1; i!=8; i++)
        prob[i] = prob[i] + prob[i-1];
}

//*****************************************************************************
// Bounce solution based on bounce events only from the highest weight vertex
// Ws    - sorted array of weights from the heighest one to the lowest one.
//         Its indices correspond to exit legs.
// enLeg - entrance leg
//*****************************************************************************
array<float, 8> BounceSolution(const array<int,8>& Ws, const int& enLeg){
    // Probability of exiting on a particular leg
    array<float,8> prob = {0,0,0,0,0,0,0,0};

    // If the entrance leg corresponds to the highest weight vertex
    // there is a probability to switching to any other possible vertex
    if (enLeg == 0){
        float total = 0;
        for (auto i=1; i!=8; i++){
            prob[i] = Ws[i]/Ws[0];
            total += prob[i]
        }
        prob[0] = 1.0 - total;  // bounce probability
    }
    // Otherwise, switch to the highest weight leg with probability 1
    else prob[0] = 1.0;
     
    return prob;
}


//*****************************************************************************
// Bounce solution B based on Syljuasen's 
// "Direct Loop Updates for Quantum Lattice Models" paper
// Ws    - sorted array of weights from the heighest one to the lowest one.
//         Its indices correspond to exit legs.
// enLeg - entrance leg
//*****************************************************************************
array<float, 8> noBounceSolutionB(const array<int,8>& Ws, const int& enLeg){
    // Probability of exiting on a particular leg
    array<float,8> prob = {0,0,0,0,0,0,0,0};

    // If the entrance leg corresponds to the highest weight vertex
    if (enLeg == 0){
            prob[1] = (Ws[0] + Ws[1] - Ws[2] - Ws[3])/2.0/Ws[enLeg];
            prob[2] = (Ws[0] - Ws[1] + Ws[2] - Ws[3])/2.0/Ws[enLeg];
            for (auto leg = 3; leg!=7; leg++)
                prob[leg] = (Ws[leg] - Ws[leg+1])/2.0/Ws[enLeg];
            prob[[7] = Ws[7]/2.0/W[enLeg];
    }
    // If it is the second highest 
    else if (enLeg == 1){
             prob[0] = ( Ws[0] + Ws[1] - Ws[2] - Ws[3])/2.0/Ws[enLeg];
             prob[2] = (-Ws[0] + Ws[1] + Ws[2] + Ws[3])/2.0/Ws[enLeg];
    }
    // If it is the third highest 
    else if (enLeg == 2){
             prob[0] = ( Ws[0] - Ws[1] + Ws[2] - Ws[3])/2.0/Ws[enLeg];
             prob[1] = (-Ws[0] + Ws[1] + Ws[2] + Ws[3])/2.0/Ws[enLeg];
    }
    // If it is neither of the special previous cases
    else{
             // Probability to jump to the leg with the highest weight
             if (enLeg == 7) prob[0] = 0.5;
             else            prob[0] = 0.5 - Ws[enLeg+1]/2.0/Ws[enLeg];
             
             // Probability to jump to the leg with the next highest weight
             prob[enLeg-1] = 0.5;

             // Probability to jump to the leg with the next smallest weight
             if (enLeg !== 7) prob[enLeg+1] = Ws[enLeg+1]/2.0/Ws[enLeg];
    }

    return prob;
}

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
