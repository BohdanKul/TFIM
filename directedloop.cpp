// lower_bound/upper_bound example
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector
#include <array>
#include <cmath>
#include "helper.cpp"   // flipBit() and sort_indexes()

using namespace std;
//*****************************************************************************
// Solution implementing the heat bath construction of the loop 
// Ws    - sorted array of weights from the heighest one to the lowest one.
//         Its indices correspond to exit legs.
// enLeg - entrance leg
//*****************************************************************************
array<float, 8> HeatBathSolution(const array<float,8>& Ws, const int& enLeg){
    // Probability of exiting on a particular leg
    array<float,8> prob = {0,0,0,0,0,0,0,0};

    // Compute the normalization first
    float total = 0.0;
    for (auto i=0; i!=8; i++)
        total += Ws[i];
    
    // Probability of exiting on a particular leg is proportional to its weight
    for (auto i=0; i!=8; i++)
        prob[i] = Ws[i]*Ws[enLeg]/total;
     
    return prob;
}

//*****************************************************************************
// Bounce solution based on bounce events only from the highest weight vertex
// Ws    - sorted array of weights from the heighest one to the lowest one.
//         Its indices correspond to exit legs.
// enLeg - entrance leg
//*****************************************************************************
array<float, 8> BounceSolution(const array<float,8>& Ws, const int& enLeg){
    // Probability of exiting on a particular leg
    array<float,8> prob = {0,0,0,0,0,0,0,0};

    // If the entrance leg corresponds to the highest weight vertex
    // there is a probability to switching to any other possible vertex
    if (enLeg == 0){
        float total = 0;
        for (auto i=1; i!=8; i++){
            prob[i] = Ws[i];
            total += prob[i];
        }
        prob[0] = Ws[0] - total;  // bounce probability
    }
    // Otherwise, switch to the highest weight leg with probability 1
    else prob[0] = Ws[enLeg];
     
    return prob;
}


//*****************************************************************************
// Bounce solution B based on Syljuasen's 
// "Direct Loop Updates for Quantum Lattice Models" paper
// Ws    - sorted array of weights from the heighest one to the lowest one.
//         Its indices correspond to exit legs.
// enLeg - entrance leg
//*****************************************************************************
array<float, 8> noBounceSolutionB(const array<float,8>& Ws, const int& enLeg){
    // Probability of exiting on a particular leg
    array<float,8> prob = {0,0,0,0,0,0,0,0};

    // If the entrance leg corresponds to the highest weight vertex
    if (enLeg == 0){
            prob[1] = (Ws[0] + Ws[1] - Ws[2] - Ws[3])/2.0;
            prob[2] = (Ws[0] - Ws[1] + Ws[2] - Ws[3])/2.0;
            prob[3] = Ws[3] - Ws[4]/2.0;
            for (auto leg = 4; leg!=7; leg++)
                prob[leg] = (Ws[leg] - Ws[leg+1])/2.0;
            prob[7] = Ws[7]/2.0;
    }
    // If it is the second highest 
    else if (enLeg == 1){
             prob[0] = ( Ws[0] + Ws[1] - Ws[2] - Ws[3])/2.0;
             prob[2] = (-Ws[0] + Ws[1] + Ws[2] + Ws[3])/2.0;
    }
    // If it is the third highest 
    else if (enLeg == 2){
             prob[0] = ( Ws[0] - Ws[1] + Ws[2] - Ws[3])/2.0;
             prob[1] = (-Ws[0] + Ws[1] + Ws[2] + Ws[3])/2.0;
    }
    else if (enLeg == 3){
        prob[0] = Ws[enLeg] - Ws[enLeg+1]/2.0;
        prob[enLeg+1] = Ws[enLeg+1]/2.0;
    }
    // If it is neither of the special previous cases
    else{
             // Probability to jump to the leg with the highest weight
             if (enLeg == 7) prob[0] = Ws[enLeg]/2.0;
             else            prob[0] = (Ws[enLeg] - Ws[enLeg+1])/2.0;
             
             // Probability to jump to the leg with the next highest weight
             prob[enLeg-1] = Ws[enLeg]/2.0;

             // Probability to jump to the leg with the next smallest weight
             if (enLeg != 7) prob[enLeg+1] = Ws[enLeg+1]/2.0;
    }

    return prob;
}

//*****************************************************************************
// Compute cumulative probability of exiting on one of 8 possible legs 
// on leg-switch move during the off-diagonal update 
//*****************************************************************************
int getSwitchLegP(const int& enleg, int vtype, array<bool,4>& clamped, array<float, 16>& VWeights, array<float,8>& prob){
   
    bool sdebug = false;
    prob = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    //cout << endl << "     getSwitchLegP: leg-"<< enleg << " vertex-"<< vtype << endl;
    if (VWeights[vtype] == 0.0)
        return -1;
    

    // If all legs are clamped, the vertex cannot be changed 
    bool allClamped = true;
    for (int i=0; i!=4; i++)
        if (clamped[i]==false) allClamped = false;
    if (allClamped) {
        prob = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
        return -2;
    }

    // If the entrance leg is clamped, the leg cannot be flipped
    if ( (enleg<4) and clamped[enleg])
        return -3;

    float totalW = 0;  // sum of all vertex weigths 

    // Construct all possible resulting vertices 
    array<float, 8> Ws = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    

    // Keep the id of the original vertex 
    int ovtype = vtype; 

    // Flip or not the entrance leg 
    if (enleg<4)   
        vtype = flipLeg(vtype, enleg);  

    // First add vertices with no change of the exit leg
    if (sdebug) cout << "      Legs' vertices: ";
    for (auto ileg = 0; ileg != 4; ileg++){
        Ws[4+ileg] = VWeights[vtype];
        totalW    += VWeights[vtype];
        if (sdebug) cout << vtype << " ";
    }

    // Continue with vertices having the exit leg flipped
    int ntype;
    for (auto ileg = 0; ileg != 4; ileg++){
        if (not clamped[ileg]){
            ntype = flipLeg(vtype, ileg); 
            Ws[ileg] = VWeights[ntype];
            totalW  += VWeights[ntype];
            if (sdebug) cout << ntype << " ";
        }
        else{
            Ws[ileg] = 0.0;
            if (sdebug) cout << " clamp ";
            }
    }
    
    if (sdebug){ 
        cout << endl << "      Their weights:  ";
        for (int i=0; i!= 8; i++){
           cout << Ws[i] << " ";
        }
        cout << endl; 
    }
    // Convert weights and legs to the format required by the directed loop solutions
    array<int,   8> reIndex = sort_indexes(Ws); // reindex legs such based on their weights
    array<float, 8> SortedWs; // sorted weights
    int enLegSort = -1;       // leg index in the sorted array that corresponds the entrance leg
    
    // Preserved state of the entrance leg requires a special treatement

    // Sort weights and determine the entrance leg in the new weights arrangement
    for (auto leg = 0; leg != 8; leg++){
        if (reIndex[leg] == enleg) enLegSort = leg;
        SortedWs[leg] = Ws[reIndex[leg]];
    }
    
    if (sdebug){ 
        cout << "      Sorted weights: ";
        for (int i=0; i!= 8; i++){
           cout << SortedWs[i] << " ";
        }
        cout << "      Entrance leg: " << enLegSort << endl; 
    } 
   
    // If possible, adopt the bounce-free solution B1 
    array<float, 8> tprob;
    int solType = -1;
    //if (!(SortedWs[1]+SortedWs[2]+SortedWs[3]< SortedWs[0])) {
    //    if (sdebug) cout << "     no bounce solution B1: " << SortedWs[1]+SortedWs[2]+SortedWs[3] << ">" << SortedWs[0] << endl; 
    //    tprob = noBounceSolutionB(SortedWs, enLegSort);
    //    solType = 0;
    //}
    //// Otherwise, adopt one of the two possible bounce solutions
    //else{
    //    if (2*SortedWs[0] > totalW){
    //        if (sdebug) cout << "     bounce solution: " << SortedWs[0] << ">" << totalW - SortedWs[0] << endl; 
    //        tprob = BounceSolution(SortedWs, enLegSort);
    //        solType = 1;
    //    }
    //    else{
    //        if (sdebug) cout << "     heat-bath solution" << endl; 
    //        tprob = HeatBathSolution(SortedWs, enLegSort);
    //        solType = 2;
    //    } 
    //} 
    tprob = HeatBathSolution(SortedWs, enLegSort);
    solType = 2;
 
    // Convert weights to probabilities
    for (int i=0; i!=8; i++)
        tprob[i] = tprob[i]/VWeights[ovtype];
    
    if (sdebug){ 
        cout << "      Probs: ";
        for (int i=0; i!=8; i++){
           cout << tprob[i] << " ";
        }
        cout << endl; 
    }

    // Restore the original legs labelling 
    for (auto leg=0; leg!=8; leg++)
        prob[reIndex[leg]] = tprob[leg];
    
    // Convert probability density to a cumulative distribution 
    for (auto i=1; i!=8; i++)
        prob[i] = prob[i] + prob[i-1];

    // Due to the algorithmic implementation particularity (lower_bound function),
    // all first zero-valued entries must be set to a negative number. This is done
    // in order to avoid the possibility of switching to the corresponding leg 
    // upon generation of a random number which exactly equals to 0. 
    int i=0;
    while (prob[i]==0.0){
        prob[i] = -1.0;
        i += 1;
    }

    if (sdebug){ 
        cout << "      Probs: ";
        for (int i=0; i!= 8; i++){
           cout << prob[i] << " ";
        }
        cout << endl << endl; 
    }
    return solType;
}





