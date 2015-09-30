#include "tfim.h"

TFIM::TFIM(Spins* const _spins, Bonds* const _bonds, 
           vector<float>* _xfield, vector<float>* _zfield,
           long _seed, float _beta, long _binSize):
    RandomBase(_seed),  // initialize ancestor
    communicator(_bonds, _beta,  _seed), 
    spins(*_spins), bonds(*_bonds),
    estimator(_spins, _bonds)
{
    // Set up physical variables
    beta   = _beta;
    Nspins = spins.getSize();
    Nbonds = bonds.getBondsN();
    xfield = _xfield;

    float  J12;
    float  h1; float h2;
    float  d1; float d2; 
    int    nA; int nB;
    int siteA; int siteB;
    tDiagOffset = 0;  
    float diagOffset; 
    debug = false;

    cout << endl <<  "=== Computing SSE vertices ===" << endl;
    vector<int>  solTypes;
    array<float, 16> VWeights;         // vertex weigths for each bond
    
    int l0; int l1; int l2; int l3;
    array<bool, 4> clampedLegs;
    // For each bond
    for (auto ibond=0; ibond!=Nbonds; ibond++){
        cout << "---Bond " << ibond << ":"<< endl; 
        // Increament main data structures
        array<array<array<float, 8>, 16*8>, 16> initP;
        VProb.push_back(initP);
        array<float, 4> initD;
        dWeights.push_back(initD);

        // Unpack bond interactions
        J12   = _bonds->getBond(ibond)->getStrength();
        
        siteA = _bonds->getBond(ibond)->getSiteA();
        siteB = _bonds->getBond(ibond)->getSiteB();
        
        nA    = _bonds->countNeighbors(siteA);
        nB    = _bonds->countNeighbors(siteB);

        h1 = 0; h2 = 0;
        if (not _zfield->empty()){
            h1  = _zfield->at(siteA)/(1.0*nA);
            h2  = _zfield->at(siteB)/(1.0*nB);
        }

        d1 = 0; d2 = 0;
        if (not _xfield->empty()){
            d1  = _xfield->at(siteA)/(1.0*nA);
            d2  = _xfield->at(siteB)/(1.0*nB);
        }

        cout << "   A: " << siteA << " B: " << siteB << " J: " << J12 << "; h: (" << h1 << ", " << h2 << "); d: (" << d1 << ", " << d2 << "); Z: (" << nA << ", " << nB <<")" << endl;

        // Redefine diagonal operators by flipping their sign
        J12 = -J12;
        h1  = -h1;
        h2  = -h2;
        
        // Get the list of possible vertices on the current bond  
        getVertices(J12, h1, h2, d1, d2, diagOffset, VWeights, dWeights[ibond]);
        
        cout << "   all  Ws: ";
        for (int j=0; j!=16; j++)
            printf("%4.2f ", VWeights[j]) ;
        cout << endl;

        cout << "   diag Ws: ";
        for (int j=0; j!=4; j++)
            printf("%4.2f ", dWeights[ibond][j]) ;
        cout << endl;

        cout << "   C: " << diagOffset << endl;
        
        // Accumulate the total diagonal offset
        tDiagOffset += diagOffset;      

        // Get the leg-switch probabilities look-up table for each vertex and leg
        if (debug) cout << "=== Computing off-diagonal moves probabilities ===" << endl;
        
        // For each combination of clamped legs
        for (auto iclamp=0; iclamp!=16; iclamp++){
            getAllLegs(iclamp, l0, l1, l2, l3);
            clampedLegs = { l0 == 1, l1 == 1, l2 == 1, l3 == 1};
            for (auto vtype=0; vtype!=16; vtype++)
                for (auto enleg=0; enleg!=8; enleg++)
                    solTypes.push_back(getSwitchLegP(enleg, vtype, clampedLegs, VWeights, VProb[ibond][iclamp][vtype*8+enleg]));
        }
    }

    cout << endl <<  "=== Solutions to the directed loop equations ===" << endl;
    float stopP;
    float startP;
    int iclamp = 15;
    for (auto ibond=0; ibond!=1; ibond++){
        cout << endl << "===Bond " << ibond << "=== " << endl;
        for (auto vtype=0; vtype!=16; vtype++)
            if (VWeights[vtype]!=0.0){
                cout << "---Vertex " << vtype << ": " << endl;
                for (auto enleg=0; enleg!=8; enleg++){
                    cout << "   Leg " << enleg;
                        if (solTypes[vtype*8+enleg]==0) cout << " (N-B): ";  // no-bounce solution
                        if (solTypes[vtype*8+enleg]==1) cout << " (L-B): ";  // bounce from the largest weight solution
                        if (solTypes[vtype*8+enleg]==2) cout << " (H-B): ";  // heat-bath bounce solution
                    
                    for (auto exleg=0; exleg!=8; exleg++)
                        printf("%3.2f; ", VProb[ibond][iclamp][vtype*8+enleg][exleg]);
                    
                    if (enleg<4){
                        if (VProb[ibond][iclamp][vtype*8+enleg][3] != -1.0) stopP = 1.0 - VProb[ibond][iclamp][vtype*8+enleg][3];
                        else                                       stopP = 1.0;
                        printf(" Stop  P: %3.2f \n", stopP);
                    }
                    else{
                        if (VProb[ibond][iclamp][vtype*8+enleg][3] != -1.0) startP = VProb[ibond][iclamp][vtype*8+enleg][3];
                        else                                       startP = 0.0;
                        printf(" Start P: %3.2f \n", startP);
                    }
                }
            }
    } 
    // Initialize algorithmic variables
    n = 0;
    M = round(((float) Nbonds)*beta*1.5);
    if (M<1) M = 1;
    lOper.resize(M);  
    lFirst.resize(Nspins,-1);
    lLast.resize(Nspins,-1);
    binSize = _binSize;
    ID = communicator.getId();
    nLoops = max(1, (int) M/10);

    // Set up measurement variables
    bMperSite = false;
    bM        = true;
    aMperSite.resize(Nspins,0);
    tMperSite.resize(Nspins,0);
    an    = 0;
    aTM   = 0;
    nMeas = 0;


    // Compute cumulative probabilities of diagonal operators insertation
    computeDiagProb();

    // Define the following map:
    // LegSpin[vertex_type][leg_index] -> spin state
    int tmp[8][4]{  {-1,-1,-1,-1},
                    { 1, 1, 1, 1},
                    { 1,-1,-1, 1},
                    {-1, 1, 1,-1},

                    {-1, 0, 0,-1},
                    { 1, 0, 0, 1},
                    {-1, 0, 0, 1},
                    { 1, 0, 0,-1}
                };
    memcpy(LegSpin,tmp,sizeof tmp);

    //Write headers for a new estimator file
    string eHeader = "";
    eHeader += boost::str(boost::format("#%15s%16s")%"nT"%"ET"); 
    if (bM)        eHeader += boost::str(boost::format("%16s")%"MT^2"); 
    if (bMperSite)
        for (int j=0; j!=Nspins; j++)
            eHeader += boost::str(boost::format("%12s%4i") %"M"%j);

    *communicator.stream("estimator")<< eHeader; 

    // complete it with the headers for local averages and correlators
    *communicator.stream("estimator")<< estimator.Xheader 
                                     << estimator.Zheader
                                     << estimator.Iheader
                                     << estimator.Xheader0s 
                                     << estimator.Zheader0s
                                     << estimator.Iheader0s
                                     << endl;    
    cout << estimator.Xheader 
         << estimator.Zheader
         << estimator.Iheader
         << endl;    

}

/*****************************************************************************
 * Compute the look-up table for diagonal vertices insertion probability
 *****************************************************************************/
void TFIM::computeDiagProb()
{
    // Commulative insertion probability 
    diagProb.clear();

    // Compute the normalization
    totalWeight = 0.0;
    for (auto b=0; b!=Nbonds; b++)
        for (auto i=0; i!=4; i++)
            totalWeight += dWeights[b][i];
        
    // Compute the cummulative distribution
    float runTotal = 0;
    for (auto b=0; b!=Nbonds; b++)
        for (auto i=0; i!=4; i++){
            runTotal += dWeights[b][i];
            diagProb.push_back(runTotal/totalWeight);
        }
   
    
    cout << endl << "=== Diagonal insertion probabilities ===" << endl << "   ";
    for (auto dp=diagProb.begin(); dp!=diagProb.end(); dp++)
        printf("%4.3f; ", *dp);
    cout << endl << endl;

}
 


/**************************************************************
* MC diagonal update 
**************************************************************/
int TFIM::DiagonalMove()
{
    int   index;        // array index
    int   ibond;        // bond index
    int   vtype;        // vertex type
    Spins ap = spins;   // propagated in imagenary time spins state
    int s0;   int s1;   // two spins state  
    int is0;  int is1;  // two spins indices 
    
    bool binsert;        // flag indicating the state of an operator insertation 
    long ioper = 0;      // counter of non-null operators in the list

    bool ldebug = false;
    if (debug)  ldebug = true;
   
    vector<vector<long>> worldLines;  // kinks in the spin worldlines caused
    worldLines.resize(Nspins);        // by off-diagonal operator sigma_x  

    //vector<int> SliceVertices; // list of compatible diagonal vertices for the zero's slice 
    //float SliceWeight = 0;     // the total weight of those vertices.
    //
    //// Compute it
    //for (auto b=0; b!=Nbonds; b++){
    //    is0 = bonds.getBond(b)->getSiteA();
    //    is1 = bonds.getBond(b)->getSiteB();
    //    s0 = ap.getSpin(is0);
    //    s1 = ap.getSpin(is1); 
    //    SliceVertices.push_back(getCompatibleDiagVertex(s0, s1));
    //    SliceWeight += dWeights[b][SliceVertices[b]];
    //}

    if (ldebug){
       cout << "---Diagonal move" << endl;
       cout << "   Spins: " << endl << "   ";
       for (int i=0; i!=Nspins; i++)
          cout << spins.getSpin(i) << " ";
       cout << endl;

       //cout << "   Allowed vertices: total weight=" << SliceWeight << endl << "   ";
       //for (int i=0; i!=Nbonds; i++)
       //    cout << SliceVertices[i] << " ";
       //cout << endl;

       printOperators();
    }

    

    // Parse through each operator in the list
    for (auto oper=lOper.begin(); oper!=lOper.end(); oper++){

        // If it is a null operator, try to insert a diagonal one
        if (oper->type == 2){
            // Compute the probability of inserting a diagonal operator
            //if (uRand() < (beta*SliceWeight/(float) (M-n)) ){
            if (uRand() < (beta*totalWeight/(float) (M-n)) ){
                // Repeat until a compatible vertex is selected 
                binsert = false;
                //do{
                    // Randomly choose bond and type of the vertex to be inserted 
                    index = lower_bound(diagProb.begin(), diagProb.end(), uRand()) - diagProb.begin();
                    ibond = (int) index / 4; 
                    vtype = index % 4;

                    // Get the spins state at the potential bond
                    s0 = ap.getSpin(bonds.getBond(ibond)->getSiteA()); 
                    s1 = ap.getSpin(bonds.getBond(ibond)->getSiteB()); 

                    // Insert new operator only if the local spin configuration is favourable
                    if (getCompatibleDiagVertex(s0,s1) == vtype){
                       oper->set(1,ibond);
                       n++;
                       binsert = true;
                    }
                //} while( !binsert );
                
                if (binsert) ioper += 1;

                // If there are too many non null operators in the list, return an error code
                if ((float)(M)/(float)(n)<1.25){     
                  cout << "M/n = " << (float) M/n << " < " << 1.25 << endl;
                  return 1;   
                }
            }
        }
        // If it is a diagonal operator, try to remove it
        else if (oper->type == 1) {
            // Compute the probability of the removal process
            //if (uRand() < ( (float)(M-n+1)/(beta*SliceWeight) )){
            if (uRand() < ( (float)(M-n+1)/(beta*totalWeight) )){
                oper->set(2);
                n--;
            }
            else{
                ioper += 1;
            }
        }
        // Otherwise it must be an off-diagonal operator. It is left as it is.
        // However the propagated spins state need to be modified.
        else{
            
            // Get coordinates of the effected spins
            is0 = bonds.getBond(oper->index)->getSiteA(); 
            is1 = bonds.getBond(oper->index)->getSiteB(); 

            ioper += 1;
            
            // Update the propagated spins state
            if (oper->type == 0 ){
                ap.flip(is0); 
                worldLines[is0].push_back(ioper);
            }
            else{
                ap.flip(is1); 
                worldLines[is1].push_back(ioper);
            } 
            // Get spins' state
            s0 = ap.getSpin(is0);
            s1 = ap.getSpin(is1);

            // Update the list of compatible vertices at the updated slice and their total weight
            //SliceWeight -= dWeights[oper->index][SliceVertices[oper->index]];
            //SliceVertices[oper->index] = getCompatibleDiagVertex(s0, s1);
            //SliceWeight += dWeights[oper->index][SliceVertices[oper->index]];
        }
    }
    estimator.Accumulate(worldLines, n);

    if (ldebug) printOperators();
    return 0;
}

/**************************************************************
* Increase length of the operator list 
**************************************************************/
void TFIM::AdjustM()
{
    M = M+7*long(sqrt(n));
    lOper.resize(M);
    resetMeas();
}

/**************************************************************
* Increase length of the operator list 
**************************************************************/
bool TFIM::AdjustLoopsN()
{
    nFlippedLegs = 0;
    nBounces     = 0;
    float visitedFrac = 0.0;

    for (int i=0; i!=binSize; i++){
        while (DiagonalMove()==1) AdjustM();
              
        ConstructLinks();         // linked list and vertex list construction 
        OffDiagonalMove();        // loop update
        MapStateBack();           // mapping back the updated state 
        visitedFrac += (1.0*nFlippedLegs)/(0.5*M);
        nFlippedLegs = 0;
    }

    if (visitedFrac/(binSize*1.0) > 1.0) return false;
    else{
        nLoops += max((int) M/20, 1);
        cout << ID << ": Flipped fraction=" << setw(1) << setprecision(2) << setfill('0') << visitedFrac/(binSize*1.0) << " 0.5xM=" << 0.5*M ;
        cout << " Increasing loops # " << nLoops << endl;
        return true;
    }
}


/******************************************************************************
* Accumulate estimator measurements 
******************************************************************************/
void TFIM::Measure()
{
    nMeas += 1;
    Accumulate(n, an);

    long temp = 0;
    for (int i=0; i!=Nspins; i++){
        temp += spins.getSpin(i);
    }
    aTM += pow((float)(temp)/(float)(Nspins),2);
    // Once sufficient number of measurements are taken, record their average
    if (nMeas == binSize){

        // Record operators list length and average energy per site
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(an/(1.0*binSize)));

        float EperSite = -(float)(an)/(float)( binSize*Nspins)/beta + tDiagOffset/(float)(Nspins);
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %EperSite);
       //cout << an/(1.0*binSize) << " " << EperSite << " "; 

        // Record the total magnetization
        if (bM)
            *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(aTM/(1.0*binSize)));

        //cout << aTM/(1.0*binSize*Nspins) << " ";
        // Record magnetization per site
        if (bMperSite)
            for (int j=0; j!=Nspins; j++){
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(aMperSite[j]/(1.0*binSize)));
                //cout << aMperSite[j]/(1.0*binSize) << " ";
            }
        
        // Record local expectation values and spin-spin correlation functions
        int i=0;
        for (auto &Sx: *(estimator.getSxs())){
             *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(-1.0*Sx/(1.0*estimator.getCount()*xfield->at(i)*beta)));
             i++;
        } 
        
        for (auto &Sz: *(estimator.getSzs())){
            *communicator.stream("estimator")  << boost::str(boost::format("%16.8E") %(Sz/(1.0*estimator.getCount())));
        } 
       
        for (auto &SzSz: *(estimator.getSzSzs())){
            *communicator.stream("estimator")  << boost::str(boost::format("%16.8E") %(SzSz/(1.0*estimator.getCount())));
        } 
       
        // Record local expectation values and spin-spin correlation functions for the first slice only
        i=0;
        for (auto &Sx: *(estimator.getSxs0s())){
             *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(-1.0*Sx/(1.0*estimator.getCount()*xfield->at(i)*beta)));
             i++;
        } 
        
        for (auto &Sz: *(estimator.getSzs0s())){
            *communicator.stream("estimator")  << boost::str(boost::format("%16.8E") %(Sz/(1.0*estimator.getCount())));
        } 
       
        for (auto &SzSz: *(estimator.getSzSzs0s())){
            *communicator.stream("estimator")  << boost::str(boost::format("%16.8E") %(SzSz/(1.0*estimator.getCount())));
        } 
 
        //cout << endl;
        *communicator.stream("estimator") << endl;    
        
        // Reset measurement variables to 0
        resetMeas();
        estimator.Reset();
        
        cout << ID << ": Measurement taken" << endl;
    }
}

/******************************************************************************
* Accumulate estimator measurements
******************************************************************************/
void TFIM::Measure(int sampleInd, double* aEnergy, double* aMagnetization)
{
    nMeas += 1;
    Accumulate(n, an);

    long temp = 0;
    for (int i = 0; i != Nspins; i++) {
        temp += spins.getSpin(i);
    }
    aTM += pow((float)(temp) / (float)(Nspins), 2);
    // One sufficient number of measurements are taken, record their average
    if (nMeas == binSize) {
        float EperSite = -(float)(an) / (float)(binSize * Nspins) / beta + tDiagOffset / (float)(Nspins);
        aEnergy[sampleInd] = EperSite;

        // Record the total magnetization
        if (bM)
            aMagnetization[sampleInd] = aTM / (1.0 * binSize);
        
        resetMeas();
    }
}

/**************************************************************
* Reset measurement variables  
**************************************************************/
void TFIM::resetMeas()
{
   nMeas = 0;
   an = 0;
   if (bMperSite) aMperSite.assign(Nspins,0);
   if (bM)        aTM =0;
    nFlippedLegs = 0;   
    nBounces     = 0;   

}


/******************************************************************************
* Construct linked list and vertices in preparation for the off-diagonal move
******************************************************************************/
void TFIM::ConstructLinks()
/* 
   Coordiantes of a leg on a operator are encoded  by one number: 
   4*operator_index+leg_index (1 of 4). Both, the index and the value,
   of an element in the linked list carry an information: the index
   encodes leg coordinates while its value encodes coordiantes of
   another leg it points to. 
*/

{
    long p=0;        // index of current non-null operator            
    Spins ap = spins; // propagated in imagenary time spins state
    array<int,2> is;  // indices of two spins

    // Reset main data structures
    fill(lFirst.begin(),lFirst.end(),-1);    // spin's 1st  leg
    fill(lLast.begin(),lLast.end(),-1);      // spin's last leg
    lLinks.assign(4*n,-1);                   // linked vertices
    lVtx.assign(n,-1);                       // vertex type
    lBond.assign(n,-1);                      // bond indes associated with a vertex
    lClamped.assign(n, {false, false, false, false, false});

    bool ldebug = false;

    if (ldebug){
        cout << "---Construct links " << endl;
        ap.print();
        printOperators();
    }
    // Start linked list construction
    for (auto oper=lOper.begin(); oper!=lOper.end(); oper++) {
        //if (ldebug) cout << "(" << oper->type << "," << oper->index << ")" ;
        
        // Ignore null operators
        if (oper->type != 2){     
            
            // Get spins' indices affected by the operator
            is[0] = bonds.getBond(oper->index)->getSiteA(); 
            is[1] = bonds.getBond(oper->index)->getSiteB(); 
            
            // Go through legs associated with those spins at the current slice
            for (int i=0; i!=2; i++){
                // If this is not a first slice leg
                if (lLast[is[i]] != -1){
                   
                    //Construct a new link  
                   lLinks[4*p + i]      = lLast[is[i]];    
                   lLinks[lLast[is[i]]] = 4*p + i;
                }
                // Otherwise, remember its coordinate 
                else{
                    lFirst[is[i]] = 4*p+i; 
                    // If clamped, mark the vertex
                    if (spins.isClamped(is[i])){
                        lClamped[p][0]   = true;    // mark vertex as clamped
                        lClamped[p][i+1] = true;    // mark leg as clamped
                    }
                }

                // Always update the last leg
                lLast[is[i]]  = 4*p+3-i ;   
            }        
            
            //Record the vertex type and associated bond index
            lVtx[p]  = getVertexID(ap.getSpin(is[0]), ap.getSpin(is[1]), oper->type);            
            lBond[p] = oper->index;

            // Propagate the spins state if it is an off-diagonal operator
            if      (oper->type ==  0) ap.flip(is[0]); 
            else if (oper->type == -1) ap.flip(is[1]); 
            
        p++;    
        }
    }
    
    //Construct links across the periodic boundary in the p-expansion
    int l;
    for (int lspin=0; lspin!=spins.getSize(); lspin++){
        if (lLast[lspin] != -1){
            // If clamped, mark the vertex
            if (spins.isClamped(lspin)){
                p = (long) lLast[lspin]/4;  // vertex coordinate
                l =        lLast[lspin]%4;  // leg coordinate
                lClamped[p][0]   = true;    // mark vertex as clamped
                lClamped[p][l+1] = true;    // mark leg as clamped
            }
            // If clamped, there is no need to construct a PBC link
            else{
                lLinks[lLast[lspin]]  = lFirst[lspin];
                lLinks[lFirst[lspin]] = lLast[lspin];
            }
        }
    }
    //printIntVector(&lLinks, "Links");
    //printIntVector(&lLast,  "Last");
    //printIntVector(&lFirst, "First");
    //printIntVector(&lVtx,   "Vertices");

    if (ldebug){
        //cout << endl << "   Construct links. " << endl;
        ap.print();
        //ham.bonds.print();
    }
    //cout << "    Construct links. " << endl;
}


/******************************************************************************
* Map the updated state of vertices back to spins-operators state 
******************************************************************************/
int TFIM::OffDiagonalMove()
{
    long enleg; // entrance leg
    long exleg; // exit leg
    int  l;     // full leg's coordinates 
    long p;     // operator index
    float accP = 0;
    bool ldebug = false;
    if (debug)  ldebug = true;
    int clampid;

    if (ldebug) {
        cout << "---Off-diagonal update" << endl;
        cout << "   n = " << n << " loops # = " << nLoops << endl;
        printIntVector(&lVtx,   "   vert: ");
        printIntVector(&lLinks, "   link: ");
        printIntVector(&lBond,  "   bond: ");
        cout << "   clamped: "; 
        for (int i=0; i!=n; i++){
            cout << "(";
            for (int j=0; j!=5; j++)
                cout << lClamped[i][j] << " ";
            cout << ") ";
        }
        cout << endl;
    }
    if (n==0) return 1;

    array<float, 8>* SwitchLegP; 
  
    // Build a pre-set number of loops  
    for (long i=0; i!=nLoops; i++){ 
        // Randomly choose the starting leg out of all legs 
        l = uRandInt()%(4*n);

        // Make a vertex move
        p   = (long) l/4;
        enleg = l%4;
        if (lClamped[p][0]){
           clampid = getVertexID(((int) lClamped[p][1])*2-1, 
                                 ((int) lClamped[p][2])*2-1, 
                                 ((int) lClamped[p][3])*2-1, 
                                 ((int) lClamped[p][4])*2-1);
            SwitchLegP = &(VProb[lBond[p]][clampid][lVtx[p]*8+enleg+4]);  // the entrance leg is not flipped (+4)
        } 
        else
            SwitchLegP = &(VProb[lBond[p]][0][lVtx[p]*8+enleg+4]);  // the entrance leg is not flipped (+4)
        
        
        accP = uRand();
        exleg = lower_bound(SwitchLegP->begin(), SwitchLegP->end(), accP) - SwitchLegP->begin();
        if (ldebug){
           cout << "   Loop " << i << endl;
           cout << "   l=" << l << " p=" << p << " bond=" << lBond[p] << " clamped=" << lClamped[p][0] << " enleg=" << enleg << " exleg=" << exleg<< " vtx: " << lVtx[p] << " -> "; 
        }
        // Modify the vertex id
        lVtx[p] = FlipVertex(lVtx[p], enleg+4, exleg);
        if ((lVtx[p] == 3) or (lVtx[p] == 5) or (lVtx[p] == 10) or (lVtx[p] == 12)){
           cout << "ERROR: non-allowed vertex is encountered" << endl;
           exit(0);
        }
        if (ldebug){
            cout << lVtx[p] << endl << "   RN=" << accP << " Switch P: ";
            for (int i=0; i!=8; i++)
                cout << SwitchLegP->at(i)<< " ";
            cout << endl;
        }
        if ((lVtx[p] == 3) or (lVtx[p] == 5) or (lVtx[p] == 10) or (lVtx[p] == 12)){
           cout << endl << "ERROR: non-allowed vertex is encountered" << endl;
           exit(0);
        }

        // Unless the new leg is unchanged, keep building the loop
        while (exleg<4){
            // Go to the link connected leg on the next vertex
            l = p*4 + exleg;
            l = lLinks[l];

            // Make a vertex move
            p   = (long) l/4;
            enleg = l%4;
            if (lClamped[p][0]){
               clampid = getVertexID(((int) lClamped[p][1])*2-1, 
                                     ((int) lClamped[p][2])*2-1, 
                                     ((int) lClamped[p][3])*2-1, 
                                     ((int) lClamped[p][4])*2-1); 
                SwitchLegP = &(VProb[lBond[p]][clampid][lVtx[p]*8+enleg]);  // the entrance leg is not flipped (+4)
            } 
            else
                SwitchLegP = &(VProb[lBond[p]][0][lVtx[p]*8+enleg]);  // the entrance leg is not flipped (+4)
            
            accP = uRand();
            exleg = lower_bound(SwitchLegP->begin(), SwitchLegP->end(), accP) - SwitchLegP->begin();
            if (ldebug)
               cout << "      l=" << l << " p=" << p << " bond=" << lBond[p] << " clamped=" << lClamped[p][0]  << " enleg=" << enleg << " exleg=" << exleg<< " vtx: " << lVtx[p] << " -> "; 

            // Modify the vertex id
            lVtx[p] = FlipVertex(lVtx[p], enleg, exleg);

            if (ldebug){
                cout << lVtx[p] << endl << "      RN=" << accP << " Switch P: ";
                for (int i=0; i!=8; i++)
                    cout << SwitchLegP->at(i)<< " ";
                cout << endl;
            }

            if ((lVtx[p] == 3) or (lVtx[p] == 5) or (lVtx[p] == 10) or (lVtx[p] == 12)){
               cout << endl << "ERROR: non-allowed vertex is encountered" << endl;
               exit(0);
            }
           // Gather loop  statistics;
            if (exleg != enleg) nFlippedLegs += 2;
            else                nBounces     += 1; 
        } 
    } 
    return 0;
}


/******************************************************************************
* Map the updated state of vertices back to spins-operators state 
******************************************************************************/
void TFIM::MapStateBack()
{
    int p   = -1;     //non null operator index
    int leg = -1;     //leg index
    // Map back the vertices state to the operators state
    for (auto oper=lOper.begin(); oper!=lOper.end(); oper++){
        if (oper->type != 2){
           p++;                                         
           oper->type = VertexToOperator(lVtx[p]);
        } 
    }
    //cout << "    Map spins state back " << endl;
    // Map back the spins state     
    for (long i=0; i!=Nspins; i++){
        // Try to flip the spin if it is not acted by an operator nor is clamped
        if  (lFirst[i]==-1){ 
             if (not spins.isClamped(i))
                if (uRand()<0.5) spins.flip(i);
        }                                    
        else{        
            // Get the coordinates of the first leg over this site            
            p   = (long) lFirst[i]/4; 
            leg = lFirst[i]%4;
            //cout << "Spin " << i << " Bond: " << p << " leg: " << leg << " vtx: " << lVtx[p] << " leg spin: " << LegSpin[lVtx[p]-1][leg] << endl;
            spins.setSpin(i, getLeg(lVtx[p], leg));
        }
    }
    //cout << "    Map spins state back. " << endl;
}



/******************************************************************************
* Print vector of integers
******************************************************************************/
void TFIM::printIntVector(vector<int>* vect, string name )
{
    cout << name;
    for (auto item=vect->begin(); item!=vect->end(); item++)
        cout << *item << " ";
    cout << endl;
}

/******************************************************************************
* Print operators list
******************************************************************************/
void TFIM::printOperators()
{
    cout << "   Oper: ";
    for (auto oper=lOper.begin(); oper!=lOper.end(); oper++)
        cout << "(" << oper->type << "," << oper->index << ") ";
    cout << endl;
}

/******************************************************************************
* Accumulate long type variable
******************************************************************************/
void TFIM::Accumulate(long from, long& to)
{
    to += from;
}

/******************************************************************************
* Accumulate float type variable
******************************************************************************/
void TFIM::Accumulate(float from, float& to)
{
    to += from;
}

/******************************************************************************
* Accumulate vector type variable
******************************************************************************/
void TFIM::Accumulate(vector<float>& from, vector<float>& to)
{
    for (int i=0; i!=to.size(); i++){
        to[i] += from[i];
    }
}
