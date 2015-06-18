#include <string.h>
#include <math.h>
#include "tfim.h"

TFIM::TFIM(Lattice* const _lattice, long _seed, float _beta, float _h, long _binSize):
    RandomBase(_seed),  // initialize ancestor
    communicator(_lattice, _beta, _h, _seed), 
    lattice(*_lattice) // storing a reference, no need to copy 
    //ap(*_lattice)       // slicing _lattice into its ancestor ap, a copy is created
{
    // Set up physical variables
    beta   = _beta;
    h      = _h;
    J      = 1.0;
    Nsites = lattice.getSize();
    Nbonds = lattice.getBondsN();

    // Initialize algorithmic variables
    n = 0;
    M = round(((float) Nbonds)*beta*1.5);
    lOper.resize(M);  
    lFirst.resize(Nsites,-1);
    lLast.resize(Nsites,-1);
    binSize = _binSize;
    debug = false;
    ID = communicator.getId();
    
    // Set up measurement variables
    bMperSite = false;
    bM        = true;
    aMperSite.resize(Nsites,0);
    tMperSite.resize(Nsites,0);
    an    = 0;
    aTM   = 0;
    nMeas = 0;


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
        for (int j=0; j!=Nsites; j++)
            eHeader += boost::str(boost::format("%12s%4i") %"M"%j);

    *communicator.stream("estimator")<< eHeader << endl;    

    cout << "---Initialization---" << endl;
    cout << "   J =  " << J << " h = " << h << endl;
    cout << "   Bonds #: " << Nbonds << endl;
    cout << "   Sites #: " << Nsites << endl;

}

/******************************************************************************
* Accumulate estimator measurements 
******************************************************************************/
void TFIM::Measure()
{
    nMeas += 1;
    Accumulate(n, an);

    long temp = 0;
    for (int i=0; i!=Nsites; i++){
        temp += lattice.getSpin(i);
    }
    aTM += pow((float)(temp)/(float)(Nsites),2);
    // One sufficient number of measurements are taken, record their average
    if (nMeas == binSize){

        // Record operators list length and average energy per site
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(an/(1.0*binSize)));

        float EperSite = -(float)(an)/(float)( binSize*Nsites)/beta + J*(float)(Nbonds)/(float)(Nsites) + h;
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %EperSite);
       //cout << an/(1.0*binSize) << " " << EperSite << " "; 

        // Record the total magnetization
        if (bM)
            *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(aTM/(1.0*binSize)));

        //cout << aTM/(1.0*binSize*Nsites) << " ";
        // Record magnetization per site
        if (bMperSite)
            for (int j=0; j!=Nsites; j++){
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(aMperSite[j]/(1.0*binSize)));
                //cout << aMperSite[j]/(1.0*binSize) << " ";
            }
        //cout << endl;
        *communicator.stream("estimator") << endl;    
        
        // Reset measurement variables to 0
        resetMeas();
        cout << ID << ": Measurement taken" << endl;
    }
}


/**************************************************************
* Reset measurement variables  
**************************************************************/
void TFIM::resetMeas()
{
   nMeas = 0;
   an = 0;
   if (bMperSite) aMperSite.assign(Nsites,0);
   if (bM)        aTM =0;

}


/**************************************************************
* MC diagonal update 
**************************************************************/
int TFIM::DiagonalMove()
{
    long  ibond = -1;     // bond index
    float AccP  =  0;     // acceptance probability
    int   index;          // site or bond index
    Spins ap = lattice;   // propagated in imagenary time spins state
    int s0;               // two spin states  
    int s1;
    pair<int,int> cSites;
    long nlast;
    bool ldebug = debug; 
    // if measuring magnetization per site, reset assisting variables
    if (bMperSite) {
        tMperSite.assign(Nsites,0); 
        nlast = 0;
    }

    if (ldebug){
        cout << "---Diagonal Update" << endl;
        printOperators();
        ap.printSpins();
        lattice.printSpins();
        lattice.printBonds();
    }
    //printIntVector(&lLast,  "Last");
    //printIntVector(&lFirst, "First");
    //printIntVector(&lVtx,   "Vertices");
    // Parse through each operator in the list
    for (auto oper=lOper.begin(); oper!=lOper.end(); oper++){
        nlast++;

        // If it is a null operator, try insert a diagonal one
        if (oper->type == -2){
            // Compute the probability of inserting a diagonal operator
            AccP = beta*(h*(float)(Nsites) + (2.0*J)*(float)(Nbonds))/(float) (M-n);
            if (uRand()<AccP){
                
                // If accepted, we can insert either one or two sites operator
                AccP = (h*(float)(Nsites))/(h*(float)(Nsites)+2.0*J*(float)(Nbonds));
                
                // If it is a one site operator, insert it on a random site
                if (uRand()<AccP){
                    index = uRandInt()%Nsites;    
                    if (ldebug) cout << "   (" << oper->type << "," << oper->index << ") -> (0," << index << ")" << endl;  
                    oper->set(0,index);
                    n++;
                }
                // Otherwise, insert a two sites operator on a random bond
                else{
                    index = uRandInt()%Nbonds; 
                     
                     //check whether the local spin configuration is favourable
                     cSites = lattice.getSites(index);
                     s0 = ap.getSpin(cSites.first);
                     s1 = ap.getSpin(cSites.second);
                     if (s0 + s1 != 0){
                        if (ldebug) cout << "   (" << oper->type << "," << oper->index << ") -> (1," << index << ")" << endl;  
                        oper->set(1,index);
                        n++;
                     } 
                }
                 // If there are too many non null operators in the list, 
                 // return an error code
                 if ((float)(M)/(float)(n)<1.25){     
                   cout << "M/n = " << (float) M/n << " < " << 1.25 << endl;
                   return 1;   
               }
            }
        }
        // If it is a diagonal operator, try to remove it
        else if ((oper->type == 0) or (oper->type == 1)){
            // Compute the probability of the removal process
            AccP = (float)(M-n+1)/(beta*(h*(float)(Nsites)+(2.0*J)*(float)(Nbonds)));
            if (ldebug)
                cout << "   (" << oper->type << "," << oper->index << ") -> (-2,-1)" << endl;  
            if (uRand()<AccP){
                oper->set(-2,-1);
                n--;
            }
        }
        // Otherwise it must be an off-diagonal operator. It is left as it is.
        // However the propagated spins state need to be modified.
        else{
            ap.flipSiteSpin(oper->index);
            if (ldebug){
                cout << "   offd: " << oper->index << endl;
                ap.printSpins();
                lattice.printSpins();
            }
           
            // Accumulate magnetization per site measurement 
            if (bMperSite){
                for (int i=0; i!=Nsites; i++)
                    tMperSite[i] += ap.getSpin(i)*nlast;
                tMperSite[oper->index] += -2*ap.getSpin(oper->index)*(nlast-1);
                nlast = 0;
            }
        }
    }
    // Accumulate magnetization measurements
    if (bMperSite){
        float tTM = 0;
        for (int i=0; i!=Nsites; i++){
            aMperSite[i] += (float)(tMperSite[i] + ap.getSpin(i)*nlast)/(float)(M);
            tTM += aMperSite[i];
        }
        //aTM += (1.0*tTM*tTM)/(1.0*Nsites*Nsites);
    }
    //cout << "M = " << M << " n = " << n << endl;
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
* Determine the vertex type for a particular operator acting 
* on a given bond based on its type and index 
**************************************************************/
int TFIM::VertexType(int otype, int index, Spins& ap)
{
    pair<int,int> cSites;
    int s0 = -3;
    int s1 = -3;
    bool ldebug = debug;

    if (debug){
        cout << endl << "---Vertex Type" << endl;
        //printOperators();
        ap.printSpins();
        lattice.printSpins();
        lattice.printBonds();
        cout << "   type: " << otype << "   index: " << index << endl << endl;
    }
    // If it is a two sites operator
    if (otype == 1){
        
        // Determine the spins state
        cSites = lattice.getSites(index);
        s0 = ap.getSpin(cSites.first);
        s1 = ap.getSpin(cSites.second);
        
        // Determine vertex type
        if  ((s0 ==-1) and (s1 ==-1))
            return 1;
        if  ((s0 == 1) and (s1 == 1))
            return 2;
        if  ((s0 == 1) and (s1 ==-1)){
            cout << " Vertex type: impossible configuration " << endl;
            exit(0);
            return 3;
        }
        if  ((s0 == -1) and (s1 == 1)){
            cout << " Vertex type: impossible configuration " << endl;
            exit(0);
            return 4;
        }
    }
    //if it is a one site operator
    else{
        // Determine the spin state
        s0 = ap.getSpin(index);

        // Deterimine vertex type
        if ((otype == 0) and (s0 == -1))
            return 5;
        if ((otype == 0) and (s0 == 1))
            return 6;
        if ((otype == -1) and (s0 == -1))
            return 7;
        if ((otype == -1) and (s0 == 1))
            return 8;
    }
      
    cout << " Vertex type: impossible configuration " << endl;
    exit(0);
    
    return -1;
} 

/******************************************************************************
* Vertex type to operator type conversion 
******************************************************************************/
int TFIM::VtxToOperator(int vtxType)
{
    // This convention corresponds to the one defined in VertexType method
    if ((vtxType>0) and (vtxType<5)) return  1;  
    if ((vtxType>4) and (vtxType<7)) return  0;   
    if ((vtxType>6) and (vtxType<9)) return -1;   
}


/******************************************************************************
* Determine new vertex type when flipped 
******************************************************************************/
int TFIM::VtxFlip(int vtxType, int leg)
{
    int fvtxType = -1;

    switch (vtxType){
        case 1: fvtxType = 2; break;
        case 2: fvtxType = 1; break;
        case 3: fvtxType = 4; break;
        case 4: fvtxType = 3; break;
        case 5:
            if (leg==0) fvtxType = 8;
            else        fvtxType = 7;
            break;
        case 6:
            if (leg==0) fvtxType = 7;
            else        fvtxType = 8;
            break;
        case 7:
            if (leg==0) fvtxType = 6;
            else        fvtxType = 5;
            break;
        case 8:
            if (leg==0) fvtxType = 5;
            else        fvtxType = 6;
            break;
    }
    return fvtxType;    
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
    //Reinitiate necessary data structures    
    long p=-1;            // index of current non null operator            
    int   index;          // site or bond index
    Spins ap = lattice;   // propagated in imagenary time spins state
    pair<int,int> cSites; // coordinates of 2 spins
    int cSite;            // coordinates of 1 spin    
    int nbSites;          // number of sites in a bond (2 or 1)
    bool ldebug = debug;

    fill(lFirst.begin(),lFirst.end(),-1);    //spin's 1st  leg
    fill(lLast.begin(),lLast.end(),-1);      //spin's last leg
    lLinks.assign(4*n,-1);                   //linked vertices
    lVtx.assign(n,-1);                       //vertex type

    if (ldebug){
        cout << "---Construct links " << endl;
        ap.printSpins();
        lattice.printSpins();
        printOperators();
    }
    // Start linked list construction
    for (auto oper=lOper.begin(); oper!=lOper.end(); oper++) {
        // Ignore null operators
        if (ldebug) cout << "(" << oper->type << "," << oper->index << ")" ;
        if (oper->type != -2){     
            p++;    
            index = oper->index;    
            
            // Check wether it is a one or two sites bond
            if (oper->type == 1){
                cSites = lattice.getSites(index);
                nbSites = 2;
                //cout << "Site A: " << cSites.first << " Site B: " << cSites.second << endl;
            }
            else{
                cSites = pair<int,int> (index, -1);
                nbSites = 1;
            }
            
            //Go through those sites 
            for (int a=0; a!=nbSites; a++){           
                if (a==0) cSite = get<0>(cSites);
                else      cSite = get<1>(cSites);
                
                //Has this site been already acted upon?
                if (lLast[cSite] != -1){   
                    lLinks[4*p+a]    = lLast[cSite];    //Construct a link between 2 legs
                    lLinks[lLast[cSite]] = 4*p+a;
                }
                
                //Record this first occurence 
                else{
                    lFirst[cSite] = 4*p+a; 
                }
                
                //Update the last acted leg on the site
                lLast[cSite]  = 4*p+3-a ;       
            }
           
            //Record the vertex type
            lVtx[p] = VertexType(oper->type, oper->index, ap);            
            
            // Propagate the spins state if it is an off-diagonal operator
            if  (oper->type == -1){
                ap.flipSiteSpin(index);                        
            }
            
        }

    }
    //Construct links across the periodic boundary in the p-expansion
    for (int lspin=0; lspin!=lattice.getSize(); lspin++){
        if (lLast[lspin] != -1){
            lLinks[lLast[lspin]]  = lFirst[lspin];
            lLinks[lFirst[lspin]] = lLast[lspin];
        }
    }
    //printIntVector(&lLinks, "Links");
    //printIntVector(&lLast,  "Last");
    //printIntVector(&lFirst, "First");
    //printIntVector(&lVtx,   "Vertices");

    if (ldebug){
        cout << endl << "   Construct links. " << endl;
        ap.printSpins();
        lattice.printSpins();
    }
    //cout << "    Construct links. " << endl;
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
        if (oper->type != -2){
           p++;                                         
           oper->type = VtxToOperator(lVtx[p]);
        } 
    }
    //cout << "    Map spins state back " << endl;
    // Map back the spins state     
    for (long i=0; i!=lattice.getSize(); i++){
        // Try to flip the spin if it is not acted by an operator 
        if  (lFirst[i]==-1){
            if (uRand()<0.5)
                lattice.flipSiteSpin(i);
        }                                    
        else{        
            // Get the coordinates of the first leg over this site            
            p   = (long) lFirst[i]/4; 
            leg = lFirst[i]%4;
            //cout << "Spin " << i << " Bond: " << p << " leg: " << leg << " vtx: " << lVtx[p] << " leg spin: " << LegSpin[lVtx[p]-1][leg] << endl;
            lattice.setSpin(i, LegSpin[lVtx[p]-1][leg]);
        }
    }
    //cout << "    Map spins state back. " << endl;
}


/******************************************************************************
* Map the updated state of vertices back to spins-operators state 
******************************************************************************/
int TFIM::OffDiagonalMove()
{
    stack<int> cluster;     // encountered legs during a cluster construction
    lMarked.assign(4*n,-1); // already marked legs 
    int Nclusters = 0;      // number of constructed clusters
    long  leg;               // leg full index
    long p;                 // bond index
    int  l;                 // leg index on a vertex
    bool flip;
    bool ldebug = debug;

    if (ldebug) {
        cout << "---Off-diagonal update" << endl;
        cout << "   n = " << n << endl;
        printIntVector(&lVtx, "   vert: ");
        printIntVector(&lLinks, "   link: ");
    }
    // Go through each leg
    for (int i=0; i!=4*n; i++){

        // If it belongs to a valid link 
        if (lLinks[i]!=-1){

            // And it hasn't been marked yet
            if (lMarked[i]==-1){

                // Initiate a new cluster from this leg
                Nclusters++;
                cluster.push(i);
                // Should we flip it?
                if (uRand()<0.5) flip = true;
                else             flip = false;
                
                if (ldebug){
                    cout << "   Cluster# " << Nclusters << " flipping: " << (bool) flip << endl;//<< " Legs: [" << endl;
                }
                //cout << "-- New cluster " << Nclusters << endl << " flip: " << flip << endl; 

                // Repeat while the stack of all encountered lack is not exausted
                while(!cluster.empty()){
                    leg = cluster.top();
                    cluster.pop();
                    if (ldebug) cout << "   cleg: " << leg << endl;
                    // If a previously unmarked vertex is encountered
                    if (lMarked[leg]==-1){
                        lMarked[leg] = Nclusters;
                        //if (ldebug){
                        //    cout << "   mark: "; 
                        //    for (auto mark=lMarked.begin(); mark!=lMarked.end(); mark++){
                        //        cout << *mark << " ";
                        //    }
                        //    cout << endl;
                        //}
                        // Flip it
                        p = (long) leg/4;
                        l = leg%4;
                        if (ldebug) cout << "    flipping vtx " << p << " (" <<lVtx[p] << "->" << VtxFlip(lVtx[p],l)<< ")" << endl;
                        if (flip) lVtx[p] = VtxFlip(lVtx[p],l);
                        // Mark all vertex legs if it is a two-sites diagonal operator
                        if (VtxToOperator(lVtx[p])==1)
                            for (int j=0; j!=4; j++){
                                if (j!=l){
                                   lMarked[4*p+j] = Nclusters;
                                   cluster.push(4*p+j); 
                                    //if (ldebug){
                                    //    cout << "   mark: "; 
                                    //    for (auto mark=lMarked.begin(); mark!=lMarked.end(); mark++){
                                    //        cout << *mark << " ";
                                    //    }
                                    //cout << endl;
                                    //}
                                //cout << "    adding leg " << 4*p+j << endl;
                                }
                            }
                    }
                    // Move to the next leg connected with a link
                    if (lMarked[lLinks[leg]]==-1)
                        cluster.push(lLinks[leg]);
                }
                //cout << "]" <<endl;
            }
        }
    }
    if (ldebug) printIntVector(&lVtx, "   vert: ");
    
    return Nclusters;
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
