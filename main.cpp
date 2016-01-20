#include <string.h>
#include <boost/format.hpp> 
#include <boost/program_options.hpp>
#include "tfim.h" 
#include <iostream>

using namespace std;
namespace po = boost::program_options;



string getOptionsList(const vector<string> options) {
    ostringstream optionList;
    copy(options.begin(), options.end() - 1, ostream_iterator<string>(optionList, ", "));
    optionList << options.back();

    return optionList.str();
}

void fReadBond(fstream* bfile, int& siteA, int& siteB, float& strength)
{
    string         sbuf;
    istringstream ssbuf;
    getline(*bfile, sbuf);
    ssbuf.str(sbuf);
    ssbuf >> siteA;
    ssbuf >> siteB;
    ssbuf >> strength;
}



int main(int argc, char *argv[])
{
    vector<string> lLatticeName;
        lLatticeName.push_back("chimera");
        lLatticeName.push_back("rectangle");
        lLatticeName.push_back("general");
    string latticeNames = getOptionsList(lLatticeName);


    po::options_description cmdLineOptions("Command line options"); 
    po::options_description simulationOptions("Simulation options"); 
    po::options_description physicalOptions("Physical parameters options"); 
    po::variables_map params;
    simulationOptions.add_options()
            ("help,h", "produce help message")
            ("state,s",      po::value<string>(),"path to the state file")
            ("process_id,p", po::value<long>()->default_value(0),"random generator seed")
            ("equil,e",      po::value<long>()->default_value(100),"number of equilibration sweeps to take")
            ("meas,m",       po::value<long>(),"number of measurements to take")
            ("binsize",      po::value<long>()->default_value(100),"number of MC sweeps in one measurement (default = 100)")
            ;
    physicalOptions.add_options()
            ("temp,T",       po::value<float>(),"temperature")
            ("beta,B",       po::value<float>(),"inverse temperature")
            ("lattice,L",    po::value<string>()->default_value("rectangle"), boost::str(boost::format("interaction potential type {%s}") % latticeNames).c_str())
            ("OBC",          po::value<bool>()->default_value(0), "set OBC on a rectangular lattice (default false)")
            ("clamped,C",    po::value<string>(),"path to the file with the state of clamped qubits")
            ("inter,I",      po::value<string>(),"path to the file with interaction values")
            ("width,X",      po::value<int>(),"lattice width")
            ("height,Y",     po::value<int>(),"lattice height")
            ("unity",        po::value<int>(),"lattice unit cell height")
            ("long,H",       po::value<float>(),"longitudinal field strength")
            ("trans,D",      po::value<float>(),"transverse field strength (must be positive)")
            ("ising,J",      po::value<float>(),"Ising term interaction strength")
            ;
    cmdLineOptions.add(simulationOptions).add(physicalOptions);
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), params);
    po::notify(params);

    if (params.count("help")) {
       cout << cmdLineOptions << "\n";
       return 1;
    } 
  
    float beta = -1.0;
    if  ((params.count("temp")) && (params.count("beta"))){
        cerr << "Error: simultanious definition of temperature via T and beta parameters" << endl;
        return 1; 
    }
    else if  (!(params.count("temp")) && !(params.count("beta"))){
        cerr << "Error: define simulation temperature (T,beta)" << endl;
        return 1;
    } 
    else if (params.count("temp"))  beta = 1.0/params["temp"].as<float>();
    else                            beta = params["beta"].as<float>();

    if  (!(params.count("meas"))){
        cerr << "Error: define the number of measurements to take (m)" << endl;
        return 1; 
    }


    // ------------------------------------------------------------------------ 
    // Initialize lattice geometry
    // ------------------------------------------------------------------------ 
    if (find(lLatticeName.begin(), lLatticeName.end(), params["lattice"].as<string>()) == lLatticeName.end()) {
       cerr << endl << "ERROR: invalid lattice type" << endl << endl;
       cerr << "Action: set it (L) to one of:" << "\t[" << latticeNames << "]" <<  endl;
       return 1;
    }
    
    int width  = 0;
    int height = 0;
    if  (!(params["lattice"].as<string>() == "general")){
        if  (!(params.count("width"))){
            cerr << "Error: define lattice (x)" << endl;
            return 1; 
        }
        if  (!(params.count("height"))){
            cerr << "Error: define lattice height (y)" << endl;
            return 1; 
        }
        width = params["width"].as<int>();
        height = params["height"].as<int>();
    }
    else 
        if (!(params.count("inter"))){
            cout << "Error: for a general lattice type, a file with interactions must be specified (I)" << endl;
            return 1;
        }

 


    int unitHeight = 0;
    int unitWidth  = 0;
    if (params["lattice"].as<string>() == "chimera"){
        if  (!(params.count("unity"))){
            cerr << "Error: define unit cell height (unity)" << endl;
            return 1; 
        }
        unitHeight = params["unity"].as<int>(); 
        unitWidth  = 2;
    }
    else if (params["lattice"].as<string>() == "rectangle"){
        unitHeight = 1; 
        unitWidth  = 1;
    }

    // ------------------------------------------------------------------------ 
    // Initialize interactions
    // ------------------------------------------------------------------------ 
    //if  (!(params.count("trans")) and !(params.count("inter"))){
    //    cerr << "Error: specify the transverse field (D or I)" << endl;
    //    return 1; 
    //}
    //else if (params.count("trans") and params.count("inter")){
    //    cerr << "Error: conflicting flags, remove one of them (D and I)" << endl;
    //    return 1; 
    //}

    
    if  (!(params.count("ising")) and !(params.count("inter"))){
        cerr << "Error: specify the Ising interaction strength (J or I)" << endl;
        return 1; 
    }
    else if  (params.count("ising") and params.count("inter")){
        cerr << "Error: conflicting flags, remove one of them (J or I)" << endl;
        return 1; 
    }
 

    Bonds* bonds;
    vector<float> xfield; 
    vector<float> zfield; 
    // If interactions are specified in a file, read them off
    if (params.count("inter")){
        
        bonds = new Bonds();
        cout << "Loading interaction values from: " << params["inter"].as<string>() << endl;
    
        fstream fInter (params["inter"].as<string>(), ios_base::in);
        // Buffer variables
        string         sbuf;

        // Read the first line containing no information
        getline(fInter, sbuf);

        // Read the second line with the header
        getline(fInter, sbuf);
        istringstream  ssbuf(sbuf.substr(1));

        // Record the number of various fields, bonds
        int nSz = 0;       int nSx=0;    int nSzSz=0;
        ssbuf >> nSz; ssbuf >> nSx; ssbuf >> nSzSz;
        
        // Read off interaction values line by line
        int   siteA;
        int   siteB;
        float strength;
        cout << "---Sz: " << endl;
        for (int i=0; i!=nSz; i++){
            fReadBond(&fInter, siteA, siteB, strength);
            zfield.push_back(strength);
            cout << "   Site " << siteA << " = " << strength << endl;
        }
        
        cout << "---Sx: " << endl;
        for (int i=0; i!=nSx; i++){
            fReadBond(&fInter, siteA, siteB, strength);
            if (strength < 0){
                cout << "Error: no negative values for the transverse field are allowed" << endl;
                return 1;
            }
            xfield.push_back(strength);
            cout << "   Site " << siteA << " = " << strength << endl;
        }

        cout << "---SzSz: " << endl;
        for (int i=0; i!=nSzSz; i++){
            fReadBond(&fInter, siteA, siteB, strength);
            bonds->setBond(siteA, siteB, strength);
            cout << "   Bond (" << siteA << "," << siteB << ") = " << strength << endl;
        }
    }
    // Otherwise, initialize them from pre-set nearest neighbours interactions
    else{
        int Nspins = 0;
        if (params["lattice"].as<string>() == "rectangle"){
            bonds = new Rectangle(width, height, params["OBC"].as<bool>(), params["ising"].as<float>());
            Nspins = width*height;
        }
        else{                    
            bonds = new Chimera(width, height, unitHeight, params["ising"].as<float>());
            Nspins = width*height*unitWidth*unitHeight;
        }
        if (params["trans"].as<float>()<0){
           cout << "Error: no negative values for the transverse field are allowed" << endl;
           return 1;
        } 
        if (params.count("trans")) xfield.resize(Nspins, params["trans"].as<float>());
        if (params.count("long"))  zfield.resize(Nspins, params["long"].as<float>());
    } 

    // ------------------------------------------------------------------------ 
    // Initialize spins 
    // ------------------------------------------------------------------------ 
    Spins * spins;
    
    // Determine the number of spins
    int Nspins = -1;
    if (params["lattice"].as<string>() == "rectangle") Nspins = width*height;
    if (params["lattice"].as<string>() == "chimera")   Nspins = width*height*unitWidth*unitHeight;
    if (params["lattice"].as<string>() == "general"){
        Nspins = xfield.size();
        if (Nspins==0) Nspins = zfield.size();
        if (Nspins==0){
           cout << "Error: cannot extract the size of the general lattice from " << params["inter"].as<string>() << endl;
           return 1;
        }
    }
    
    // Initialize them from a pre-set configuration if needed
    if (not params.count("clamped")){
        spins = new Spins(Nspins, params["process_id"].as<long>());
    }
    else{     
        cout << "Loading the state of clamped qubits from: " << params["clamped"].as<string>() << endl;
    
        fstream fClamped(params["clamped"].as<string>(), ios_base::in);
        string         sbuf;
        istringstream  ssbuf;
    
        // Read off the first line only 
        getline(fClamped, sbuf);
        ssbuf.str(sbuf);
    
        vector<int> clamped;
        int spin;
        while (ssbuf >> spin){
              if (spin!=1 and spin!=-1 and spin!=0){ 
                 cout << "Error: impossible qubit state " << spin << " encountered" << endl;
                 return 1;
              }
              clamped.push_back(spin);
        }
        if (clamped.size() != Nspins){
           cout << "Error: number of qubits "<< clamped.size() << " is inconsistent with the hamiltonian" << Nspins << endl;
           return 1;
        }
        spins = new Spins(Nspins, params["process_id"].as<long>(), clamped);
    }
        


    // ------------------------------------------------------------------------ 
    // Initialize Monte Carlo class 
    // ------------------------------------------------------------------------ 
    TFIM tfim(spins, bonds, &xfield, &zfield, params["process_id"].as<long>(), beta, params["binsize"].as<long>());
    //------------------------------------------------------------------------ 
    //Run pre-equilibration only if starting from scratch
    //------------------------------------------------------------------------ 
    //if  (params.count("state")){
        cout << endl << "=== Equilibration stage ===" << endl << endl;
        int NoAdjust;
        NoAdjust=0;
        for (int i=0; i!=params["binsize"].as<long>()*params["equil"].as<long>(); i++){
            if (!tfim.AdjustLoopsN()) NoAdjust += 1;
            else                      NoAdjust  = 0;

            if (NoAdjust == 10){
                cout << "  The chain seems to be equilibrated after " << i << " MC updates" << endl;
                cout << "  Time to be productive! " << endl;
                break;
            }        
        }
    //}
    //else{ 
    //    cout << "No pre-equilibration. Loading from: " << params["state"].as<string>() << endl;
    //}
        
    // ------------------------------------------------------------------------ 
    // Run the main loop 
    // ------------------------------------------------------------------------ 
    cout << endl << "=== Measurement stage ===" << endl << endl;
    long binSize = params["binsize"].as<long>();
    for (long i=0; i!=binSize*params["meas"].as<long>(); i++){
        // Perform  a full MC sweep
        
        while (tfim.DiagonalMove()==1)               // diagonal update
              tfim.AdjustM();
              
        tfim.ConstructLinks();                       // linked list and vertex list construction 
        tfim.OffDiagonalMove();                      // cluster update
        tfim.MapStateBack();                         // mapping back the updated state 
        if ((i % binSize) == 0){ 
            tfim.Record();
            tfim.resetMeas();
        }
    }

    delete spins;
    delete bonds;
    return 0;
}
