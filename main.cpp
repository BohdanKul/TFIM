#include <string.h>
#include <boost/format.hpp> 
#include <boost/program_options.hpp>
#include "tfim.h" 

using namespace std;
namespace po = boost::program_options;



string getOptionsList(const vector<string> options) {
    ostringstream optionList;
    copy(options.begin(), options.end() - 1, ostream_iterator<string>(optionList, ", "));
    optionList << options.back();

    return optionList.str();
}



int main(int argc, char *argv[])
{
    vector<string> lLatticeName;
        lLatticeName.push_back("chimera");
        lLatticeName.push_back("rectangle");
    string latticeNames = getOptionsList(lLatticeName);


    po::options_description cmdLineOptions("Command line options"); 
    po::options_description simulationOptions("Simulation options"); 
    po::options_description physicalOptions("Physical parameters options"); 
    po::variables_map params;
    simulationOptions.add_options()
            ("help,h", "produce help message")
            ("state,s",      po::value<string>()->default_value(""),"path to the state file")
            ("process_id,p", po::value<long>()->default_value(0),"random generator seed")
            ("equil,e",      po::value<long>()->default_value(0),"number of equilibration sweeps to take")
            ("meas,m",       po::value<long>(),"number of measurements to take")
            ("binsize",      po::value<long>()->default_value(100),"number of MC sweeps in one measurement (default = 100)")
            ;
    physicalOptions.add_options()
            ("temp,T",       po::value<float>(),"temperature")
            ("beta,B",       po::value<float>(),"inverse temperature")
            ("lattice,L",    po::value<string>()->default_value("rectangle"), boost::str(boost::format("interaction potential type {%s}") % latticeNames).c_str())
            ("width,X",      po::value<int>(),"lattice width")
            ("height,Y",     po::value<int>(),"lattice height")
            ("unity",        po::value<int>(),"lattice unit cell height")
            ("trans,H",      po::value<float>(),"strength of the transverse field")
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


    if  (!(params.count("trans"))){
        cerr << "Error: specify the transverse field (h)" << endl;
        return 1; 
    }
    
    if  (!(params.count("meas"))){
        cerr << "Error: define the number of measurements to take (m)" << endl;
        return 1; 
    }



    if (find(lLatticeName.begin(), lLatticeName.end(), params["lattice"].as<string>()) == lLatticeName.end()) {
       cerr << endl << "ERROR: invalid lattice type" << endl << endl;
       cerr << "Action: set it (L) to one of:" << "\t[" << latticeNames << "]" <<  endl;
       return 1;
    }



    if (params["lattice"].as<string>() == "chimera"){
        if  (!(params.count("unity"))){
            cerr << "Error: define unit cell height (unity)" << endl;
            return 1; 
        }
    }



    if  (!(params.count("width"))){
        cerr << "Error: define lattice (x)" << endl;
        return 1; 
    }
    


    if  (!(params.count("height"))){
        cerr << "Error: define lattice height (y)" << endl;
        return 1; 
    }
 

    // Initialize lattice object to a rectangle or chimera lattice  
    Lattice*   pLattice; 
    if (params["lattice"].as<string>() == "rectangle"){
       pLattice = new Rectangle(params["width"].as<int>(),params["height"].as<int>(),
                                false,                    params["process_id"].as<long>());
    }
    else{
       pLattice = new Chimera(params["width"].as<int>(), params["height"].as<int>(),
                              params["process_id"].as<long>(),   params["unity"].as<int>());
    }


    TFIM tfim(pLattice, params["process_id"].as<long>(),
              beta,  params["trans"].as<float>(), params["binsize"].as<long>());
    
    // Run pre-equilibration only if starting from scratch
    if  (params["state"].as<string>() == ""){
        cout << endl << "Equilibration stage" << endl << endl;
        int NoAdjust;
        NoAdjust=0;
        for (int i=0; i!=params["binsize"].as<long>()*params["equil"].as<long>(); i++){
        
            // Perform  a full MC sweep
            if (tfim.DiagonalMove()==1) tfim.AdjustM();  // diagonal update
            tfim.ConstructLinks();                       // linked list and vertex list construction 
            tfim.OffDiagonalMove();                      // cluster update
            tfim.MapStateBack();                         // mapping back the updated state 
        }
    }
    else{ 
        cout << "No pre-equilibration. Loading from: " << params["state"].as<string>() << endl;
    }
        
    cout << endl << "Measurement stage" << endl << endl;
    for (long i=0; i!=params["binsize"].as<long>()*params["meas"].as<long>(); i++){
        // Perform  a full MC sweep
        
        while (tfim.DiagonalMove()==1)               // diagonal update
              tfim.AdjustM();
              
        tfim.ConstructLinks();                       // linked list and vertex list construction 
        tfim.OffDiagonalMove();                      // cluster update
        tfim.MapStateBack();                         // mapping back the updated state 
        tfim.Measure();
    }

    delete pLattice;
    return 0;
}
