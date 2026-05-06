#include "particle_tree.h"
#include "pion.h"
#include "kt_parser.h"

#include <TROOT.h>
#include <TMath.h>
#include <TString.h>
#include <TH2D.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <chrono>

std::chrono::time_point<std::chrono::high_resolution_clock> progressbar(
    Int_t event, Int_t nEvents, Int_t percentTime,
    std::chrono::time_point<std::chrono::high_resolution_clock> t_prev,
    Bool_t printBar
);
std::string get_current_time();

int main(int argc, char** argv) {
    auto t_start = std::chrono::high_resolution_clock::now();
    
    const char* env_p = std::getenv("HEPQC_ROOT");
    if (!env_p) {
        std::cerr << "Error: HEPQC_ROOT not set!\n";
    }
    std::string projRoot = env_p;

    // --- Parse CLI arguments ---
    std::string inTreeFile = projRoot+"/config/rootFiles.txt";
    Int_t Nmaxevt = 0;
    TString outHistFilename = (
        projRoot+"/data/zphihists_"+get_current_time()+".root"
    );
    Int_t pionCharge = 1;
    Bool_t applyCut = kTRUE;
    Int_t percentInterval {1};
    Bool_t printProgressbar {kFALSE};
    Int_t startFromEvent {0};
    Int_t writeAfterEach {-1};

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [-n|--nevent Nmaxevt=0] "
                      << "[-i|--input "
                      << "inTreeFile=HEPQC_ROOT/config/rootFiles.txt] "
                      << "[-o|--output outHistFilename="
                      << "HEPQC_ROOT/data/zphihists_CURRENTTIME.root] "
                      << "[-c|--charge pionCharge=1] "
                      << "[--cut applyCut=1] "
                      << "[--no-cut applyCut=0] "
                      << "[--percent-interval percentInterval=1] "
                      << "[-p|--progressbar printProgressbar=kFALSE] "
                      << "[-s|--start-from startFromEvent=0] "
                      << "[-w|--write-after writeAfterEach=-1] "
                      << std::endl;
            return -1;
        } else if ((arg == "-n" || arg == "--nevent") && i+1 < argc) {
            Nmaxevt = std::atoi(argv[++i]);
        } else if ((arg == "-i" || arg == "--input") && i+1 < argc) {
            inTreeFile = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i+1 < argc) {
            outHistFilename = argv[++i];
        } else if ((arg == "-c" || arg == "--charge") && i+1 < argc) {
            pionCharge = std::atoi(argv[++i]);
            if (pionCharge == -1) pionCharge = 0;
        } else if (arg == "--percent-interval" && i+1 < argc) {
            percentInterval = std::atoi(argv[++i]);
        } else if (arg == "-p" || arg == "--progressbar") {
            printProgressbar = kTRUE;
        } else if (arg == "-s" || arg == "--start-from") {
            startFromEvent = std::atoi(argv[++i]);
        } else if (arg == "-w" || arg == "--write-after") {
            writeAfterEach = std::atoi(argv[++i]);
        }
    }
    // --- Parsing finished ---

    // Create output file and check wether it is valid
    TFile* outHistFile = TFile::Open(outHistFilename, "RECREATE");

    if (!outHistFile || outHistFile->IsZombie()) {
        std::cerr << "FATAL: Could not open output file " << outHistFilename
                  << " for writing\n";
        return 1;
    }

    // Initializing KTParser
    std::string chargeStr = pionCharge > 0 ? "plus" : "minus";
    KTParser<TH1D> KThists(
        projRoot+"/config/KTbins.txt",
        projRoot+"/config/KTcuts_"+chargeStr+".txt",
        "Q", "(Q); Q [GeV]; counts",
        0, 30, 10, -15.0, 15.0, 10
    );

    KThists.initialize_ABhists(1000, 0, 0.5);

    // Initializing particle_tree
    std::cerr << '\n';
    particle_tree ptree(inTreeFile.c_str());
    if(ptree.fChain) std::cerr << "Tree initialized" << std::endl;
    else { std::cerr << "No tree found." << std::endl;}

    // Preparing for event loop
    int Nevents = ptree.fChain->GetEntries()-startFromEvent;
    
    if(Nmaxevt>0 && Nmaxevt<Nevents) {Nevents=Nmaxevt;}
    std::cerr << "\nWill run on " << Nevents << " events (out of "
    << ptree.fChain->GetEntries()  << ").\n"
    << "Will start from event " << startFromEvent << '\n';
    
    std::cerr << "Will analyze " << (pionCharge==0 ? "negative" : "positive")
    << " pions.\n";
    
    auto t_init = std::chrono::high_resolution_clock::now();
    auto t_prev = t_init;
    std::cerr << "\nPre initialization took "
    << std::chrono::duration_cast<std::chrono::seconds>(
        t_init-t_start
    ).count() << " seconds\n";
    
    // Running event loop
    for(int i = startFromEvent; i < startFromEvent+Nevents; i++) {
        
        t_prev = progressbar(
            i-startFromEvent+1, Nevents, percentInterval,
            t_prev, printProgressbar
        );

        if (writeAfterEach>0 && (i-startFromEvent) % writeAfterEach == 0) {
            std::cerr << "Analyzed " << i-startFromEvent << " events. ";
            KThists.write_hists(outHistFile);
        }

        ptree.GetEntry(i);        

        KThists.clear_current_event();
        KThists.current_event.centrality = ptree.Centrality;
        KThists.current_event.Zvertex = ptree.Zvertex;

        // Apply event cut
        if (KThists.get_pool_index() == -1) {continue;}

        for(int ipart = 0; ipart < ptree.Ntracks; ipart++) {
            // Search for pions with specific charge
            if(ptree.isPi[ipart]==0 && ptree.ch[ipart]==pionCharge) {
                Pion pi(
                    ptree.px[ipart], ptree.py[ipart], ptree.pz[ipart],
                    ptree.detz[ipart], ptree.detp[ipart], ptree.ch[ipart]
                );
                // Apply particle cuts
                if(
                    TMath::Abs(pi.eta) < 0.35
                    && pi.p_abs > 0.2 && pi.p_abs < 1.0
                ) { KThists.add_pion(pi); }
            }
        }

        if (KThists.current_event.pions.empty()) {
            KThists.Nemptypionvects++;
            continue;
        }
 
        KThists.fill_Ahists(applyCut);
        KThists.fill_Bhists(applyCut);
    }
    
    auto t_endanal = std::chrono::high_resolution_clock::now();
    std::cerr << "Analyzing events took "
              << std::chrono::duration_cast<std::chrono::seconds>(
                t_endanal-t_init
              ).count() << " seconds\n";
              
    std::cerr << "\nNumber of empty pion vectors: "
              << KThists.Nemptypionvects << '\n';
    std::cerr << "\nNumber of rejected pairs of pions: "
              << KThists.NrejectedPairs << '\n';
    
    KThists.create_Chists();
    KThists.write_hists(outHistFile);

    outHistFile->Close();
    delete outHistFile;

    auto t_write = std::chrono::high_resolution_clock::now();
    std::cerr << "\nWriting histograms took "
              << std::chrono::duration_cast<std::chrono::seconds>(
                t_write-t_endanal
              ).count() << " seconds\n";
    std::cerr << "\nWhole analysis took "
              << std::chrono::duration_cast<std::chrono::seconds>(
                t_write-t_start
              ).count() << " seconds\n";

    return 0;
}