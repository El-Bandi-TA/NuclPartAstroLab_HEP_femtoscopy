#include "kt_pair_cuts.h"

#include <TROOT.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>

void KTPairCuts::parse_config_file(std::string pairConfigFile) {
    std::ifstream config_file(pairConfigFile);

    if(!config_file) {
        std::cerr << "\nFile not found: " << pairConfigFile << '\n';
        return;
    }

    std::cerr << "\nParsing KT config file " << pairConfigFile << "...\n";

    std::string config_line {};
    Bool_t parseCutEdges {kFALSE};
    Bool_t parseCutRects {kFALSE};
    Int_t cutBin {0};
    while(getline(config_file, config_line)) {
        std::stringstream ss(config_line);
        
        if (parseCutEdges) {
            Float_t cut_edge {};
            while (!ss.eof()) {
                ss >> cut_edge;
                KT_cut_edges.push_back(cut_edge);
            }
            std::sort(KT_cut_edges.begin(), KT_cut_edges.end());
            parseCutEdges = kFALSE;
            
            cuts.reserve(KT_cut_edges.size()-1);
            cuts.resize(KT_cut_edges.size()-1);
            continue;
        }
        
        if (parseCutRects) {
            Float_t z1, z2, p1, p2;
            while (!ss.eof()) {
                ss >> z1 >> z2 >> p1 >> p2;
                cuts[cutBin].emplace_back(z1,z2,p1,p2);
            }
            cutBin++;
            continue;
        }

        std::string param {};
        ss >> param;
        if (param=="#") continue;
        else if (param=="KTCutBinInfo") parseCutEdges=kTRUE;
        else if (param=="ZPhiCuts") parseCutRects=kTRUE;
    }

    config_file.close();

    std::cerr << "Parsing finished\n";
}

void KTPairCuts::initialize_bin_map(std::vector<Float_t> KT_bin_edges) {
    if (
        KT_bin_edges.at(0)!=KT_cut_edges.at(0)
        || KT_bin_edges.back() != KT_cut_edges.back()
    ) {
        std::cerr << "Cant initialize bin_map: Cuts defined between "
                  << KT_cut_edges.at(0) << " GeV and " << KT_cut_edges.back()
                  << "GeV, but K_T bins are between " << KT_bin_edges.at(0)
                  << "GeV and " << KT_bin_edges.back() << "GeV.";
        return;
    }
    
    Int_t cutBin = 0;
    for(Int_t ktbin = 1; ktbin < KT_bin_edges.size(); ktbin++) {
        if (KT_bin_edges.at(ktbin)<=KT_cut_edges.at(cutBin+1)) {
            bin_map.push_back(cutBin);
        } else {
            bin_map.push_back(++cutBin);
        }
    }
}

void KTPairCuts::print_cuts(std::vector<Float_t> KT_bin_edges) {
    std::cerr << "The cuts are:\n"
              << "KT\tDelta z\tDelta Phi\n";
    for(Int_t ktbin = 0; ktbin < KT_bin_edges.size()-1; ktbin++) {
        for(
            Int_t icut = 0;
            icut < cuts.at(bin_map.at(ktbin)).size();
            icut++
        ) {
            std::cerr << std::setw(2);
            std::cerr << KT_bin_edges.at(ktbin)
                      << "-"
                      << KT_bin_edges.at(ktbin+1)
                      << "GeV\t"
                      << cuts.at(bin_map.at(ktbin)).at(icut).dzmin
                      << "-"
                      << cuts.at(bin_map.at(ktbin)).at(icut).dzmax
                      << " cm\t"
                      << cuts.at(bin_map.at(ktbin)).at(icut).dphimin
                      << "-"
                      << cuts.at(bin_map.at(ktbin)).at(icut).dphimax
                      << "\n";
        } 
    }
    std::cerr << '\n';
}