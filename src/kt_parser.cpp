#include "kt_parser.h"

#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <type_traits>

void KTParserBase::parse_config_file() {
    std::ifstream config_file(KTconfig_file);

    if(!config_file) {
        std::cerr << "\nFile not found: " << KTconfig_file << '\n';
        return;
    }

    std::cerr << "\nParsing KT config file " << KTconfig_file << "...\n";

    std::string config_line {};
    while(getline(config_file, config_line)) {
        std::stringstream ss(config_line);
        std::string param {};
        ss >> param;
        if (param=="#") {
            continue;
        } else if (param=="NumBins") {
            ss >> num_KT_bins;
        } else if (param=="KTMin") {
            ss >> KT_min;
        } else if (param=="KTMax") {
            ss >> KT_max;
        } else if (param=="BinEdges") {
            Float_t bin_edge {};
            while(!ss.eof()) {
                ss >> bin_edge;
                KT_bin_edges.push_back(bin_edge);
            }
            std::sort(KT_bin_edges.begin(), KT_bin_edges.end());
        }
    }

    config_file.close();

    if (KT_min > KT_max) {
        Float_t temp {KT_min};
        KT_min = KT_max;
        KT_max = temp;
    }

    if ((KT_min < 0 || KT_max < 0 || num_KT_bins < 0) && KT_bin_edges.empty()) {
        std::cerr << "Config file is incomplete. Can't define KT bins.\n";
        return;
    } else if (!KT_bin_edges.empty()) {
        std::cerr << "Using provided KT bin edges.\n";
        KT_min = KT_bin_edges.at(0);
        KT_max = KT_bin_edges.back();
        num_KT_bins = KT_bin_edges.size();
    } else {
        std::cerr << "Creating " << num_KT_bins << " bins between " << KT_min
                  << " GeV and " << KT_max << " GeV.\n";
        for (Int_t i = 0; i < num_KT_bins+1; i++) {
            KT_bin_edges.push_back(KT_min+i*(KT_max-KT_min)/num_KT_bins);
        }
    }

    std::cerr << "Parsing finished. The bin edges are [GeV]:\n";
    for(double edge: KT_bin_edges) {
        std::cerr << edge << '\t';
    }
    std::cerr << '\n';
}

Int_t KTParserBase::get_KT_index(Float_t KT) {
    Int_t idx {-1};
    for(double bin_edge: KT_bin_edges) {
        if (KT > bin_edge) {idx++;}
        else {break;}
    }
    if (idx == num_KT_bins) {idx = -1;}
    return idx;
}

template <typename T>
void KTParser<T>::add_pion(const Pion& pi) {
    current_event.add_pion(pi);
}

template <typename T>
Int_t KTParser<T>::get_pool_index() {
    Int_t pool_idx = pool.get_index(
        current_event.centrality, current_event.Zvertex
    );
    return pool_idx;
}

template <typename T>
void KTParser<T>::clear_current_event() {
    current_event.clear();
}

template <typename T>
void KTParser<T>::print_num_pions() {
    std::cerr << "Number of pions in the current event:"
              << current_event.pions.size() << '\n';
}

template <typename T>
void KTParser<T>::print_kt_bin_edges() {
    std::cerr << "K_T bin edges [GeV]:\n";
    for(double edge: KT_bin_edges) {
        std::cerr << edge << '\t';
    }
    std::cerr << '\n';
}

template <typename T>
void KTParser<T>::initialize_ABhists(
    Int_t num_bins_x, Double_t x_min, Double_t x_max,
    Int_t num_bins_y, Double_t y_min, Double_t y_max
) {
    std::cerr << "KTParser<T>::initialize_ABhists(): No implementation for type "
              << typeid(T).name() << '\n';
}

template <>
void KTParser<TH1D>::initialize_ABhists(
    Int_t num_bins_x, Double_t x_min, Double_t x_max,
    Int_t num_bins_y, Double_t y_min, Double_t y_max
) {
    std::cerr << "\nInitializing A and B histograsm...\n";

    for (Int_t i = 0; i < num_KT_bins; i++) {
        TH1D* Ahist = new TH1D(
            "A"+name+"_"+i, "A"+axis_titles, num_bins_x, x_min, x_max
        );
        Ahist->Sumw2();
        Ahist->SetDirectory(nullptr);
        Ahists.at(i) = std::move(Ahist);
        
        TH1D* Bhist = new TH1D(
            "B"+name+"_"+i, "B"+axis_titles, num_bins_x, x_min, x_max
        );
        Bhist->Sumw2();
        Bhist->SetDirectory(nullptr);
        Bhists.at(i) = std::move(Bhist);
    }
    
    std::cerr << Ahists.size() << " pieces of A and "
              << Bhists.size() << " pieces of B histograms initialized.\n";
}

template <>
void KTParser<TH2D>::initialize_ABhists(
    Int_t num_bins_x, Double_t x_min, Double_t x_max,
    Int_t num_bins_y, Double_t y_min, Double_t y_max
) {
    std::cerr << "\nInitializing A and B histograsm...\n";

    for (Int_t i = 0; i < num_KT_bins; i++) {
        TH2D* Ahist = new TH2D(
            "A"+name+"_"+i, "A"+axis_titles,
            num_bins_x, x_min, x_max,
            num_bins_y, y_min, y_max
        );
        Ahist->Sumw2();
        Ahist->SetDirectory(nullptr);
        Ahists.at(i) = std::move(Ahist);

        TH2D* Bhist = new TH2D(
            "B"+name+"_"+i, "B"+axis_titles,
            num_bins_x, x_min, x_max,
            num_bins_y, y_min, y_max
        );
        Bhist->Sumw2();
        Bhist->SetDirectory(nullptr);
        Bhists.at(i) = std::move(Bhist);
    }

    std::cerr << Ahists.size() << " pieces of A and "
              << Bhists.size() << " pieces of B histograms initialized.\n";
}

// -- Helper --
Float_t calc_KT(const Pion& p1, const Pion& p2) {
    Float_t KT = TMath::Sqrt(
        (p1.px+p2.px) * (p1.px+p2.px) + (p1.py+p2.py) * (p1.py+p2.py)
    ) / 2;

    return KT;
}

Float_t calc_moment_diff(const Pion& p1, const Pion& p2) {
    Float_t q_long_LCMS2 = 4*(p1.pz*p2.E-p2.pz*p1.E)*(p1.pz*p2.E-p2.pz*p1.E) / (
        (p1.E+p2.E)*(p1.E+p2.E) - (p1.pz+p2.pz)*(p1.pz+p2.pz)
    );

    Float_t Q = TMath::Sqrt(
        (p1.px-p2.px)*(p1.px-p2.px) + (p1.py-p2.py)*(p1.py-p2.py) + q_long_LCMS2
    );

    return Q;
}

Float_t calc_dz(const Pion& p1, const Pion& p2) {
    return TMath::Abs(p1.det_z - p2.det_z);
}

Float_t calc_dphi(const Pion& p1, const Pion& p2) {
    return TMath::Abs(p1.det_phi - p2.det_phi);
}
// --...---

template <>
void KTParser<TH1D>::fill_Ahists(Bool_t cut) {
    if(current_event.pions.size() < 2) return;
    
    for(
        std::size_t ipi = 0;
        ipi < current_event.pions.size();
        ipi++
    ) {
        for(
            std::size_t jpi = ipi+1;
            jpi < current_event.pions.size();
            jpi++
        ) {
            const auto& p1 = current_event.pions.at(ipi);
            const auto& p2 = current_event.pions.at(jpi);
            Float_t dz = calc_dz(p1, p2);
            Float_t dphi = calc_dphi(p1, p2);
            Float_t KT = calc_KT(p1, p2);
            Int_t ktbin = get_KT_index(KT);
            if (cut && ZPhiCuts.isRejected(ktbin, dz, dphi)) {
                NrejectedPairs++;
                continue;
            }
            
            Float_t Q = calc_moment_diff(p1, p2);
            Ahists.at(ktbin)->Fill(Q);
        }
    }
}

template <>
void KTParser<TH2D>::fill_Ahists(Bool_t cut) {
    if(current_event.pions.size() < 2) return;
    
    for(
        std::size_t ipi = 0;
        ipi < current_event.pions.size();
        ipi++
    ) {
        for(
            std::size_t jpi = ipi+1;
            jpi < current_event.pions.size();
            jpi++
        ) {
            const auto& p1 = current_event.pions.at(ipi);
            const auto& p2 = current_event.pions.at(jpi);
            Float_t dz = calc_dz(p1, p2);
            Float_t dphi = calc_dphi(p1, p2);
            Float_t KT = calc_KT(p1, p2);
            Int_t ktbin = get_KT_index(KT);
            if (cut && ZPhiCuts.isRejected(ktbin, dz, dphi)) {
                NrejectedPairs++;
                continue;
            }

            Ahists.at(ktbin)->Fill(dz, dphi);
        }
    }
}

template <>
void KTParser<TH1D>::fill_Bhists(Bool_t cut) {
    Int_t pool_idx = get_pool_index();
    if(pool.pool.at(pool_idx).size()==5) {
        for(const auto& p1: current_event.pions) {
            for(const auto& bkg_evt: pool.pool.at(pool_idx)) {
                for(const auto& p2: bkg_evt.pions) {
                    Float_t dz = calc_dz(p1, p2);
                    Float_t dphi = calc_dphi(p1, p2);
                    Float_t KT = calc_KT(p1, p2);
                    Int_t ktbin = get_KT_index(KT);
                    if (cut && ZPhiCuts.isRejected(ktbin, dz, dphi)) {
                        NrejectedPairs++;
                        continue;
                    }
                    
                    Float_t Q = calc_moment_diff(p1, p2);
                    Bhists.at(ktbin)->Fill(Q);
                }
            }
        }

        pool.pool.at(pool_idx).pop_front();
        pool.pool.at(pool_idx).push_back(current_event);
    } else {
        pool.pool.at(pool_idx).emplace_back(std::move(current_event));
    }
}

template <>
void KTParser<TH2D>::fill_Bhists(Bool_t cut) {
    Int_t pool_idx = get_pool_index();
    if(current_event.pions.empty()) return;

    if(pool.pool.at(pool_idx).size()==5) {
        for(const auto& p1: current_event.pions) {
            for(const auto& bkg_evt: pool.pool.at(pool_idx)) {
                for(const auto& p2: bkg_evt.pions) {
                    Float_t dz = calc_dz(p1, p2);
                    Float_t dphi = calc_dphi(p1, p2);
                    Float_t KT = calc_KT(p1, p2);
                    Int_t ktbin = get_KT_index(KT);
                    if (cut && ZPhiCuts.isRejected(ktbin, dz, dphi)) {
                        NrejectedPairs++;
                        continue;
                    }

                    Bhists.at(ktbin)->Fill(dz, dphi);
                }
            }
        }

        pool.pool.at(pool_idx).pop_front();
        pool.pool.at(pool_idx).push_back(current_event);
    } else {
        pool.pool.at(pool_idx).emplace_back(std::move(current_event));
    }
}

template <typename T>
void KTParser<T>::scale_ABhists() {
    std::cerr << "KTParser<T>::scale_ABhists(): No definition for type"
              << typeid(T).name() << '\n';
}

template <>
void KTParser<TH1D>::scale_ABhists() {
    // Scaling only with the integral of non-femtoscopic region
    Int_t start_bin = Bhists.at(0)->FindBin(0.15);
    Int_t end_bin = Bhists.at(0)->GetNbinsX();
    for (std::size_t i = 0; i < Ahists.size(); i++) {
        Bhists.at(i)->Scale(
            Ahists.at(i)->Integral(start_bin, end_bin, "width")
            / Bhists.at(i)->Integral(start_bin, end_bin, "width")
        );
    }
}

template <>
void KTParser<TH2D>::scale_ABhists() {    
    for(std::size_t i = 0; i < Ahists.size(); i++) {
        Double_t aint = Ahists.at(i)->Integral("width");
        Double_t bint = Bhists.at(i)->Integral("width");
        if (bint == 0.0) {
            std::cerr << "Warning: B histogram " << i << " has zero integral. "
                      << "Skipping scale.\n";
            continue;
        }
        Bhists.at(i)->Scale( aint / bint );
    }
}

template <typename T>
void KTParser<T>::scale_hists(Double_t scaleFactor, char which) {    
    std::vector<T*>* hists {nullptr};
    switch (which)
    {
    case 'A':
        hists = &Ahists;
        break;
    case 'B':
        hists = &Bhists;
        break;
    case 'C':
        hists = &Chists;
        break;
    
    default:
        break;
    }
    if (!hists) {
        std::cerr << "Choose histograms A,B or C.\n";
        return;
    }
    for(Int_t ktbin = 0; ktbin < num_KT_bins; ktbin++) {
        hists->at(ktbin)->Scale( 1.0 / scaleFactor );
    }
}

template <typename T>
void KTParser<T>::scale_hists(
    const std::vector<Double_t>& scaleFactors, char which
) {    
    std::vector<T*>* hists {nullptr};
    switch (which)
    {
    case 'A':
        hists = &Ahists;
        break;
    case 'B':
        hists = &Bhists;
        break;
    case 'C':
        hists = &Chists;
        break;
    
    default:
        break;
    }
    if (!hists) {
        std::cerr << "Choose histograms A,B or C.\n";
        return;
    }
    for(Int_t ktbin = 0; ktbin < num_KT_bins; ktbin++) {
        hists->at(ktbin)->Scale( 1.0 / scaleFactors.at(ktbin) );
    }
}

template <typename T>
void KTParser<T>::scale_Chists() {
    for(std::size_t i = 0; i < Chists.size(); i++) {
        Double_t cint = Chists.at(i)->Integral("width");
        if (cint == 0.0) {
            std::cerr << "Warning: C histogram " << i << " has zero integral. "
                      << "Skipping scale.\n";
            continue;
        }
        Chists.at(i)->Scale(1.0 / cint);
    }
}

template <typename T>
void KTParser<T>::create_Chist(Int_t kt_idx) {
    T* Chist = (T*)Ahists.at(kt_idx)->Clone("C"+name+"_"+kt_idx);
    Chist->Divide(Bhists.at(kt_idx));
    Chist->SetTitle("C"+axis_titles);
    Chist->SetDirectory(nullptr);
    Chists.at(kt_idx) = std::move(Chist);
}

template <typename T>
void KTParser<T>::create_Chists() {
    std::cerr << "\nCreating C histograms...\n";
    for (std::size_t i = 0; i < Ahists.size(); i++) {
        create_Chist(i);
    }
    std::cerr << Chists.size() << " pieces of C histograms created\n";
}

template <typename T>
void KTParser<T>::write_hists(TFile* out_file) {

    if (!out_file || out_file->IsZombie() || !out_file->IsWritable()) {
        std::cerr << "Error: Provided TFile is invalid or not writable!"
                  << std::endl;
        return;
    }

    out_file->cd();

    std::cerr << "\nWriting histograms to " << out_file->GetName() << "...\n";
    for (auto h: Ahists) { if (h) h->Write(nullptr, TObject::kWriteDelete); }
    for (auto h: Bhists) { if (h) h->Write(nullptr, TObject::kWriteDelete); }
    for (auto h: Chists) { if (h) h->Write(nullptr, TObject::kWriteDelete); }

    out_file->SaveSelf();
}

template <typename T>
void KTParser<T>::load_hists(TString hist_file, TString hist_name) {
    TFile* file = new TFile(hist_file, "READ");

    std::cerr << "\nLoading histograms from " << hist_file << "...\n";
    for(Int_t i = 0; i < num_KT_bins; i++) {
        Ahists.push_back((T*)file->Get("A"+hist_name+"_"+i));
        Bhists.push_back((T*)file->Get("B"+hist_name+"_"+i));
        Chists.push_back((T*)file->Get("C"+hist_name+"_"+i));
    }
    std::cerr << "# A histograms loaded: " << Ahists.size() << '\n'
              << "# B histograms loaded: " << Bhists.size() << '\n'
              << "# C histograms loaded: " << Chists.size() << '\n';
}

// Drafted by copilot, but examined thoroughly
// Rewritten since to handle non-single merges as well
template <typename T>
void KTParser<T>::merge_hists(Int_t idx1, Int_t idx2) {
    // Check if indices are valid
    if (idx1 < 0 || idx1 >= num_KT_bins || idx2 < 0 || idx2 >= num_KT_bins) {
        std::cerr << "Invalid indices for merging: "
                  << idx1 << ", " << idx2 << '\n';
        return;
    }
    
    std::cerr << "Merging histograms from "
              << idx1 << " to " << idx2 << "...\n";
    
    // Ensure idx1 < idx2
    if (idx1 > idx2) {
        std::swap(idx1, idx2);
    }

    for(std::size_t i = idx1; i < idx2; i++) {
        // Merge histograms
        if (Ahists.at(idx1) && Ahists.at(idx1+1)) {
            Ahists.at(idx1)->Add(Ahists.at(idx1+1));
            delete Ahists.at(idx1+1);
        }
        if (Bhists.at(idx1) && Bhists.at(idx1+1)) {
            Bhists.at(idx1)->Add(Bhists.at(idx1+1));
            delete Bhists.at(idx1+1);
        }
        // Remove the merged elements from vectors
        Ahists.erase(Ahists.begin() + idx1+1);
        Bhists.erase(Bhists.begin() + idx1+1);
        
        create_Chist(idx1);
        delete Chists.at(idx1+1);
        Chists.erase(Chists.begin() + idx1+1);
        
        // Update KT_bin_edges: remove the edge between the merged bins
        KT_bin_edges.erase(KT_bin_edges.begin() + idx1+1);
        // Update num_KT_bins
        num_KT_bins--;
    }
    std::cerr << "Merging finished. Number of K_T bins remaining: "
              << num_KT_bins << '\n';
}

template <>
void KTParser<TH1D>::merge_hist_bins(Int_t nMergeX, Int_t nMergeY) {
    Int_t nBins = Ahists.at(0)->GetNbinsX();
    if (nBins % nMergeX != 0) {
        std::cerr << "Can not merge every " << nMergeX
                  << " bins from " << nBins << " bins\n";
        return;
    }

    for(auto h: Ahists) {if (h) h->Rebin(nMergeX);}
    for(auto h: Bhists) {if (h) h->Rebin(nMergeX);}

    for(Int_t ktbin = 0; ktbin < num_KT_bins; ktbin++) {
        create_Chist(ktbin);
    }
}

template <>
void KTParser<TH2D>::merge_hist_bins(Int_t nMergeX, Int_t nMergeY) {
    Int_t nBinsX = Ahists.at(0)->GetNbinsX();
    Int_t nBinsY = Ahists.at(0)->GetNbinsX();
    if (nBinsX % nMergeX != 0) {
        std::cerr << "Can not merge every " << nMergeX
                  << " bins from " << nBinsX << " bins\n";
        return;
    }
    if (nBinsY % nMergeY != 0) {
        std::cerr << "Can not merge every " << nMergeY
                  << " bins from " << nBinsY << " bins\n";
        return;
    }

    for(auto h: Ahists) {if (h) h->Rebin2D(nMergeX, nMergeY);}
    for(auto h: Bhists) {if (h) h->Rebin2D(nMergeX, nMergeY);}

    for(Int_t ktbin = 0; ktbin < num_KT_bins; ktbin++) {
        create_Chist(ktbin);
    }
}

template class KTParser<TH1D>;
template class KTParser<TH2D>;