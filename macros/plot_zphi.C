#include "kt_parser.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TString.h>

#include <iostream>
#include <vector>

void plot_hist(
    TCanvas* c, TH2D* hist, TString filename,
    Double_t histMin=-1, Double_t histMax=-1
);

Double_t calc_hist_region_avg(TH2D* hist, Int_t numXBins, Int_t numYBins);

void plot_zphi(
    TString hist_file, TString addToFilename, TString filetype=".png"
) {
    KTParser<TH2D> kt_hists(
        "../config/KTbins.txt", "../data/"+hist_file,
        "zphi", "(#Delta z, #Delta #phi); #Delta z [cm]; #Delta #phi"
    );

    kt_hists.scale_ABhists();
    kt_hists.merge_hists(0,4);
    kt_hists.merge_hists(1,5);
    kt_hists.merge_hists(2,6);
    kt_hists.merge_hists(3,7);
    kt_hists.merge_hist_bins(10,10);

    for(Int_t i = 0; i < kt_hists.num_KT_bins; i++) {
        TString title = Form(
            "C(#Delta z, #Delta #phi), %.2f < K_{T} [GeV] < %.2f",
            kt_hists.KT_bin_edges.at(i), kt_hists.KT_bin_edges.at(i+1)
        );
        kt_hists.Chists.at(i)->SetTitle(title);
    }

    std::vector<Double_t> scaleFactors(kt_hists.num_KT_bins, 0);
    for(Int_t i = 0; i < kt_hists.num_KT_bins; i++) {
        scaleFactors.at(i) = calc_hist_region_avg(
            kt_hists.Chists.at(i), 10, 10
        );
    }
    kt_hists.scale_hists(scaleFactors, 'C');

    kt_hists.print_kt_bin_edges();

    TCanvas *c = new TCanvas("c","c",1200,900);
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    for(Int_t i = 0; i < kt_hists.num_KT_bins; i++) {
        std::cerr << "A integral = " << kt_hists.Ahists.at(i)->Integral("width")
                  << '\n';
        plot_hist(
            c, kt_hists.Ahists.at(i),
            "../figures/Ahists/Azphi_"+addToFilename+"_"+
            std::to_string(i)+filetype
        );
        std::cerr << "B integral = " << kt_hists.Bhists.at(i)->Integral("width")
                  << '\n';
        plot_hist(
            c, kt_hists.Bhists.at(i),
            "../figures/Bhists/Bzphi_"+addToFilename+"_"
            +std::to_string(i)+filetype
        );
        std::cerr << "C integral = " << kt_hists.Chists.at(i)->Integral("width")
                  << '\n';
        plot_hist(
            c, kt_hists.Chists.at(i),
            "../figures/Chists/Czphi_"+addToFilename+"_"
            +std::to_string(i)+filetype
            , 0.0, 1.6
        );
    }
}

void plot_hist(
    TCanvas* c, TH2D* hist, TString filename,
    Double_t histMin, Double_t histMax
) {
    c->cd();
    if (histMin<0) {
        hist->SetMinimum(histMin);
    } else {
        hist->SetMinimum(histMin);
    }
    if (histMax<0) {
        hist->SetMinimum(histMax);
    } else {
        hist->SetMaximum(histMax);
    }
    hist->Draw("COLZ");
    c->Print(filename);
    c->Clear();
}

Double_t calc_hist_region_avg(TH2D* hist, Int_t numXBins, Int_t numYBins) {
    Double_t integral = hist->Integral(
        hist->GetNbinsX()-numXBins, hist->GetNbinsX(),
        hist->GetNbinsY()-numYBins, hist->GetNbinsY()
    );

    Int_t num_bins = (numXBins+1)*(numYBins+1);
    return integral/num_bins;
}