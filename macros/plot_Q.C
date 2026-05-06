#include "kt_parser.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>

#include <iostream>
#include <vector>

void plot_ABhist(
    TCanvas* c, TH1D* Ahist, TH1D* Bhist, TString filename,
    Double_t histMin=-1111, Double_t histMax=-1111,
    Double_t Qmin=0, Double_t Qmax=-1
);
void plot_Chist(
    TCanvas* c, TH1D* Chist, TString filename,
    Double_t histMin=-1, Double_t histMax=-1,
    Double_t QMin=0, Double_t QMax=-1
);

void plot_Q(TString hist_file, TString addToFilename) {
    KTParser<TH1D> kt_hists(
        "../config/KTbins.txt", "../data/"+hist_file,
        "Q", "(Q); Q [GeV]; Counts"
    );

    kt_hists.merge_hist_bins(4);
    kt_hists.scale_ABhists();

    for(Int_t i = 0; i < kt_hists.num_KT_bins; i++) {
        TString title = Form(
            "A(Q) and B(Q), %.2f < K_{T} [GeV] < %.2f",
            kt_hists.KT_bin_edges.at(i), kt_hists.KT_bin_edges.at(i+1)
        );
        kt_hists.Ahists.at(i)->SetTitle(title);
        title = Form(
            "C(Q), %.2f < K_{T} [GeV] < %.2f",
            kt_hists.KT_bin_edges.at(i), kt_hists.KT_bin_edges.at(i+1)
        );
        kt_hists.Chists.at(i)->SetTitle(title);
    }

    kt_hists.print_kt_bin_edges();

    TCanvas *c = new TCanvas("c","c",1200,900);
    
    gStyle->SetOptStat(0);

    for(Int_t i = 0; i < kt_hists.num_KT_bins; i++) {
        plot_ABhist(
            c, kt_hists.Ahists.at(i), kt_hists.Bhists.at(i),
            "../figures/ABhists/ABQ_"+addToFilename+"_"+
            std::to_string(i)+".pdf",
            -1111, -1111, 0.0, 0.5
        );
        plot_Chist(
            c, kt_hists.Chists.at(i),
            "../figures/CQhists/CQ_"+addToFilename+"_"
            +std::to_string(i)+".pdf",
            0, 2, 0.0, 0.5
        );
    }
}

void plot_ABhist(
    TCanvas* c, TH1D* Ahist, TH1D* Bhist, TString filename,
    Double_t histMin, Double_t histMax,
    Double_t Qmin, Double_t Qmax
) {
    c->cd();
    c->SetLogx();
    c->SetLogy();
    Ahist->SetMinimum(histMin);
    Ahist->SetMaximum(histMax);
    if (Qmin > 0 || Qmax > 0) {
        Ahist->GetXaxis()->SetRangeUser(Qmin,Qmax);
    }
    
    Ahist->SetLineColor(kBlue);
    Ahist->Draw();
    Bhist->SetLineColor(kRed);
    Bhist->Draw("SAME");
    
    TLegend *leg = new TLegend(0.1, 0.7, 0.2, 0.9); // (x1, y1, x2, y2) in NDC coordinates    
    leg->AddEntry(Ahist, "A(Q)", "l"); // "l" means use the line style
    leg->AddEntry(Bhist, "B(Q)", "l");
    leg->Draw();
    
    c->Print(filename);
    c->SetLogx(0);
    c->SetLogy(0);
    c->Clear();
}

void plot_Chist(
    TCanvas* c, TH1D* Chist, TString filename,
    Double_t histMin, Double_t histMax,
    Double_t Qmin, Double_t Qmax
) {
    c->cd();
    Chist->SetMinimum(histMin);
    Chist->SetMaximum(histMax);
    if (Qmin > 0 || Qmax > 0) {
        Chist->GetXaxis()->SetRangeUser(Qmin,Qmax);
    }

    Chist->Draw();

    c->Print(filename);
    c->Clear();
}
