#include "kt_parser.h"
#include "CoulCorrCalc.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TLine.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

CoulCorrCalc *cccinstance;

void plot_ABhist(
    TCanvas* c, TH1D* Ahist, TH1D* Bhist, TString filename,
    Double_t histMin=-1111, Double_t histMax=-1111,
    Double_t Qmin=0, Double_t Qmax=-1
);
void plot_Chist(
    TCanvas* c, TH1D* Chist, TString filename,
    Double_t histMin=-1, Double_t histMax=-1,
    Double_t QMin=0, Double_t QMax=-1,
    std::vector<Double_t> fitpars = {}
);
void plot_Levyfit(
    TCanvas* c, TLegend* leg, Double_t Qmin, Double_t Qmax,
    std::vector<Double_t> fitpars
);
void plot_fitRange(
    TCanvas* c, TLegend* leg, Double_t Qmin, Double_t Qmax
);

void plot_Q(
    TString hist_file, TString addToFilename,
    std::string fitparam_file="", TString filetype=".png"
) {
    KTParser<TH1D> kt_hists(
        "../config/KTbins.txt", "../data/"+hist_file,
        "Q", "(Q); Q [GeV]; Counts"
    );

    kt_hists.merge_hist_bins(4);
    // kt_hists.merge_hists(0,3);
    kt_hists.scale_ABhists();
    kt_hists.create_Chists();

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

    // for(Int_t i = kt_hists.num_KT_bins-1; i >= 0 ; i--) {
    for(Int_t i = 0; i < kt_hists.num_KT_bins; i++) {
        // Extract fit parameters
        std::vector<double> fitpars {};
        std::string fitstatus {"no status"};
        if (fitparam_file!="") {
            double N, lambda, R, alpha, Qfitmin, Qfitmax;
            std::ifstream fitparams("../data/"+fitparam_file);
            if (!fitparams.is_open()) {
                std::cerr << "Could not open fit parameter file " << fitparam_file
                          << ". Omitting plot of fitted function.\n";
                fitparam_file="";
            }
            std::string line {};
            double KTmin {-1};
            double KTmax {-1};
            while(getline(fitparams,line)) {
                if(line[0]=='#') continue;
                std::stringstream ss {line};
                ss >>  KTmin >> KTmax;
                if (
                    (std::abs(KTmin-kt_hists.KT_bin_edges.at(i))<1e-5)
                    && (std::abs(KTmax-kt_hists.KT_bin_edges.at(i+1))<1e-5)
                ) {
                    double buffer;
                    ss >> Qfitmin >> Qfitmax >> buffer
                       >> N >> buffer >> buffer >> buffer
                       >> lambda >> buffer >> buffer >> buffer
                       >> R >> buffer >> buffer >> buffer
                       >> alpha >> buffer >> buffer >> buffer
                       >> buffer >> buffer >> buffer >> fitstatus;
                    if (fitstatus != "converged") break;
                    fitpars.push_back(N);
                    fitpars.push_back(lambda);
                    fitpars.push_back(R);
                    fitpars.push_back(alpha);
                    fitpars.push_back(Qfitmin);
                    fitpars.push_back(Qfitmax);
                    break;
                }
            }
            // Initialize cccinstance
            cccinstance = new CoulCorrCalc();
        }
        
        plot_ABhist(
            c, kt_hists.Ahists.at(i), kt_hists.Bhists.at(i),
            "../figures/ABhists/ABQ_"+addToFilename+"_"+
            std::to_string(i)+filetype,
            -1111, -1111, 0.0, 0.1
        );
        plot_Chist(
            c, kt_hists.Chists.at(i),
            "../figures/CQhists/CQ_"+addToFilename+"_"
            +std::to_string(i)+filetype,
            0, 2, 0.0, 0.4, fitpars
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
    
    TLegend *leg = new TLegend(0.1, 0.8, 0.2, 0.9); // (x1, y1, x2, y2) in NDC coordinates
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
    Double_t Qmin, Double_t Qmax,
    std::vector<Double_t> fitpars
) {
    c->cd();
    Chist->GetXaxis()->SetNdivisions(20,kFALSE);
    c->SetGrid();
    Chist->SetMinimum(histMin);
    Chist->SetMaximum(histMax);
    if (Qmin > 0 || Qmax > 0) {
        Chist->GetXaxis()->SetRangeUser(Qmin,Qmax);
    }
    
    for (int i = 0; i < 20; i++) {
        Chist->GetXaxis()->ChangeLabel(i, 45.0, -1.0, 32, -1, -1, "");
    }
    Chist->Draw();

    TLegend *leg = new TLegend(0.1, 0.75, 0.25, 0.9); // (x1, y1, x2, y2) in NDC coordinates
    leg->AddEntry(Chist, "measured C(Q)", "l");

    if (!fitpars.empty()) {
        plot_Levyfit(c, leg, Qmin, Qmax, fitpars);
        plot_fitRange(c, leg, fitpars.at(4), fitpars.at(5));
    }
    leg->Draw();

    c->Print(filename);
    c->Clear();
}

void plot_Levyfit(
    TCanvas* c, TLegend* leg, Double_t Qmin, Double_t Qmax,
    std::vector<Double_t> fitpars
) {
    std::cerr << "Plottin Lévy fit...\t";
    Int_t nPoints = 400;
    std::vector<Double_t> xvals(nPoints-1);
    std::vector<Double_t> yvals(nPoints-1);
    Double_t step = (Qmax-Qmin) / (nPoints-2);

    std::cerr << "Calculating Lévy function points with\n"
              << "Norm = " << fitpars.at(0) << '\t'
              << "alpha = " << fitpars.at(3) << '\t'
              << "R = " << fitpars.at(2) << '\t'
              << "lambda = " << fitpars.at(1) << '\n';
    for (int i = 1; i < nPoints; i++) {
        double x = Qmin + i*step;
        xvals.at(i-1) = x;
        yvals.at(i-1) = fitpars.at(0) * cccinstance->FullCorrFuncValueLambda(
            fitpars.at(3), fitpars.at(2), fitpars.at(1), x
        );
    }
    
    TGraph* gr = new TGraph(nPoints-1, xvals.data(), yvals.data());
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);
    gr->Draw("L");
    leg->AddEntry(gr, "Levy fit", "l");
}

void plot_fitRange(
    TCanvas* c, TLegend* leg, Double_t Qmin, Double_t Qmax
) {
    gPad->Update();
    Double_t ymin = gPad->GetUymin();
    Double_t ymax = gPad->GetUymax();

    TLine* line1 = new TLine(Qmin, ymin, Qmin, ymax);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw();

    TLine* line2 = new TLine(Qmax, ymin, Qmax, ymax);
    line2->SetLineColor(kRed);
    line2->SetLineWidth(2);
    line2->SetLineStyle(2);
    line2->Draw();

    leg->AddEntry(line2, "Fit range", "l");
}