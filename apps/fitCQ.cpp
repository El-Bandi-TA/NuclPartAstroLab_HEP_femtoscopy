#include "CoulCorrCalc.h"
#include "kt_parser.h"

#include <TH1D.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

CoulCorrCalc *cccinstance;

const int NPARS = 4;

const char *statuses[6] = {
    "converged",
    "cov. made pos.def.",
    "Hesse invalid",
    "Edm above max",
    "call lim. reached",
    "other failure"
};
const char *covstatuses[5] = {
    "not available",
    "not pos.def.",
    "approximate",
    "forced pos.def.",
    "accurate"
};

struct FitResults {
    bool success;
    std::vector<double> params;
    std::vector<double> errors;
    std::vector<double> minosPlus;
    std::vector<double> minosMinus;
    int ndf;
    double chi2;
    double probability;
    int fitstatus;
    int covstatus;

    FitResults() : success(0), chi2(-1), probability(-1), fitstatus(-1),
                   covstatus(-1), ndf(-1), params(NPARS, 0), errors(NPARS, 0),
                   minosPlus(NPARS, 0), minosMinus(NPARS, 0) {}

    void print_results();
    void print_status();
    void write_results(
        std::ostream& outfile,
        double KTmin, double KTmax,
        double Qmin, double Qmax
    );
};

std::string get_current_time();

double FitFunction(const double *x, const double *par);
double MyChi2(
    const double *par, const std::vector<double>& Qvals,
    const std::vector<double>& expvals, const std::vector<double>& experrs
);
FitResults fit_CoulCorrLevy(
    TH1D* hist, double Qmin, double Qmax, bool calcMinos
);

void write_header(std::ostream& out_file);

int main(int argc, char** argv) {

    const char* env_p = std::getenv("HEPQC_ROOT");
    if (!env_p) {
        std::cerr << "Error: HEPQC_ROOT not set!\n";
    }
    std::string projRoot = env_p;

    // --- Parse CLI arguments ---
    std::string inHistFile = projRoot+"/data/Q_full/Q_plus.root";
    std::string outFilename = (
        projRoot+"/data/fitparams_"+get_current_time()+".txt"
    );
    int mergeQBins = 1;
    int startFrom = 0;
    int fitN = 0;
    bool calcMinos = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " "
                      << "[-i|--input "
                      << "inHistFile=HEPQC_ROOT/data/Q_full/Q_plus.root] "
                      << "[-o|--output outFilename="
                      << "HEPQC_ROOT/data/fitparams_CURRENTTIME.txt] "
                      << "[-b|--merge-bins mergeQBins=1]"
                      << "[-n|--fitN fitN=0] "
                      << "[-s|--start-from startFrom=0] "
                      << "[-m|--calc-minos calcMinos=false]"
                      << std::endl;
            return -1;
        } else if ((arg == "-i" || arg == "--input") && i+1 < argc) {
            inHistFile = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i+1 < argc) {
            outFilename = argv[++i];
        } else if ((arg == "-b" || arg == "--merge-bins")) {
            mergeQBins = std::atoi(argv[++i]);
        } else if ((arg == "-s" || arg == "--start-from") && i+1 < argc) {
            startFrom = std::atoi(argv[++i]);
            if (startFrom < 0) startFrom = 0;
        } else if ((arg == "-n" || arg == "--fitN")) {
            fitN = std::atoi(argv[++i]);
            if (fitN < 0) fitN = 0;
        } else if ((arg == "-m" || arg == "--calc-minos")) {
            calcMinos = true;
        }
    }

    std::vector<double> Qmin {};
    std::vector<double> Qmax {};

    // --- Parsing finished ---

    // Create output file and check wether it is valid
    std::ofstream outParamFile(outFilename);

    if (!outParamFile.is_open()) {
        std::cerr << "FATAL: Could not open output file " << outFilename
                  << "for writing\n";
        return 1;
    }
    write_header(outParamFile);

    // Initializing KTParser
    KTParser<TH1D> KThists(
        projRoot+"/config/KTbins.txt",
        inHistFile,
        "Q", "(Q); Q [GeV]; Counts"
    );
    // KThists.merge_hists(0,3);
    if (mergeQBins>1) {
        KThists.merge_hist_bins(mergeQBins);
        KThists.scale_ABhists();
        KThists.create_Chists();
    };

    // Initialize Coulomb Fourier calculation instance
    cccinstance = new CoulCorrCalc();

    // Calculate how many to fit
    int fit = KThists.num_KT_bins-startFrom;
    if (fitN>0 && fitN<fit) fit = fitN;

    // Parse Q fit range config file
    std::string QfitrangeConfigFilename = projRoot+"/config/Qfitrange.txt";
    std::ifstream QfitrangeConfig(QfitrangeConfigFilename);
    if (!QfitrangeConfig.is_open()) {
        std::cerr << "FATAL: Could not open Q fit range config file "
                  << QfitrangeConfigFilename << ". Aborting...";
        return 1;
    }
    double qmin, qmax;
    int bin;
    std::string line {};
    while (std::getline(QfitrangeConfig,line)) {
        if (line[0]=='#') continue;
        std::stringstream ss {line};
        ss >> bin >> qmin >> qmax;
        Qmin.push_back(qmin);
        Qmax.push_back(qmax);
    }

    // Doing fits
    for(int i = startFrom; i < fit+startFrom; i++) {
        std::cerr << "Fitting hist " << i << " between " << Qmin.at(i) << " GeV"
                  << " and " << Qmax.at(i) << " GeV...\n";
        FitResults res = fit_CoulCorrLevy(
            KThists.Chists.at(i), Qmin.at(i), Qmax.at(i), calcMinos
        );
        res.print_results();
        std::cerr << "-------------------------------------------" << std::endl;
        res.write_results(
            outParamFile,
            KThists.KT_bin_edges.at(i),
            KThists.KT_bin_edges.at(i+1),
            Qmin.at(i), Qmax.at(i)
        );
    }

    return 0;
}

double FitFunction(const double *x, const double *par) {
    double N      = par[0];
    double lambda = par[1];
    double R      = par[2];
    double alpha  = par[3];
    double Q      = x[0];
    double corrfunc = cccinstance->FullCorrFuncValueLambda(alpha, R, lambda, Q);
    return N*corrfunc;
}

double MyChi2(
    const double *par, const std::vector<double>& Qvals,
    const std::vector<double>& expvals, const std::vector<double>& experrs
) {
    double chi2 = 0;
    for(std::size_t i = 0; i < Qvals.size(); i++)
    {
        double Q = Qvals.at(i);
        double exp = expvals.at(i);
        double err = experrs.at(i);
        double theor = FitFunction(&Q,par);
        if(err==0) continue;
        double chi = (exp-theor)/err;
        chi2 += chi*chi;
    }
    return chi2;
}

FitResults fit_CoulCorrLevy(
    TH1D* hist, double Qmin, double Qmax, bool calcMinos
) {
    std::vector<double> Q_values;
    std::vector<double> exp_values;
    std::vector<double> exp_errors;

    // 1. Find the bin indices for your range
    int binStart = hist->GetXaxis()->FindBin(Qmin);
    int binEnd   = hist->GetXaxis()->FindBin(Qmax);
    int nBins    = binEnd - binStart + 1;

    // 2. Pre-allocate memory (good for performance)
    Q_values.reserve(nBins);
    exp_values.reserve(nBins);
    exp_errors.reserve(nBins);

    // 3. Fill the vectors
    for (int i = binStart; i <= binEnd; ++i) {
        double err = hist->GetBinError(i);

        if (err > 1e-14) {
            Q_values.push_back(hist->GetXaxis()->GetBinCenter(i));
            exp_values.push_back(hist->GetBinContent(i));
            exp_errors.push_back(err);
        }
    }
    
    std::cerr << "Initializing minimizer...\n";
    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kCombined );

    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(0.001);

    auto chi2Lambda = [Q_values, exp_values, exp_errors](const double *par) {
        return MyChi2(par, Q_values, exp_values, exp_errors);
    };

    ROOT::Math::Functor f(chi2Lambda,NPARS);

    min.SetFunction(f);

    // Set the free variables to be minimized
    // start from completely out of context values to test fit accuracy
    min.SetVariable(0,"N",     0.95,0.01);
    min.SetVariable(1,"lambda",1.00,0.01);
    min.SetVariable(2,"R",     9.20,0.01);
    min.SetVariable(3,"alpha", 1.20,0.01);
    // min.SetLowerLimitedVariable(0,"N",0.95,0.01,0);
    // min.SetVariable(1,"lambda",1.00,0.01);
    // min.SetLowerLimitedVariable(2,"R", 7.50,0.01,0);
    // min.SetLimitedVariable(3,"alpha", 1.30,0.01,1,2);

    // min.SetFixedVariable(0,"x0",91);
    
    FitResults res;

    std::cerr << "Run minimizer...\n";
    res.success = min.Minimize();
    std::cerr << "Summarizing results...\n";
    // min.PrintResults();
    min.ProvidesError();
    min.Hesse();

    res.ndf = Q_values.size()-NPARS;
    res.chi2 = chi2Lambda(min.X());
    res.probability = TMath::Prob(res.chi2, res.ndf);
    for(int i = 0; i < NPARS; i++) {
        res.params.at(i) = min.X()[i];
        res.errors.at(i) = min.Errors()[i];
        
        if (calcMinos) {
            std::cerr << "Calculating minos error for parameter "
                      << i << "...\n";
            double errdn;
            double errup;
            if (min.GetMinosError(i,errdn,errup)) {
                res.minosMinus.at(i) = errdn;
                res.minosPlus.at(i) = errup;
            }
        }
    }

    res.fitstatus = min.Status();
    res.covstatus = min.CovMatrixStatus();
    if(res.fitstatus<0 || res.fitstatus>5) res.fitstatus=5;
    if(res.covstatus<-1 || res.covstatus>3) res.covstatus=-1;

    return res;
}

void write_header(std::ostream& out_file) {
    out_file << "#Params: N, lambda, R, alpha\n";
    out_file << "#KTmin\tKTmax\tQmin\tQmax\tSuccess\t";
    for(int i = 0; i < NPARS; i++) {
        out_file << "Par" << i << "\t";
        out_file << "Err" << i << "\t";
        out_file << "MinosUP" << i << "\t";
        out_file << "MinosDN" << i << "\t";
    }
    out_file << "NDF\tchi2\tprob\tfitstatus\tcovstatus";
}

void FitResults::print_status() {
    if(fitstatus<0 || fitstatus>5) fitstatus=5;
    if(covstatus<-1 || covstatus>3) covstatus=-1;
    if(fitstatus == 0 && covstatus == 3) {
        std::cout << "Fit converged, full accurate cov. matrix";
    } else {
        std::cout << "fit status: " << statuses[fitstatus] << std::endl;
        std::cout << "cov. matrix " << covstatuses[covstatus+1] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "(fitstatus=" << fitstatus << ",covstatus=" << covstatus
              << ")" << std::endl;
}

void FitResults::print_results() {
    std::cout << "Probability: " << chi2 << "/" << ndf << "->"
              << probability << endl;
    std::cout << "Parameters:" << std::endl;
    for(int ipar=0;ipar<NPARS;ipar++) {
        std::cout << "par" << ipar << "=" << params[ipar] << "+-" 
                  << errors[ipar] << std::endl;
    }

    print_status();

    std::cout << "Minos errors: " << std::endl;
    for(unsigned int ipar=0; ipar<NPARS; ipar++) {
        std::cout << "err" << ipar << ": +" << minosPlus.at(ipar) << " -"
                  << minosMinus.at(ipar) << std::endl;
    }
}

void FitResults::write_results(
    std::ostream& outfile, double KTmin, double KTmax, double Qmin, double Qmax
) {
    outfile << '\n' << KTmin << '\t' << KTmax << '\t'
            << Qmin << '\t' << Qmax << '\t' << success << '\t';
    for(int i = 0; i < NPARS; i++) {
        outfile << params.at(i) << "\t";
        outfile << errors.at(i) << "\t";
        outfile << minosPlus.at(i) << "\t";
        outfile << minosMinus.at(i) << "\t";
    }
    outfile << ndf << '\t' << chi2 << '\t' << probability << '\t';
    outfile << statuses[fitstatus] << '\t' << covstatuses[covstatus+1];
}
