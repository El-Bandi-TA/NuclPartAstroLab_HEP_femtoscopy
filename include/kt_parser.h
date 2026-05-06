#ifndef KT_PARSER_H
#define KT_PARSER_H

#include "pool.h"
#include "event.h"
#include "pion.h"
#include "kt_pair_cuts.h"

#include <TROOT.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>

#include <string>
#include <vector>

class KTParserBase
{
private:
    void parse_config_file();
public:
    std::string KTconfig_file;
    Int_t num_KT_bins;
    Float_t KT_min;
    Float_t KT_max;
    std::vector<Float_t> KT_bin_edges;
    
    KTParserBase() = default;
    KTParserBase(
        std::string config_file
    ) : KTconfig_file(config_file), num_KT_bins(-1), KT_min(-1.0), KT_max(-1.0) {
        parse_config_file();
    }
    virtual ~KTParserBase() {}
    
    Int_t get_KT_index(Float_t KT);
};

template <typename T>
class KTParser : public KTParserBase {
private:
    KTPairCuts ZPhiCuts;
    void reserve_hists() {
        Ahists.reserve(num_KT_bins);
        Ahists.resize(num_KT_bins);
        Bhists.reserve(num_KT_bins);
        Bhists.resize(num_KT_bins);
        Chists.reserve(num_KT_bins);
        Chists.resize(num_KT_bins);
    }
public:
    Event current_event;
    Int_t Nemptypionvects;
    Pool pool;
    std::vector<T*> Ahists;
    std::vector<T*> Bhists;
    std::vector<T*> Chists;
    TString name;
    TString axis_titles;

    KTParser() : KTParserBase() {
        reserve_hists();
    }
    KTParser(
        std::string config_file
    ) : KTParserBase(config_file), Nemptypionvects(0), pool() {
        reserve_hists();
    }
    KTParser(
        std::string KT_config_file, std::string pair_config_file
    ) : KTParserBase(KT_config_file), Nemptypionvects(0), pool(),
        ZPhiCuts(pair_config_file, KT_bin_edges) {
            reserve_hists();
    }
    KTParser(
        std::string config_file, TString hname, TString haxis_titles,
        Int_t pminc, Int_t pmaxc, Int_t pnumc,
        Float_t pminz, Float_t pmaxz, Int_t pnumz
    ) : KTParserBase(config_file), Nemptypionvects(0),
        pool(pminc, pmaxc, pnumc, pminz, pmaxz, pnumz),
        name(hname), axis_titles(haxis_titles) {
            reserve_hists();
    }
    KTParser(
        std::string KT_config_file, std::string pair_config_file,
        TString hname, TString haxis_titles,
        Int_t pminc, Int_t pmaxc, Int_t pnumc,
        Float_t pminz, Float_t pmaxz, Int_t pnumz
    ) : KTParserBase(KT_config_file), Nemptypionvects(0),
        pool(pminc, pmaxc, pnumc, pminz, pmaxz, pnumz),
        name(hname), axis_titles(haxis_titles),
        ZPhiCuts(pair_config_file, KT_bin_edges) {
            reserve_hists();
    }
    KTParser(
        std::string config_file, TString hist_file,
        TString hname, TString haxis_titles
    ) : KTParserBase(config_file), name(hname), axis_titles(haxis_titles) {
        load_hists(hist_file, hname);
    }

    ~KTParser() {
        for (auto h: Ahists) { if (h) { delete h; } }
        Ahists.clear();
        for (auto h: Bhists) { if (h) { delete h; } }
        Bhists.clear();
        for (auto h: Chists) { if (h) { delete h; } }
        Chists.clear();
    };

    Int_t get_pool_index();
    void add_pion(const Pion& pi);
    void clear_current_event();
    void print_num_pions();
    void print_kt_bin_edges();

    void initialize_ABhists(
        Int_t num_bins_x, Double_t x_min, Double_t x_max,
        Int_t num_bins_y=-1, Double_t y_min=-1, Double_t y_max=-1
    );
    void create_Chist(Int_t kt_idx);
    void create_Chists();
    void fill_Ahists();
    void fill_Bhists();
    void scale_ABhists();
    void scale_hists(Double_t scaleFactor, char which);
    void scale_hists(const std::vector<Double_t>& scaleFactor, char which);
    void scale_Chists();
    void write_hists(TFile* out_file);
    void load_hists(TString hist_file, TString hist_name);
    void merge_hists(Int_t idx1, Int_t idx2);
    void merge_hist_bins(Int_t nMergeX, Int_t nMergeY=1);
};

#endif //KT_PARSER_H