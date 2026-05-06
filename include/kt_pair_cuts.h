#ifndef KT_PAIR_CUTS_H
#define KT_PAIR_CUTS_H

#include <TROOT.h>

#include <string>
#include <vector>

struct CutRect {
    Float_t dzmin = 0.0f;
    Float_t dzmax = 0.0f;
    Float_t dphimin = 0.0f;
    Float_t dphimax = 0.0f;

    CutRect(
        Float_t z1, Float_t z2, Float_t p1, Float_t p2
    ) : dzmin(z1), dzmax(z2), dphimin(p1), dphimax(p2) {}

    Bool_t contains(Float_t dz, Float_t dphi) const {
        return !(dz < dzmin || dz > dzmax || dphi < dphimin || dphi > dphimax);
    }
};

class KTPairCuts
{
private:
    void parse_config_file(std::string pairConfigFile);
    void initialize_bin_map(std::vector<Float_t> KT_bin_edges);
    std::vector<Int_t> bin_map;
    std::vector<std::vector<CutRect>> cuts;
    std::vector<Float_t> KT_cut_edges;
public:
    KTPairCuts() = default;
    KTPairCuts(std::string pairConfigFile, std::vector<Float_t> KT_bin_edges) {
        parse_config_file(pairConfigFile);
        initialize_bin_map(KT_bin_edges);
        print_cuts(KT_bin_edges);
    }
    ~KTPairCuts() = default;

    Bool_t isRejected(Int_t ktbin, Float_t dz, Float_t dphi) {
        Int_t cutbin = bin_map.at(ktbin);
        const auto& active_cuts = cuts.at(cutbin);

        for(const auto& cutrect: active_cuts) {
            if (cutrect.contains(dz, dphi)) return kTRUE;
        }
        return kFALSE;
    }

    void print_cuts(std::vector<Float_t> KT_bin_edges);
};

#endif //KT_PAIR_CUTS_H
