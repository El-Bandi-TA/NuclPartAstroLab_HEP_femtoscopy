#ifndef POOL_H
#define POOL_H

#include "event.h"

#include <TROOT.h>
#include <vector>
#include <list>

class Pool
{
private:
    /* data */
public:
    Int_t min_centrality;
    Int_t max_centrality;
    Int_t num_centrality_class;
    Float_t min_zvertex;
    Float_t max_zvertex;
    Int_t num_zvertex_class;
    std::vector<std::list<Event>> pool {};
    
    Int_t get_index(Int_t centrality, Float_t zvertex);

    Pool() = default;
    Pool(
        Int_t minc, Int_t maxc, Int_t numc,
        Float_t minz, Float_t maxz, Int_t numz
    ) : min_centrality(minc), max_centrality(maxc), num_centrality_class(numc),
        min_zvertex(minz), max_zvertex(maxz), num_zvertex_class(numz) {
            pool.reserve(numc*numz);
            pool.resize(numc*numz);
        }
    ~Pool() = default;
};

inline Int_t Pool::get_index(Int_t centrality, Float_t zvertex) {
    if(
        centrality < min_centrality || centrality >= max_centrality
        || zvertex < min_zvertex || zvertex >= max_zvertex
    ) {return -1;}
    Int_t idx = (
        num_centrality_class * std::min(
            (int)((zvertex-min_zvertex) * num_zvertex_class
            / (max_zvertex-min_zvertex)), 9
        )
        + (centrality-min_centrality) * num_centrality_class
        / (max_centrality-min_centrality)
    );

    return idx;
}

#endif //POOL_H