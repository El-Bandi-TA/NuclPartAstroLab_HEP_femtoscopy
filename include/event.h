#ifndef EVENT_H
#define EVENT_H

#include "pion.h"

#include <TROOT.h>
#include <vector>

class Event
{
public:
    std::vector<Pion> pions;
    Int_t centrality;
    Float_t Zvertex;

    void add_pion(const Pion& pi);
    void clear() {
        centrality=-1;
        Zvertex=-1e10;
        pions.clear();
    }
    
    Event() : centrality(-1), Zvertex(-1e10) {
        pions.reserve(30);
    }
    Event(
        Int_t pcentrality, Float_t pZvertex
    ) : centrality(pcentrality), Zvertex(pZvertex) {}
    ~Event() = default;
};

#endif //EVENT_H