#include "pion.h"
#include <TROOT.h>
#include <TMath.h>

Pion::Pion(
    Float_t ppx, Float_t ppy, Float_t ppz, Int_t pch
) : px(ppx), py(ppy), pz(ppz), det_z(-1e10), det_phi(-1e10), charge(pch) {
    calc_pion_energy();
    calc_eta();
    calc_p_abs();
    calc_KT();
}

Pion::Pion(
    Float_t ppx, Float_t ppy, Float_t ppz,
    Float_t pdet_z, Float_t pdet_phi, Int_t pch
) : px(ppx), py(ppy), pz(ppz), det_z(pdet_z), det_phi(pdet_phi), charge(pch) {
    calc_pion_energy();
    calc_eta();
    calc_p_abs();
    calc_KT();
}