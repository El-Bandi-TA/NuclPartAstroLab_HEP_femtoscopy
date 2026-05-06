#ifndef PION_H
#define PION_H

#include <TROOT.h>
#include <TMath.h>

class Pion
{
public:
    static constexpr Float_t PI_MASS = 0.13957039;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t det_z;
    Float_t det_phi;
    Float_t p_abs;
    Float_t KT;
    Float_t E;
    Float_t eta;
    Int_t charge;
    
    Pion() = default;
    Pion(
        Float_t ppx, Float_t ppy, Float_t ppz, Int_t pch
    );
    Pion(
        Float_t ppx, Float_t ppy, Float_t ppz,
        Float_t pdet_z, Float_t pdet_phi, Int_t pch
    );
    ~Pion() = default;

    void calc_pion_energy() {
        E = TMath::Sqrt(
            px*px + py*py + pz*pz + PI_MASS*PI_MASS
        );
    }

    void calc_eta() {
        eta = 0.5 * TMath::Log((E+pz)/(E-pz));
    }

    void calc_p_abs() {
        p_abs = TMath::Sqrt(px*px + py*py + pz*pz);
    }

    void calc_KT() {
        KT = TMath::Sqrt(px*px+py*py);
    }
};

#endif //PION_H