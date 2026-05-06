#include "event.h"
#include "pion.h"

void Event::add_pion(const Pion& pi) {
    pions.push_back(pi);
}