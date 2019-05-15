#include "container.h"

VecContainer::VecContainer(const int _n) {
    n = _n;
    vectors = new LargeVec[_n]; 
    list = new int[_n];
    for (int i = 0; i < n; i++) {
        list[i] = i;
    }
}

void VecContainer::rotate(const LargeVec &v) {
    for (int i = 0; i < n; i++) {
        list[i] = list[(i+1) % n];
    }
    vectors[n-1] = v;
}