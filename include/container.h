#pragma once

#include <cmath>
#include "largevec.h"

class VecContainer {
private:
    LargeVec *vectors;
    int *list;
    int n;
public:
    VecContainer(const int _n);
    ~VecContainer() { delete[] vectors; delete[] list; }

    inline LargeVec operator[](int i) const { return vectors[list[i]]; }
    inline LargeVec& operator[](int i) { return vectors[list[i]]; }

    void rotate(const LargeVec &v);
};