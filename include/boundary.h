#pragma once

#include "mesh.h"

inline void applyBC(LargeVec &R, int nVerts) {
    double rad = 11.21;
    double rad_sq = rad * rad;

    Vec3 newpos;
    for (int i = 0; i < nVerts; i++) {
        if (R[i]*R[i] + R[i+nVerts]*R[i+nVerts] - rad_sq >= 0) {
            newpos = Vec3(R[i], R[i+nVerts], R[i+nVerts*2]);
            newpos.Make_Unit();
            newpos *= rad;
            
            R[i] = newpos.x();
            R[i+nVerts] = newpos.y();
            R[i+nVerts*2] = newpos.z();
        }
    }
}