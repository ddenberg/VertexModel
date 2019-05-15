#pragma once

#include "mesh.h"

int calcForces(LargeVec &R, int nVerts, std::vector<std::vector<int>> &cells, 
    std::vector<double> &A_0, double k_V, double k_P, double p_0, int yolk_id, LargeVec &F);

int calcForces_parallel(LargeVec &R, int nVerts, std::vector<int> &vertPerCell_cumul, std::vector<std::vector<int>> &cells, 
    std::vector<double> &A_0, double k_V, double k_P, double p_0, int yolk_id, LargeVec &F_ind, LargeVec &F_vec, LargeVec &F);
inline Vec3 grad_A_Fun(int d_jp, int d_jq, Vec3 cross_vec1, Vec3 cross_vec2, int S_m, double norm_v);
inline Vec3 grad_L_Fun(int d_jp, int d_jq, Vec3 v, double norm_v);