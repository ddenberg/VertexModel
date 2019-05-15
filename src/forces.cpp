#include "forces.h"

int calcForces(LargeVec &R, int nVerts, std::vector<std::vector<int>> &cells, 
    std::vector<double> &A_0, double k_V, double k_P, double p_0, int yolk_id, LargeVec &F) {

    int nCells = cells.size();

    for (int m = 0; m < nCells; m++) {
        // Get center of cell m
        Vec3 c_m;
        int nVerts_m = cells[m].size();
        for (int i = 0; i < nVerts_m; i++) {
            // c_m += R_vec(R, cells[m][i], nVerts);
            c_m += Vec3(R[cells[m][i]], R[cells[m][i] + nVerts], R[cells[m][i] + nVerts * 2]);
        }
        c_m /= nVerts_m;

        
        //Create perimeter and areas for cell m
        double p_m = 0, a_m = 0;

        //Create vector of gradients of perimeter and area for cell m
        std::vector<Vec3> a_m_grad_r_v(nVerts_m);
        std::vector<Vec3> p_m_grad_r_v(nVerts_m);

        for (int v1 = 0; v1 < nVerts_m; v1++) {
            int p = v1;
            int q = (v1 + 1) % nVerts_m;

            // Vec3 r_p = R_vec(R, cells[m][p], nVerts);
            // Vec3 r_q = R_vec(R, cells[m][q], nVerts);
            Vec3 r_p = Vec3(R[cells[m][p]], R[cells[m][p] + nVerts], R[cells[m][p] + nVerts * 2]);
            Vec3 r_q = Vec3(R[cells[m][q]], R[cells[m][q] + nVerts], R[cells[m][q] + nVerts * 2]);

            Vec3 temp_vec1 = Cross(r_p - c_m, r_q - c_m);
            double temp_norm_vec1 = temp_vec1.Length();
            
            Vec3 temp_vec2 = r_p - r_q;
            double temp_norm_vec2 = temp_vec2.Length();
            
            Vec3 temp_vec3 = Cross(r_q - c_m, temp_vec1);
            Vec3 temp_vec4 = Cross(r_p - c_m, temp_vec1);
            
            Vec3 temp_grad00 = (1.0 / (double)nVerts_m) * (temp_vec4 - temp_vec3)  / (2.0 * temp_norm_vec1);

            //Update Perimeter and Area
            p_m += temp_norm_vec2;
            a_m += 0.5 * temp_norm_vec1;

            //%Calc Gradients
            for (int v2 = 0; v2 < nVerts_m; v2++) {
                bool d_jp = v2 == p, d_jq = v2 == q;
                if (!(d_jp == 0 && d_jq == 0)) {
                    a_m_grad_r_v[v2] += grad_A_Fun(d_jp, d_jq, temp_vec3, temp_vec4, nVerts_m, temp_norm_vec1);
                    p_m_grad_r_v[v2] += grad_L_Fun(d_jp, d_jq, temp_vec2, temp_norm_vec2);
                } else {
                    a_m_grad_r_v[v2] += temp_grad00;
                }
            }
        }

        for (int v1 = 0; v1 < nVerts_m; v1++) {
            Vec3 F_v1 = -2.0 * k_V * (a_m -  A_0[m]) * a_m_grad_r_v[v1];
            if (m != yolk_id) {
                F_v1 -= 2.0 * k_P * (p_m - p_0) * p_m_grad_r_v[v1];
            }
            
            F[cells[m][v1]] += F_v1.x();
            F[cells[m][v1] + nVerts] += F_v1.y();
            F[cells[m][v1] + 2 * nVerts] += F_v1.z();
        }
    }
    return 0;
}

int calcForces_parallel(LargeVec &R, int nVerts, std::vector<int> &vertPerCell_cumul, std::vector<std::vector<int>> &cells, 
    std::vector<double> &A_0, double k_V, double k_P, double p_0, int yolk_id, LargeVec &F_ind, LargeVec &F_vec, LargeVec &F) {

    int nCells = cells.size();
    // LargeVec F_vec(vertPerCell_cumul[nCells] * 3);
    // LargeVec F_ind(vertPerCell_cumul[nCells] * 3);


#pragma omp parallel for shared(R, F_vec, F_ind, cells, vertPerCell_cumul)
    for (int m = 0; m < nCells; m++) {
        // Get center of cell m
        Vec3 c_m;
        int nVerts_m = cells[m].size();
        for (int i = 0; i < nVerts_m; i++) {
            c_m += Vec3(R[cells[m][i]], R[cells[m][i] + nVerts], R[cells[m][i] + nVerts * 2]);
        }
        c_m /= nVerts_m;

        
        //Create perimeter and areas for cell m
        double p_m = 0, a_m = 0;

        //Create vector of gradients of perimeter and area for cell m
        std::vector<Vec3> a_m_grad_r_v(nVerts_m);
        std::vector<Vec3> p_m_grad_r_v(nVerts_m);

        for (int v1 = 0; v1 < nVerts_m; v1++) {
            int p = v1;
            int q = (v1 + 1) % nVerts_m;

            Vec3 r_p = Vec3(R[cells[m][p]], R[cells[m][p] + nVerts], R[cells[m][p] + nVerts * 2]);
            Vec3 r_q = Vec3(R[cells[m][q]], R[cells[m][q] + nVerts], R[cells[m][q] + nVerts * 2]);

            Vec3 temp_vec1 = Cross(r_p - c_m, r_q - c_m);
            double temp_norm_vec1 = temp_vec1.Length();
            
            Vec3 temp_vec2 = r_p - r_q;
            double temp_norm_vec2 = temp_vec2.Length();
            
            Vec3 temp_vec3 = Cross(r_q - c_m, temp_vec1);
            Vec3 temp_vec4 = Cross(r_p - c_m, temp_vec1);
            
            Vec3 temp_grad00 = (1.0 / (double)nVerts_m) * (temp_vec4 - temp_vec3)  / (2.0 * temp_norm_vec1);

            //Update Perimeter and Area
            p_m += temp_norm_vec2;
            a_m += 0.5 * temp_norm_vec1;

            //%Calc Gradients
            for (int v2 = 0; v2 < nVerts_m; v2++) {
                bool d_jp = v2 == p, d_jq = v2 == q;
                if (!(d_jp == 0 && d_jq == 0)) {
                    a_m_grad_r_v[v2] += grad_A_Fun(d_jp, d_jq, temp_vec3, temp_vec4, nVerts_m, temp_norm_vec1);
                    p_m_grad_r_v[v2] += grad_L_Fun(d_jp, d_jq, temp_vec2, temp_norm_vec2);
                } else {
                    a_m_grad_r_v[v2] += temp_grad00;
                }
            }
        }

        for (int v1 = 0; v1 < nVerts_m; v1++) {
            Vec3 F_v1 = -2.0 * k_V * (a_m -  A_0[m]) * a_m_grad_r_v[v1];
            if (m != yolk_id) {
                F_v1 -= 2.0 * k_P * (p_m - p_0) * p_m_grad_r_v[v1];
            }
            
            for (int d = 0; d < 3; d++) {
                F_ind[vertPerCell_cumul[m] * 3 + v1 * 3 + d] = cells[m][v1] + nVerts * d;
                F_vec[vertPerCell_cumul[m] * 3 + v1 * 3 + d] = F_v1.e[d];
            }
        }
    }

    for (int i = 0; i < 3 * vertPerCell_cumul[nCells]; i++) {
        F[F_ind[i]] += F_vec[i];
    }

    return 0;
}


inline Vec3 grad_A_Fun(int d_jp, int d_jq, Vec3 cross_vec1, Vec3 cross_vec2, int S_m, double norm_v) {
    return ((d_jp - 1.0 / (double)S_m) * cross_vec1 + (1.0 / (double)S_m - d_jq) * cross_vec2) / (2.0 * norm_v);
}

inline Vec3 grad_L_Fun(int d_jp, int d_jq, Vec3 v, double norm_v) {
    return (d_jp - d_jq) * (v / (norm_v + 1e-2));
}