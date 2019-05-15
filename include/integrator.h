#pragma once

#include "mesh.h"
#include "forces.h"
#include "container.h"
#include "boundary.h"

class BaseSolver {
public:
    virtual int run(Mesh &cell_mesh, double t_end) = 0;

protected:
    void initVars(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc);

    LargeVec F, R;
    LargeVec F_ind, F_vec;
    int nDof;

    double k_V = 1e2;
    double k_P = 1;
    double p_0 = 0;
    double t;
    double dt;
    int step;

    int write_step;
    double write_dt;
    double write_time;
    
    int apply_bc;
    int use_parallel;
};

class ExplicitSolver : public BaseSolver {
public:
    ExplicitSolver(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc);

    int run(Mesh &cell_mesh, double t_end);
};

class RungeKutta4Solver : public BaseSolver {
public:
    RungeKutta4Solver(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc);

    int run(Mesh &cell_mesh, double t_end);

private:
    void initRK4Vars();

    LargeVec k1, k2, k3, R_temp;
};

class Implicit_LBFGS_Solver : public BaseSolver {
public:
    Implicit_LBFGS_Solver(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc);

    int run(Mesh &cell_mesh, double t_end);

private:
    int L_BFGS(double &hess_guess, Mesh &cell_mesh, int print);
    void gradient(Mesh &cell_mesh, LargeVec &x, LargeVec &x_prev, LargeVec &grad);

    int m = 10;
    double tolDiffX = 1e-6;
    double tolDiffGrad = 1e-6;
};