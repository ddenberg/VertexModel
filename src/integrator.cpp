#include "integrator.h"

//////////////////////////////////////////////////////////////////////////////
///////////////////////////////  BaseSolver   ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void BaseSolver::initVars(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc) {
    step = 0;
    t = 0.;
    dt = _dt;

    write_step = 0;
    write_time = 0;
    write_dt = _write_dt;

    nDof = cell_mesh.numVerts() * 3;
    F.Init(nDof);
    R.Init(nDof);
    F_ind.Init(cell_mesh.vertPerCell_cumul[cell_mesh.numCells()] * 3);
    F_vec.Init(cell_mesh.vertPerCell_cumul[cell_mesh.numCells()] * 3);

    use_parallel = _use_parallel;
    apply_bc = _apply_bc;
}

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////  ExplicitSolver   //////////////////////////////
//////////////////////////////////////////////////////////////////////////////
ExplicitSolver::ExplicitSolver(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc) {
    initVars(cell_mesh, _dt, _write_dt, _use_parallel, _apply_bc);
}

int ExplicitSolver::run(Mesh &cell_mesh, double t_end) {
    int nVerts = cell_mesh.numVerts();
    for (int i = 0; i < nVerts; i++) {
        R[i] = cell_mesh.r[i].x();
        R[i + nVerts] = cell_mesh.r[i].y();
        R[i + 2 * nVerts] = cell_mesh.r[i].z();
    }
    F.setZero();
    std::string filename = "output/out" + std::to_string(write_step) + ".vtk";
    cell_mesh.writeMesh(filename);

    while ((t_end - t) > 1e-6) {
        t += dt;
        step++;

        if (use_parallel)
            calcForces_parallel(R, nVerts, cell_mesh.vertPerCell_cumul, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F_ind, F_vec, F);
        else
            calcForces(R, nVerts, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F);
        

        //Update positions
        R += (dt * F);
        F.setZero();
        if (apply_bc)
            applyBC(R, nVerts);

        if ((step % 1000) == 0)
            std::cout << "t = " << t << std::endl;

        if (t >= write_time) {
            for (int i = 0; i < nVerts; i++) {
                cell_mesh.r[i] = Vec3(R[i], R[i+nVerts], R[i+2*nVerts]);
            }
            write_step++;
            write_time += write_dt;
            filename = "output/out" + std::to_string(write_step) + ".vtk";
            cell_mesh.writeMesh(filename);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////  Runge Kutta 4 Solver   ///////////////////////////
//////////////////////////////////////////////////////////////////////////////
RungeKutta4Solver::RungeKutta4Solver(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc) {
    initVars(cell_mesh, _dt, _write_dt, _use_parallel, _apply_bc);
    initRK4Vars();
}

void RungeKutta4Solver::initRK4Vars() {
    k1.Init(nDof);
    k2.Init(nDof);
    k3.Init(nDof);
    R_temp.Init(nDof);
}

int RungeKutta4Solver::run(Mesh &cell_mesh, double t_end) {
    int nVerts = cell_mesh.numVerts();
    for (int i = 0; i < nVerts; i++) {
        R[i] = cell_mesh.r[i].x();
        R[i + nVerts] = cell_mesh.r[i].y();
        R[i + 2 * nVerts] = cell_mesh.r[i].z();
    }
    F.setZero();
    std::string filename = "output/out" + std::to_string(write_step) + ".vtk";
    cell_mesh.writeMesh(filename);

    while ((t_end - t) > 1e-6) {
        t += dt;
        step++;

        //Update k1
        if (use_parallel)
            calcForces_parallel(R, nVerts, cell_mesh.vertPerCell_cumul, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F_ind, F_vec, F);
        else
            calcForces(R, nVerts, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F);
        k1 = dt * F;
        R_temp = 0.5 * k1 + R;
        F.setZero();

        //Update k2
        if (use_parallel)
            calcForces_parallel(R_temp, nVerts, cell_mesh.vertPerCell_cumul, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F_ind, F_vec, F);
        else
            calcForces(R_temp, nVerts, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F);
        k2 = dt * F;
        R_temp = 0.5 * k2 + R;
        F.setZero();

        //Update k3
        if (use_parallel)
            calcForces_parallel(R_temp, nVerts, cell_mesh.vertPerCell_cumul, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F_ind, F_vec, F);
        else
            calcForces(R_temp, nVerts, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F);
        k3 = dt * F;
        R_temp = k3 + R;
        F.setZero();

        //Update k4/R
        if (use_parallel)
            calcForces_parallel(R_temp, nVerts, cell_mesh.vertPerCell_cumul, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F_ind, F_vec, F);
        else
            calcForces(R_temp, nVerts, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F);
        R += (1. / 6.) * (k1 + 2.0 * k2 + 2.0  * k3 + dt * F); // R_{n_1} = R_{n} + 1/6*(k1 + 2*k2 + 2*k3 + k4)
        F.setZero();

        if (apply_bc)
            applyBC(R, nVerts);

        if ((step % 1000) == 0)
            std::cout << "t = " << t << std::endl;

        if (t >= write_time) {
            for (int i = 0; i < nVerts; i++) {
                cell_mesh.r[i] = Vec3(R[i], R[i+nVerts], R[i+2*nVerts]);
            }
            write_step++;
            write_time += write_dt;
            filename = "output/out" + std::to_string(write_step) + ".vtk";
            cell_mesh.writeMesh(filename);
        }
    }
}

//////////////////////////////////////////////////////////////////
Implicit_LBFGS_Solver::Implicit_LBFGS_Solver(Mesh &cell_mesh, double _dt, double _write_dt, int _use_parallel, int _apply_bc) {
    initVars(cell_mesh, _dt, _write_dt, _use_parallel, _apply_bc);
}

int Implicit_LBFGS_Solver::run(Mesh &cell_mesh, double t_end) {
    int nVerts = cell_mesh.numVerts();
    for (int i = 0; i < nVerts; i++) {
        R[i] = cell_mesh.r[i].x();
        R[i + nVerts] = cell_mesh.r[i].y();
        R[i + 2 * nVerts] = cell_mesh.r[i].z();
    }
    F.setZero();
    std::string filename = "output/out" + std::to_string(write_step) + ".vtk";
    cell_mesh.writeMesh(filename);

    double hess_guess = 1.0;
    while ((t_end - t) > 1e-6) {
        t += dt;
        step++;

        if (L_BFGS(hess_guess, cell_mesh, 0) != 0) {
            std::cout << "Failed to Converge, t = " << t << "\n";
            return -1;
        }
        if (apply_bc)
            applyBC(R, nVerts);

        // std::cout << "\n";

        if ((step % 1000) == 0)
            std::cout << "t = " << t << std::endl;

        if (t >= write_time) {
            for (int i = 0; i < nVerts; i++) {
                cell_mesh.r[i] = Vec3(R[i], R[i+nVerts], R[i+2*nVerts]);
            }
            write_step++;
            write_time += write_dt;
            filename = "output/out" + std::to_string(write_step) + ".vtk";
            cell_mesh.writeMesh(filename);
        }
    }
}

void Implicit_LBFGS_Solver::gradient(Mesh &cell_mesh, LargeVec &x, LargeVec &x_prev, LargeVec &grad) {
    if (use_parallel)
        calcForces_parallel(x, cell_mesh.numVerts(), cell_mesh.vertPerCell_cumul, cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F_ind, F_vec, F);
    else
        calcForces(x, cell_mesh.numVerts(), cell_mesh.cells, cell_mesh.A_0, k_V, k_P, p_0, cell_mesh.yolk_id, F);
    

    grad = x - x_prev - dt * F;
    F.setZero();
}

int Implicit_LBFGS_Solver::L_BFGS(double &hess_guess, Mesh &cell_mesh, int print) {

    LargeVec alpha = LargeVec(m);
    alpha.setZero();
	LargeVec rho = LargeVec(m);
    rho.setZero();
    LargeVec grad(nDof), q(nDof), grad_old(nDof), x_old(nDof), x0(nDof), x_prev(nDof);

    //Move R into x;
    x_prev = R;
    x0 = x_prev;

    gradient(cell_mesh, x0, x_prev, grad);

    VecContainer s(m), y(m);
    for (int i = 0; i < m; i++) {
        s[i].Init(nDof);
        y[i].Init(nDof);
    }

    double gamma_k = hess_guess, gradNorm = 0.0;
	double alpha_init = std::min(1.0, 1.0 / grad.norm());


	int globalIter = 0;

    int k_max = 1000;
    int k = 0;
    bool useGrad = false;
    while (k < k_max) {
        x_old = x0;
        grad_old = grad;
        q = grad;
        globalIter++;

        int iter = std::min(m, k);
		for (int i = iter - 1; i >= 0; --i) {
            rho[i] = 1.0 / s[i].dot(y[i]);
            alpha[i] = rho[i] * s[i].dot(q);
            q -= (alpha[i] * y[i]);
		}

        q = q * gamma_k;

		for (int i = 0; i < iter; i++) {
			double beta = rho[i] * q.dot(y[i]);
			q = q + (alpha[i] - beta) * s[i];
		}

		double dir = q.dot(grad);
        useGrad = false;
		if (dir < 1e-2 * grad.norm() * q.norm()){
            useGrad = true;
			q = grad;
			k_max -= k;
			k = 0;
			alpha_init = std::min(1.0, 1.0 / grad.norm());
		}

        double rate = alpha_init;
		x0 = x0 - rate * q;

        double diffNorm = (x_old - x0).norm() / rate;
		if (diffNorm < tolDiffX){
            if (print > 0) {
                if (useGrad) {
                    printf("t = %f  GD_Iteration: %d    ResNorm = DNC    DiffNorm = %e    rate = %e\n", t, globalIter, diffNorm, rate);
                } else {
                    printf("t = %f    Iteration: %d    ResNorm = DNC     DiffNorm = %e    rate = %e\n", t, globalIter, diffNorm, rate);
                }
            }
			break;
		}

		gradient(cell_mesh, x0, x_prev, grad);

		gradNorm = grad.norm();

        if (print > 0) {
            if (useGrad) {
                printf("t = %f  GD_Iteration: %d    ResNorm = %e     DiffNorm = %e    rate = %e\n", t, globalIter, gradNorm, diffNorm, rate);
            } else {
                printf("t = %f    Iteration: %d    ResNorm = %e     DiffNorm = %e    rate = %e\n", t, globalIter, gradNorm, diffNorm, rate);
            }
        }

		if (gradNorm < tolDiffGrad){
			hess_guess = gamma_k;
            break;
		}

		if (k < m) {
			s[k] = x0 - x_old;
			y[k] = grad - grad_old;

            gamma_k = s[k].dot(y[k]) / y[k].normSquared();
		} else {
            s.rotate(x0 - x_old);
            y.rotate(grad - grad_old);

            gamma_k = s[m-1].dot(y[m-1]) / y[m-1].normSquared();
		}

		alpha_init = 1.0;
        k++;
    }

    if (k == k_max) {
        return -1;
    }


    R = x0;
    return 0;
}