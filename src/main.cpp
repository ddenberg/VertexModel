#include "mesh.h"
#include "integrator.h"

#include <sys/time.h>
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp () {
    struct timeval now;
    gettimeofday (&now, NULL);
    return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main() {
    double dt = 1e-3;
    double write_dt = 1e-1;
    int use_parallel = 1; //Set to 1 and re-compile to use parallel
    int apply_bc = 1; //Only set to 1 for the small cavity example


    Mesh cell_mesh = Mesh();

    //////////////////////////////////////////////////
    //Un-comment the following code block (and re-compile) to run the small voronoi mesh example
    //////////////////////////////////////////////////
    // cell_mesh.importMesh("init_voronoi.vtk");
    // std::vector<double> cell_areas = cell_mesh.cellAreas();
    // cell_mesh.setA_0(cell_areas);

    // //Un-comment to use the explicit solver
    // ExplicitSolver solver = ExplicitSolver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    // //Un-comment to use the RK4 solver
    // // RungeKutta4Solver solver = RungeKutta4Solver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    // //Un-comment to use the L-BFGS solver
    // // Implicit_LBFGS_Solver solver = Implicit_LBFGS_Solver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    // timestamp_t t0 = get_timestamp();
    // solver.run(cell_mesh, 5);
    // timestamp_t t1 = get_timestamp();
    // double secs = (t1 - t0) / 1000000.0L;

    // std::cout << "Runtime: " << secs << " Seconds\n";

    //////////////////////////////////////////////////
    //Un-comment the following code block (and re-compile) to run the large voronoi mesh example
    //////////////////////////////////////////////////

    // cell_mesh.importMesh("init_voronoi_large.vtk");
    // std::vector<double> cell_areas = cell_mesh.cellAreas();
    // cell_mesh.setA_0(cell_areas);

    // //Un-comment to use the explicit solver
    // // ExplicitSolver solver = ExplicitSolver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    // //Un-comment to use the RK4 solver
    // // RungeKutta4Solver solver = RungeKutta4Solver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    // //Un-comment to use the L-BFGS solver
    // Implicit_LBFGS_Solver solver = Implicit_LBFGS_Solver(cell_mesh, dt, write_dt, use_parallel, apply_bc);(cell_mesh, dt, write_dt, use_parallel);

    // timestamp_t t0 = get_timestamp();
    // solver.run(cell_mesh, 1);
    // timestamp_t t1 = get_timestamp();
    // double secs = (t1 - t0) / 1000000.0L;

    // std::cout << "Runtime: " << secs << " Seconds\n";

    //////////////////////////////////////////////////
    //Un-comment the following code block (and re-compile) to run the small cavity mesh example
    //////////////////////////////////////////////////

    cell_mesh.importMesh("init_sphere_small.vtk");
    cell_mesh.randomize(1e-2);
    std::vector<double> cell_areas = cell_mesh.cellAreas();
    cell_mesh.setA_0(cell_areas);
    cell_mesh.setYolk_id(0);

    //Un-comment to use the explicit solver
    // ExplicitSolver solver = ExplicitSolver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    //Un-comment to use the RK4 solver
    // RungeKutta4Solver solver = RungeKutta4Solver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    //Un-comment to use the L-BFGS solver
    Implicit_LBFGS_Solver solver = Implicit_LBFGS_Solver(cell_mesh, dt, write_dt, use_parallel, apply_bc);

    timestamp_t t0 = get_timestamp();
    solver.run(cell_mesh, 40);
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;

    std::cout << "Runtime: " << secs << " Seconds\n";

    return 0;
}