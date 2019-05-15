#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "vec3.h"
#include "util.h"
#include "largevec.h"
#include <random>

class Mesh {
public:
    Mesh();
    ~Mesh();
    
    int importMesh(std::string filename);
    int writeMesh(std::string filename);
    int numVerts() { return r.size(); }
    int numCells() { return cells.size(); }
    std::vector<double> cellAreas();
    void setA_0(const std::vector<double> &A) { A_0 = A; }
    void setA_0_i(double A, double i) { A_0[i] = A; }
    void randomize(double factor);
    void setYolk_id(int i) { yolk_id = i; }
    
    std::vector<std::vector<int>> cells;
    std::vector<Vec3> r;
    std::vector<double> A_0;
    int yolk_id = -1;
    std::vector<int> vertPerCell;
    std::vector<int> vertPerCell_cumul;
};