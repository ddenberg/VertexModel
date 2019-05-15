#include "mesh.h"

Mesh::Mesh() {};
Mesh::~Mesh() {};

int Mesh::importMesh(std::string filename) {
    std::string line;
    std::ifstream file(filename);
    if (file.is_open()) {
        std::getline(file, line);

        //Check to make sure file has vtk header
        if (line.compare(2, 3, "vtk") != 0) {
            std::cout << "The file: " << filename << " Is not a valid .vtk formated file." << std::endl;
        }

        //Check to make sure there includes a POINTS section
        while (line.compare(0, 6, "POINTS") != 0) {
            std::getline(file, line);
        }
        int nVerts;
        sscanf(line.c_str(), "%*s %d %*s", &nVerts); //Read number of vertices
        r.resize(nVerts);

        //Read in list of vertices
        double x, y, z;
        for (int i = 0; i < nVerts; i++) {
            std::getline(file, line);
            sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
            r[i] = Vec3(x, y, z);
        }

        //Check to make sure there is a POLYGONS section
        while (line.compare(0, 8, "POLYGONS") != 0) {
            std::getline(file, line);
        }
        int nCells;
        sscanf(line.c_str(), "%*s %d %*d", &nCells); //Read in number of cells
        cells.resize(nCells);

        //Read in list of cells
        std::vector<int> temp;
        for (int i = 0; i < nCells; i++) {
            std::getline(file, line);
            temp = string_to_vecInt(line);
            temp.erase(temp.begin());
            cells[i] = temp;
        }
        
        file.close();

        std::cout << "Read in VTK file: " << filename << std::endl;
    } else {
        std::cout << "Unable to read file: " << filename << std::endl;
        return -1;
    }
    vertPerCell.resize(cells.size());
    for (int i = 0; i < cells.size(); i++) {
        vertPerCell[i] = cells[i].size();
    }
    vertPerCell_cumul.resize(cells.size()+1);
    vertPerCell_cumul[0] = 0;
    for (int i = 1; i < vertPerCell_cumul.size(); i++) {
        vertPerCell_cumul[i] = vertPerCell_cumul[i-1] + vertPerCell[i-1];
    }

    return 0;
}

int Mesh::writeMesh(std::string filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        //Write vtk header
        file << "# vtk DataFile Version 3.0\n";
        file << "Cell Mesh\n";
        file << "ASCII\n\n";
        file << "DATASET POLYDATA\n";

        //Write out array of vertices
        int nV = numVerts();
        file << "POINTS " << nV << " double\n";
        file << std::setprecision(12);
        for (int i = 0; i < nV; i++) {
            file << r[i].x() << " " << r[i].y() << " " << r[i].z() << "\n";
        }

        //Calculate number of entries in cell list for VTK formatting
        int cell_size = 0;
        int nC = numCells();
        for (int i = 0; i < nC; i++) {
            cell_size += cells[i].size() + 1;
        }
        
        //Write out cell list
        file << "POLYGONS " << nC << " " << cell_size << "\n";
        for (int i = 0; i < nC; i++) {
            file << cells[i].size() << " ";
            for (int j = 0; j < cells[i].size()-1; j++) {
                file << cells[i][j] << " ";
            }
            file << cells[i][cells[i].size()-1] << "\n";
        }

        file.close();
    } else {
        std::cout << "Unable to write to file: " << filename << std::endl;
        return -1;
    }
    return 0;
}

std::vector<double> Mesh::cellAreas() {
    std::vector<double> A(cells.size());
    for (int c = 0; c < cells.size(); c++) {
        double A_temp = 0.;
        int n = cells[c].size();
        for (int i = 0; i < n; i++) {
            int ind1 = (i + 1) % n;
            int ind2 = i == 0 ? n - 1 : i - 1;
            A_temp += r[cells[c][i]].x() * (r[cells[c][ind1]].y() - r[cells[c][ind2]].y());
        }
        A[c] = 0.5 * std::abs(A_temp);
    }
    return A;
}

void Mesh::randomize(double factor) {
    std::default_random_engine gen;
    std::normal_distribution<double> dist(0.0,1.0);

    for (int i = 0; i < r.size(); i++) {
        Vec3 p = Vec3(dist(gen), dist(gen), 0.0);
        r[i] += factor * p;
    }
}