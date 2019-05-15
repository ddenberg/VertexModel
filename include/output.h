#pragma once

#include <string>

void write_vtk(std::string filename) {    
    int cell_size = 0;
    for (int i = 1; i <= Nc; i++) {
        cell_size += basal_vertices[i][2] + 1;
    }
    
    //BASAL SIDES
    snprintf(filename, sizeof(char) * 50, "./output/outBasal_%d.vtk", kcount);
    FILE *file2; 
    file2 = fopen(filename, "w");
    fprintf(file2, "# vtk DataFile Version 3.0\n");
    fprintf(file2, "Cell Mesh\n");
    fprintf(file2, "ASCII\n\n");
    fprintf(file2, "DATASET POLYDATA\n");
    fprintf(file2, "POINTS %d double\n", Nv);
    //Write out verts
    for(int i = 1; i <= Nv; i++) {
        // if (v[i][0] > 0.5) {
            fprintf(file2, "%f %f %f\n", v[i][1], v[i][2], v[i][3]);
        // }
    }
    fprintf(file2, "\nPOLYGONS\t%d\t%d\n", Nc, cell_size);
    for (int i = 1; i <= Nc; i++) {
        fprintf(file2, "%d ", basal_vertices[i][2]);
        for(size_t j = 3; j <= 2+basal_edges[i][2]; ++j) {
            fprintf(file2, "%d ", basal_vertices[i][j]-1);
        }
        fprintf(file2, "\n");
    }
    fclose(file2);
    
    //APICAL SIDES
    snprintf(filename, sizeof(char) * 50, "./output/outApical_%d.vtk", kcount);
    FILE *file3;
    file3 = fopen(filename, "w");
    fprintf(file3, "# vtk DataFile Version 3.0\n");
    fprintf(file3, "Cell Mesh\n");
    fprintf(file3, "ASCII\n\n");
    fprintf(file3, "DATASET POLYDATA\n");
    fprintf(file3, "POINTS %d double\n", Nv);
    //Write out verts
    for(int i = 1; i <= Nv; i++) {
        // if (v[i][0] > 0.5) {
            fprintf(file3, "%f %f %f\n", v[i][1], v[i][2], v[i][3]);
        // }
    }
    fprintf(file3, "\nPOLYGONS\t%d\t%d\n", Nc, cell_size);
    for(int i = 1; i <= Nc; i++){
        fprintf(file3, "%d ", basal_vertices[i][2]);
        for(size_t j = 3; j <= 2+basal_edges[i][2]; ++j) {
            fprintf(file3, "%d ", v_partner[basal_vertices[i][j]]-1);
        }
        // fprintf(file3, "%d\n", v_partner[basal_vertices[i][2+basal_edges[i][2]]]);
        fprintf(file3, "\n");
    }
    fclose(file3);
    
    //LATERAL SIDES
    snprintf(filename, sizeof(char) * 50, "./output/outLateral_%d.vtk", kcount);
    FILE *file4;
    file4 = fopen(filename, "w");
    fprintf(file4, "# vtk DataFile Version 3.0\n");
    fprintf(file4, "Cell Mesh\n");
    fprintf(file4, "ASCII\n\n");
    fprintf(file4, "DATASET POLYDATA\n");
    fprintf(file4, "POINTS %d double\n", Nv);
    //Write out verts
    for(int i = 1; i <= Nv; i++) {
        // if (v[i][0] > 0.5) {
            fprintf(file3, "%f %f %f\n", v[i][1], v[i][2], v[i][3]);
        // }
    }
    fprintf(file4, "\nPOLYGONS\t%d\t%d\n", Ne, 5*Ne);
    for(int i = 1; i <= Ne; i++) {
        fprintf(file4, "4 %d %d %d %d\n", e[i][1]-1, e[i][2]-1, v_partner[e[i][2]]-1, v_partner[e[i][1]]-1);
    }
    fclose(file4);
}