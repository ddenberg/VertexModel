# VertexModel

## To compile:

mkdir build

cd build

cmake -DCMAKE_BUILD_TYPE=Release ..

make

### After initial compile:

Copy the mesh files in the folder 'meshes' to the build directory

Next create a folder in the build directory called 'output'

## To run the small voronoi mesh:

Uncomment the relevant section of code in 'main.cpp' (code sections are highlighted within the file)

Set apply_bc = 0

Uncomment the solver you wish to use within the code block and change dt and write_dt as you see fit.

Set use_parallel to 1 to run the multithreaded version

Recompile the code after setting your paramters and run ./VertexModel

## To run the large voronoi mesh:

Uncomment the relevant section of code in 'main.cpp' (code sections are highlighted within the file)

Set apply_bc = 0

Uncomment the solver you wish to use within the code block and change dt and write_dt as you see fit.

Set use_parallel to 1 to run the multithreaded version

Recompile the code after setting your paramters and run ./VertexModel

## To run the small cavity mesh:

Uncomment the relevant section of code in 'main.cpp' (code sections are highlighted within the file)

Set apply_bc = 1

Uncomment the solver you wish to use within the code block and change dt and write_dt as you see fit.

Set use_parallel to 1 to run the multithreaded version

Recompile the code after setting your paramters and run ./VertexModel

## For Visualization:

Use the program Paraview (https://www.paraview.org/) which can visualize .VTK files.
