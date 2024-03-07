elem = readmatrix("epeltr4.dat");
elem = elem(:, 2:5);

node = readmatrix("npeltr4.dat");
node = node(:, 2:3);

LHS_lapack = readmatrix("output/LHS_lapack.dat");
LHS_unpacked = readmatrix("output/LHS_unpacked.dat");
RHS = readmatrix("output/RHS.dat");

u_lapack = readmatrix("output/u_lapack.dat");
u_dsolve = readmatrix("output/u_dsolve.dat");

close all;
FEM_patch(elem, node, u_dsolve);
FEM_patch(elem, node, u_lapack);
