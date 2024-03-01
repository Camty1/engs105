LHS = readmatrix("output/LHS.dat");
LHS_lapack = readmatrix("output/LHS_lapack.dat");

sum(LHS ~= LHS_lapack(end-size(LHS,2) + 1, :)')
