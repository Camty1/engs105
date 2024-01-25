%% ENGS 105 HW 2
% Cameron Wolfe 1/19/24

%% Part 1 - Type 3 boundary condition
% Input data files
x_50_1 = readmatrix_fortran("output/1/x050.dat");
y_50_1 = readmatrix_fortran("output/1/y050.dat");
u_50_1 = readmatrix_fortran("output/1/u050.dat");

x_50_10 = readmatrix_fortran("output/10/x050.dat");
y_50_10 = readmatrix_fortran("output/10/y050.dat");
u_50_10 = readmatrix_fortran("output/10/u050.dat");

x_50_100 = readmatrix_fortran("output/100/x050.dat");
y_50_100 = readmatrix_fortran("output/100/y050.dat");
u_50_100 = readmatrix_fortran("output/100/u050.dat");

x_50_01 = readmatrix_fortran("output/01/x050.dat");
y_50_01 = readmatrix_fortran("output/01/y050.dat");
u_50_01 = readmatrix_fortran("output/01/u050.dat");

x_50_001 = readmatrix_fortran("output/001/x050.dat");
y_50_001 = readmatrix_fortran("output/001/y050.dat");
u_50_001 = readmatrix_fortran("output/001/u050.dat");

% Form meshgrids of data
x_50_1_meshgrid = reshape(x_50_1, 50, 50);
y_50_1_meshgrid = reshape(y_50_1, 50, 50);
u_50_1_meshgrid = reshape(u_50_1, 50, 50);

x_50_10_meshgrid = reshape(x_50_10, 50, 50);
y_50_10_meshgrid = reshape(y_50_10, 50, 50);
u_50_10_meshgrid = reshape(u_50_10, 50, 50);

x_50_100_meshgrid = reshape(x_50_100, 50, 50);
y_50_100_meshgrid = reshape(y_50_100, 50, 50);
u_50_100_meshgrid = reshape(u_50_100, 50, 50);

x_50_01_meshgrid = reshape(x_50_01, 50, 50);
y_50_01_meshgrid = reshape(y_50_01, 50, 50);
u_50_01_meshgrid = reshape(u_50_01, 50, 50);

x_50_001_meshgrid = reshape(x_50_001, 50, 50);
y_50_001_meshgrid = reshape(y_50_001, 50, 50);
u_50_001_meshgrid = reshape(u_50_001, 50, 50);

% Use common limits for all plots to make comparison easier
plotting_limits = [0, 1; 0, 1; -1.5, 1];

%% Plotting for part 1
set(0,'DefaultFigureWindowStyle','docked')
figure(1);
colormap winter;
surf(x_50_1_meshgrid, y_50_1_meshgrid, u_50_1_meshgrid, "EdgeAlpha", 0.3);
view(15,30);
title("Numerical Solution of PDE, N=50, A=1");
xlabel("x");
ylabel("y");
zlabel("Potential");
xlim(plotting_limits(1,:));
ylim(plotting_limits(2,:));
zlim(plotting_limits(3,:));
saveas(gcf, "figures/1.png")

figure(2);
colormap winter;
surf(x_50_10_meshgrid, y_50_10_meshgrid, u_50_10_meshgrid, "EdgeAlpha", 0.3);
view(15,30);
title("Numerical Solution of PDE, N=50, A=10");
xlabel("x");
ylabel("y");
zlabel("Potential");
xlim(plotting_limits(1,:));
ylim(plotting_limits(2,:));
zlim(plotting_limits(3,:));
saveas(gcf, "figures/10.png")

figure(3);
colormap winter;
surf(x_50_100_meshgrid, y_50_100_meshgrid, u_50_100_meshgrid, "EdgeAlpha", 0.3);
view(15,30);
title("Numerical Solution of PDE, N=50, A=100");
xlabel("x");
ylabel("y");
zlabel("Potential");
xlim(plotting_limits(1,:));
ylim(plotting_limits(2,:));
zlim(plotting_limits(3,:));
saveas(gcf, "figures/100.png")

figure(4);
colormap winter;
surf(x_50_01_meshgrid, y_50_01_meshgrid, u_50_01_meshgrid, "EdgeAlpha", 0.3);
view(15,30);
title("Numerical Solution of PDE, N=50, A=0.1");
xlabel("x");
ylabel("y");
zlabel("Potential");
xlim(plotting_limits(1,:));
ylim(plotting_limits(2,:));
zlim(plotting_limits(3,:));
saveas(gcf, "figures/01.png")

figure(5);
colormap winter;
surf(x_50_001_meshgrid, y_50_001_meshgrid, u_50_001_meshgrid, "EdgeAlpha", 0.3);
view(15,30);
title("Numerical Solution of PDE, N=50, A=0.01");
xlabel("x");
ylabel("y");
xlim(plotting_limits(1,:));
ylim(plotting_limits(2,:));
zlim(plotting_limits(3,:));
zlabel("Potential");
saveas(gcf, "figures/001.png")

%% Part 2 - Iterative solvers
% Read in matrices
x_10 = readmatrix_fortran('output/x010.dat');
y_10 = readmatrix_fortran('output/y010.dat');

x_20 = readmatrix_fortran('output/x020.dat');
y_20 = readmatrix_fortran('output/y020.dat');

x_40 = readmatrix_fortran('output/x040.dat');
y_40 = readmatrix_fortran('output/y040.dat');

err_jac_10_1 = readmatrix_fortran('output/err_jac010_1.0E+00.dat');
err_jac_10_100 = readmatrix_fortran('output/err_jac010_1.0E+02.dat');
err_jac_10_001 = readmatrix_fortran('output/err_jac010_1.0E-02.dat');

err_jac_20_1 = readmatrix_fortran('output/err_jac020_1.0E+00.dat');
err_jac_20_100 = readmatrix_fortran('output/err_jac020_1.0E+02.dat');
err_jac_20_001 = readmatrix_fortran('output/err_jac020_1.0E-02.dat');

err_jac_40_1 = readmatrix_fortran('output/err_jac040_1.0E+00.dat');
err_jac_40_100 = readmatrix_fortran('output/err_jac040_1.0E+02.dat');
err_jac_40_001 = readmatrix_fortran('output/err_jac040_1.0E-02.dat');

err_gs_10_1 = readmatrix_fortran('output/err_gs010_1.0E+00.dat');
err_gs_10_100 = readmatrix_fortran('output/err_gs010_1.0E+02.dat');
err_gs_10_001 = readmatrix_fortran('output/err_gs010_1.0E-02.dat');

err_gs_20_1 = readmatrix_fortran('output/err_gs020_1.0E+00.dat');
err_gs_20_100 = readmatrix_fortran('output/err_gs020_1.0E+02.dat');
err_gs_20_001 = readmatrix_fortran('output/err_gs020_1.0E-02.dat');

err_gs_40_1 = readmatrix_fortran('output/err_gs040_1.0E+00.dat');
err_gs_40_100 = readmatrix_fortran('output/err_gs040_1.0E+02.dat');
err_gs_40_001 = readmatrix_fortran('output/err_gs040_1.0E-02.dat');

err_sor_10_1 = readmatrix_fortran('output/err_sor010_1.0E+00.dat');
err_sor_10_100 = readmatrix_fortran('output/err_sor010_1.0E+02.dat');
err_sor_10_001 = readmatrix_fortran('output/err_sor010_1.0E-02.dat');

err_sor_20_1 = readmatrix_fortran('output/err_sor020_1.0E+00.dat');
err_sor_20_100 = readmatrix_fortran('output/err_sor020_1.0E+02.dat');
err_sor_20_001 = readmatrix_fortran('output/err_sor020_1.0E-02.dat');

err_sor_40_1 = readmatrix_fortran('output/err_sor040_1.0E+00.dat');
err_sor_40_100 = readmatrix_fortran('output/err_sor040_1.0E+02.dat');
err_sor_40_001 = readmatrix_fortran('output/err_sor040_1.0E-02.dat');

u_10_1 = readmatrix_fortran('output/u010_1.0E+00.dat');
u_10_100 = readmatrix_fortran('output/u010_1.0E+02.dat');
u_10_001 = readmatrix_fortran('output/u010_1.0E-02.dat');

u_20_1 = readmatrix_fortran('output/u020_1.0E+00.dat');
u_20_100 = readmatrix_fortran('output/u020_1.0E+02.dat');
u_20_001 = readmatrix_fortran('output/u020_1.0E-02.dat');

u_40_1 = readmatrix_fortran('output/u040_1.0E+00.dat');
u_40_100 = readmatrix_fortran('output/u040_1.0E+02.dat');
u_40_001 = readmatrix_fortran('output/u040_1.0E-02.dat');

u_jac_10_1 = readmatrix_fortran('output/u_jac010_1.0E+00.dat');
u_jac_10_100 = readmatrix_fortran('output/u_jac010_1.0E+02.dat');
u_jac_10_001 = readmatrix_fortran('output/u_jac010_1.0E-02.dat');

u_jac_20_1 = readmatrix_fortran('output/u_jac020_1.0E+00.dat');
u_jac_20_100 = readmatrix_fortran('output/u_jac020_1.0E+02.dat');
u_jac_20_001 = readmatrix_fortran('output/u_jac020_1.0E-02.dat');

u_jac_40_1 = readmatrix_fortran('output/u_jac040_1.0E+00.dat');
u_jac_40_100 = readmatrix_fortran('output/u_jac040_1.0E+02.dat');
u_jac_40_001 = readmatrix_fortran('output/u_jac040_1.0E-02.dat');

u_gs_10_1 = readmatrix_fortran('output/u_gs010_1.0E+00.dat');
u_gs_10_100 = readmatrix_fortran('output/u_gs010_1.0E+02.dat');
u_gs_10_001 = readmatrix_fortran('output/u_gs010_1.0E-02.dat');

u_gs_20_1 = readmatrix_fortran('output/u_gs020_1.0E+00.dat');
u_gs_20_100 = readmatrix_fortran('output/u_gs020_1.0E+02.dat');
u_gs_20_001 = readmatrix_fortran('output/u_gs020_1.0E-02.dat');

u_gs_40_1 = readmatrix_fortran('output/u_gs040_1.0E+00.dat');
u_gs_40_100 = readmatrix_fortran('output/u_gs040_1.0E+02.dat');
u_gs_40_001 = readmatrix_fortran('output/u_gs040_1.0E-02.dat');

u_sor_10_1 = readmatrix_fortran('output/u_sor010_1.0E+00.dat');
u_sor_10_100 = readmatrix_fortran('output/u_sor010_1.0E+02.dat');
u_sor_10_001 = readmatrix_fortran('output/u_sor010_1.0E-02.dat');

u_sor_20_1 = readmatrix_fortran('output/u_sor020_1.0E+00.dat');
u_sor_20_100 = readmatrix_fortran('output/u_sor020_1.0E+02.dat');
u_sor_20_001 = readmatrix_fortran('output/u_sor020_1.0E-02.dat');

u_sor_40_1 = readmatrix_fortran('output/u_sor040_1.0E+00.dat');
u_sor_40_100 = readmatrix_fortran('output/u_sor040_1.0E+02.dat');
u_sor_40_001 = readmatrix_fortran('output/u_sor040_1.0E-02.dat');

x_10_meshgrid = reshape(x_10, 10, 10);
y_10_meshgrid = reshape(y_10, 10, 10);

x_20_meshgrid = reshape(x_20, 20, 20);
y_20_meshgrid = reshape(y_20, 20, 20);

x_40_meshgrid = reshape(x_40, 40, 40);
y_40_meshgrid = reshape(y_40, 40, 40);

u_10_1_meshgrid = reshape(u_10_1, 10, 10);
u_10_100_meshgrid = reshape(u_10_100, 10, 10);
u_10_001_meshgrid = reshape(u_10_001, 10, 10);

u_20_1_meshgrid = reshape(u_20_1, 20, 20);
u_20_100_meshgrid = reshape(u_20_100, 20, 20);
u_20_001_meshgrid = reshape(u_20_001, 20, 20);

u_40_1_meshgrid = reshape(u_40_1, 40, 40);
u_40_100_meshgrid = reshape(u_40_100, 40, 40);
u_40_001_meshgrid = reshape(u_40_001, 40, 40);

u_jac_10_1_meshgrid = reshape(u_jac_10_1, 10, 10);
u_jac_10_100_meshgrid = reshape(u_jac_10_100, 10, 10);
u_jac_10_001_meshgrid = reshape(u_jac_10_001, 10, 10);

u_jac_20_1_meshgrid = reshape(u_jac_20_1, 20, 20);
u_jac_20_100_meshgrid = reshape(u_jac_20_100, 20, 20);
u_jac_20_001_meshgrid = reshape(u_jac_20_001, 20, 20);

u_jac_40_1_meshgrid = reshape(u_jac_40_1, 40, 40);
u_jac_40_100_meshgrid = reshape(u_jac_40_100, 40, 40);
u_jac_40_001_meshgrid = reshape(u_jac_40_001, 40, 40);

u_gs_10_1_meshgrid = reshape(u_gs_10_1, 10, 10);
u_gs_10_100_meshgrid = reshape(u_gs_10_100, 10, 10);
u_gs_10_001_meshgrid = reshape(u_gs_10_001, 10, 10);

u_gs_20_1_meshgrid = reshape(u_gs_20_1, 20, 20);
u_gs_20_100_meshgrid = reshape(u_gs_20_100, 20, 20);
u_gs_20_001_meshgrid = reshape(u_gs_20_001, 20, 20);

u_gs_40_1_meshgrid = reshape(u_gs_40_1, 40, 40);
u_gs_40_100_meshgrid = reshape(u_gs_40_100, 40, 40);
u_gs_40_001_meshgrid = reshape(u_gs_40_001, 40, 40);

u_sor_10_1_meshgrid = reshape(u_sor_10_1, 10, 10);
u_sor_10_100_meshgrid = reshape(u_sor_10_100, 10, 10);
u_sor_10_001_meshgrid = reshape(u_sor_10_001, 10, 10);

u_sor_20_1_meshgrid = reshape(u_sor_20_1, 20, 20);
u_sor_20_100_meshgrid = reshape(u_sor_20_100, 20, 20);
u_sor_20_001_meshgrid = reshape(u_sor_20_001, 20, 20);

u_sor_40_1_meshgrid = reshape(u_sor_40_1, 40, 40);
u_sor_40_100_meshgrid = reshape(u_sor_40_100, 40, 40);
u_sor_40_001_meshgrid = reshape(u_sor_40_001, 40, 40);

e_jac = zeros(1,9);

e_jac_10_1_meshgrid = error_map_mat(u_jac_10_1_meshgrid, u_10_1_meshgrid);
e_jac_10_100_meshgrid = error_map_mat(u_jac_10_100_meshgrid, u_10_100_meshgrid);
e_jac_10_001_meshgrid = error_map_mat(u_jac_10_001_meshgrid, u_10_001_meshgrid);
e_jac(1) = rms_error(e_jac_10_1_meshgrid);
e_jac(2) = rms_error(e_jac_10_100_meshgrid);
e_jac(3) = rms_error(e_jac_10_001_meshgrid);

e_jac_20_1_meshgrid = error_map_mat(u_jac_20_1_meshgrid, u_20_1_meshgrid);
e_jac_20_100_meshgrid = error_map_mat(u_jac_20_100_meshgrid, u_20_100_meshgrid);
e_jac_20_001_meshgrid = error_map_mat(u_jac_20_001_meshgrid, u_20_001_meshgrid);
e_jac(4) = rms_error(e_jac_20_1_meshgrid);
e_jac(5) = rms_error(e_jac_20_100_meshgrid);
e_jac(6) = rms_error(e_jac_20_001_meshgrid);

e_jac_40_1_meshgrid = error_map_mat(u_jac_40_1_meshgrid, u_40_1_meshgrid);
e_jac_40_100_meshgrid = error_map_mat(u_jac_40_100_meshgrid, u_40_100_meshgrid);
e_jac_40_001_meshgrid = error_map_mat(u_jac_40_001_meshgrid, u_40_001_meshgrid);
e_jac(7) = rms_error(e_jac_40_1_meshgrid);
e_jac(8) = rms_error(e_jac_40_100_meshgrid);
e_jac(9) = rms_error(e_jac_40_001_meshgrid);

e_gs = zeros(1,9);

e_gs_10_1_meshgrid = error_map_mat(u_gs_10_1_meshgrid, u_10_1_meshgrid);
e_gs_10_100_meshgrid = error_map_mat(u_gs_10_100_meshgrid, u_10_100_meshgrid);
e_gs_10_001_meshgrid = error_map_mat(u_gs_10_001_meshgrid, u_10_001_meshgrid);
e_gs(1) = rms_error(e_gs_10_1_meshgrid);
e_gs(2) = rms_error(e_gs_10_100_meshgrid);
e_gs(3) = rms_error(e_gs_10_001_meshgrid);

e_gs_20_1_meshgrid = error_map_mat(u_gs_20_1_meshgrid, u_20_1_meshgrid);
e_gs_20_100_meshgrid = error_map_mat(u_gs_20_100_meshgrid, u_20_100_meshgrid);
e_gs_20_001_meshgrid = error_map_mat(u_gs_20_001_meshgrid, u_20_001_meshgrid);
e_gs(4) = rms_error(e_gs_20_1_meshgrid);
e_gs(5) = rms_error(e_gs_20_100_meshgrid);
e_gs(6) = rms_error(e_gs_20_001_meshgrid);

e_gs_40_1_meshgrid = error_map_mat(u_gs_40_1_meshgrid, u_40_1_meshgrid);
e_gs_40_100_meshgrid = error_map_mat(u_gs_40_100_meshgrid, u_40_100_meshgrid);
e_gs_40_001_meshgrid = error_map_mat(u_gs_40_001_meshgrid, u_40_001_meshgrid);
e_gs(7) = rms_error(e_gs_40_1_meshgrid);
e_gs(8) = rms_error(e_gs_40_100_meshgrid);
e_gs(9) = rms_error(e_gs_40_001_meshgrid);

e_sor = zeros(1,9);

e_sor_10_1_meshgrid = error_map_mat(u_sor_10_1_meshgrid, u_10_1_meshgrid);
e_sor_10_100_meshgrid = error_map_mat(u_sor_10_100_meshgrid, u_10_100_meshgrid);
e_sor_10_001_meshgrid = error_map_mat(u_sor_10_001_meshgrid, u_10_001_meshgrid);
e_sor(1) = rms_error(e_sor_10_1_meshgrid);
e_sor(2) = rms_error(e_sor_10_100_meshgrid);
e_sor(3) = rms_error(e_sor_10_001_meshgrid);

e_sor_20_1_meshgrid = error_map_mat(u_sor_20_1_meshgrid, u_20_1_meshgrid);
e_sor_20_100_meshgrid = error_map_mat(u_sor_20_100_meshgrid, u_20_100_meshgrid);
e_sor_20_001_meshgrid = error_map_mat(u_sor_20_001_meshgrid, u_20_001_meshgrid);
e_sor(4) = rms_error(e_sor_20_1_meshgrid);
e_sor(5) = rms_error(e_sor_20_100_meshgrid);
e_sor(6) = rms_error(e_sor_20_001_meshgrid);

e_sor_40_1_meshgrid = error_map_mat(u_sor_40_1_meshgrid, u_40_1_meshgrid);
e_sor_40_100_meshgrid = error_map_mat(u_sor_40_100_meshgrid, u_40_100_meshgrid);
e_sor_40_001_meshgrid = error_map_mat(u_sor_40_001_meshgrid, u_40_001_meshgrid);
e_sor(7) = rms_error(e_sor_40_1_meshgrid);
e_sor(8) = rms_error(e_sor_40_100_meshgrid);
e_sor(9) = rms_error(e_sor_40_001_meshgrid);

% Plotting for Part 2
figure(6);
tcl = tiledlayout(3,3);

nexttile;
semilogy(0:length(err_jac_10_001)-1, err_jac_10_001);
hold on;
semilogy(0:length(err_gs_10_001)-1, err_gs_10_001);
semilogy(0:length(err_sor_10_001)-1, err_sor_10_001);
hold off;
title("N = 10, A = 0.01");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_10_1)-1, err_jac_10_1);
hold on;
semilogy(0:length(err_gs_10_1)-1, err_gs_10_1);
semilogy(0:length(err_sor_10_1)-1, err_sor_10_1);
hold off;
title("N = 10, A = 1");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_10_100)-1, err_jac_10_100);
hold on;
semilogy(0:length(err_gs_10_100)-1, err_gs_10_100);
semilogy(0:length(err_sor_10_100)-1, err_sor_10_100);
hold off;
title("N = 10, A = 100");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_20_001)-1, err_jac_20_001);
hold on;
semilogy(0:length(err_gs_20_001)-1, err_gs_20_001);
semilogy(0:length(err_sor_20_001)-1, err_sor_20_001);
hold off;
title("N = 20, A = 0.01");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_20_1)-1, err_jac_20_1);
hold on;
semilogy(0:length(err_gs_20_1)-1, err_gs_20_1);
semilogy(0:length(err_sor_20_1)-1, err_sor_20_1);
hold off;
title("N = 20, A = 1");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_20_100)-1, err_jac_20_100);
hold on;
semilogy(0:length(err_gs_20_100)-1, err_gs_20_100);
semilogy(0:length(err_sor_20_100)-1, err_sor_20_100);
hold off;
title("N = 20, A = 100");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_40_001)-1, err_jac_40_001);
hold on;
semilogy(0:length(err_gs_40_001)-1, err_gs_40_001);
semilogy(0:length(err_sor_40_001)-1, err_sor_40_001);
hold off;
title("N = 40, A = 0.01");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_40_1)-1, err_jac_40_1);
hold on;
semilogy(0:length(err_gs_40_1)-1, err_gs_40_1);
semilogy(0:length(err_sor_40_1)-1, err_sor_40_1);
hold off;
title("N = 40, A = 1");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

nexttile;
semilogy(0:length(err_jac_40_100)-1, err_jac_40_100);
hold on;
semilogy(0:length(err_gs_40_100)-1, err_gs_40_100);
semilogy(0:length(err_sor_40_100)-1, err_sor_40_100);
hold off;
title("N = 40, A = 100");
xlabel("Iteration");
ylabel("L\infty Error");
legend("Jacobi", "Gauss Sidel", "SOR");

title(tcl, "L\infty Error vs Iteration for N = 10, 20, 40, and A = 0.01, 1, 100");
saveas(gcf, 'figures/errorvsiteration.png');

N_values = [10 20 40];
A_values = [0.01 1 100];

figure(7);
tcl = tiledlayout(3,3);

% N = 10
nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_jac_10_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 0.01, rms = %f", e_jac(1)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_jac_10_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 1, rms = %f", e_jac(2)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_jac_10_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 100, rms = %f", e_jac(3)));
xlabel("x")
ylabel("y");
colorbar;

% N = 20
nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_jac_20_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 0.01, rms = %f", e_jac(4)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_jac_20_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 1, rms = %f", e_jac(5)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_jac_20_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 100, rms = %f", e_jac(6)));
xlabel("x")
ylabel("y");
colorbar;

% N = 40
nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_jac_40_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 0.01, rms = %f", e_jac(7)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_jac_40_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 1, rms = %f", e_jac(8)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_jac_40_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 1, rms = %f", e_jac(9)));
xlabel("x")
ylabel("y");
colorbar;

title(tcl, "Jacobi Error Contours for N = 10, 20, 40, and A = 0.001, 1, 100");
saveas(gcf, "figures/e_jac.png");

figure(8);
tcl = tiledlayout(3,3);

% N = 10
nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_gs_10_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 0.01, rms = %f", e_gs(1)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_gs_10_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 1, rms = %f", e_gs(2)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_gs_10_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 100, rms = %f", e_gs(3)));
xlabel("x")
ylabel("y");
colorbar;

% N = 20
nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_gs_20_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 0.01, rms = %f", e_gs(4)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_gs_20_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 1, rms = %f", e_gs(5)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_gs_20_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 100, rms = %f", e_gs(6)));
xlabel("x")
ylabel("y");
colorbar;

% N = 40
nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_gs_40_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 0.01, rms = %f", e_gs(7)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_gs_40_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 1, rms = %f", e_gs(8)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_gs_40_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 100, rms = %f", e_gs(9)));
xlabel("x")
ylabel("y");
colorbar;

title(tcl, "Gauss Sidel Error Contours for N = 10, 20, 40, and A = 0.001, 1, 100");
saveas(gcf, "figures/e_gs.png");

figure(9);
tcl = tiledlayout(3,3);

% N = 10
nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_sor_10_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 0.01, rms = %f", e_sor(1)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_sor_10_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 1, rms = %f", e_sor(2)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_10_meshgrid, y_10_meshgrid, (e_sor_10_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 10, A = 100, rms = %f", e_sor(3)));
xlabel("x")
ylabel("y");
colorbar;

% N = 20
nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_sor_20_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 0.01, rms = %f", e_sor(4)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_sor_20_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 1, rms = %f", e_sor(5)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_20_meshgrid, y_20_meshgrid, (e_sor_20_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 20, A = 100, rms = %f", e_sor(6)));
xlabel("x")
ylabel("y");
colorbar;

% N = 40
nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_sor_40_001_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 0.01, rms = %f", e_sor(7)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_sor_40_1_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 1, rms = %f", e_sor(8)));
xlabel("x")
ylabel("y");
colorbar;

nexttile;
colormap winter;
contourf(x_40_meshgrid, y_40_meshgrid, (e_sor_40_100_meshgrid), 20, 'EdgeAlpha', 0.3);
title(sprintf("N = 40, A = 100, rms = %f", e_sor(9)));
xlabel("x")
ylabel("y");
colorbar;

title(tcl, "SOR Error Contours for N = 10, 20, 40, and A = 0.001, 1, 100");
saveas(gcf, "figures/e_sor.png");

%% Spectral Radius
rho_jac_10_1_arr = err_jac_10_1(2:end) ./ err_jac_10_1(1:end-1);
rho_jac_10_100_arr = err_jac_10_100(2:end) ./ err_jac_10_100(1:end-1);
rho_jac_10_001_arr = err_jac_10_001(2:end) ./ err_jac_10_001(1:end-1);

rho_jac_20_1_arr = err_jac_20_1(2:end) ./ err_jac_20_1(1:end-1);
rho_jac_20_100_arr = err_jac_20_100(2:end) ./ err_jac_20_100(1:end-1);
rho_jac_20_001_arr = err_jac_20_001(2:end) ./ err_jac_20_001(1:end-1);

rho_jac_40_1_arr = err_jac_40_1(2:end) ./ err_jac_40_1(1:end-1);
rho_jac_40_100_arr = err_jac_40_100(2:end) ./ err_jac_40_100(1:end-1);
rho_jac_40_001_arr = err_jac_40_001(2:end) ./ err_jac_40_001(1:end-1);

rho_gs_10_1_arr = err_gs_10_1(2:end) ./ err_gs_10_1(1:end-1);
rho_gs_10_100_arr = err_gs_10_100(2:end) ./ err_gs_10_100(1:end-1);
rho_gs_10_001_arr = err_gs_10_001(2:end) ./ err_gs_10_001(1:end-1);

rho_gs_20_1_arr = err_gs_20_1(2:end) ./ err_gs_20_1(1:end-1);
rho_gs_20_100_arr = err_gs_20_100(2:end) ./ err_gs_20_100(1:end-1);
rho_gs_20_001_arr = err_gs_20_001(2:end) ./ err_gs_20_001(1:end-1);

rho_gs_40_1_arr = err_gs_40_1(2:end) ./ err_gs_40_1(1:end-1);
rho_gs_40_100_arr = err_gs_40_100(2:end) ./ err_gs_40_100(1:end-1);
rho_gs_40_001_arr = err_gs_40_001(2:end) ./ err_gs_40_001(1:end-1);

rho_sor_10_1_arr = err_sor_10_1(2:end) ./ err_sor_10_1(1:end-1);
rho_sor_10_100_arr = err_sor_10_100(2:end) ./ err_sor_10_100(1:end-1);
rho_sor_10_001_arr = err_sor_10_001(2:end) ./ err_sor_10_001(1:end-1);

rho_sor_20_1_arr = err_sor_20_1(2:end) ./ err_sor_20_1(1:end-1);
rho_sor_20_100_arr = err_sor_20_100(2:end) ./ err_sor_20_100(1:end-1);
rho_sor_20_001_arr = err_sor_20_001(2:end) ./ err_sor_20_001(1:end-1);

rho_sor_40_1_arr = err_sor_40_1(2:end) ./ err_sor_40_1(1:end-1);
rho_sor_40_100_arr = err_sor_40_100(2:end) ./ err_sor_40_100(1:end-1);
rho_sor_40_001_arr = err_sor_40_001(2:end) ./ err_sor_40_001(1:end-1);

rho_jac_10_1 = mean(quantile(rho_jac_10_1_arr(rho_jac_10_1_arr < 1), [0.25, 0.75]));
rho_jac_10_100 = mean(quantile(rho_jac_10_100_arr(rho_jac_10_100_arr < 1), [0.25, 0.75]));
rho_jac_10_001 = mean(quantile(rho_jac_10_001_arr(rho_jac_10_001_arr < 1), [0.25, 0.75]));

rho_jac_20_1 = mean(quantile(rho_jac_20_1_arr(rho_jac_20_1_arr < 1), [0.25, 0.75]));
rho_jac_20_100 = mean(quantile(rho_jac_20_100_arr(rho_jac_20_100_arr < 1), [0.25, 0.75]));
rho_jac_20_001 = mean(quantile(rho_jac_20_001_arr(rho_jac_20_001_arr < 1), [0.25, 0.75]));

rho_jac_40_1 = mean(quantile(rho_jac_40_1_arr(rho_jac_40_1_arr < 1), [0.25, 0.75]));
rho_jac_40_100 = mean(quantile(rho_jac_40_100_arr(rho_jac_40_100_arr < 1), [0.25, 0.75]));
rho_jac_40_001 = mean(quantile(rho_jac_40_001_arr(rho_jac_40_001_arr < 1), [0.25, 0.75]));

rho_gs_10_1 = mean(quantile(rho_gs_10_1_arr(rho_gs_10_1_arr < 1), [0.25, 0.75]));
rho_gs_10_100 = mean(quantile(rho_gs_10_100_arr(rho_gs_10_100_arr < 1), [0.25, 0.75]));
rho_gs_10_001 = mean(quantile(rho_gs_10_001_arr(rho_gs_10_001_arr < 1), [0.25, 0.75]));

rho_gs_20_1 = mean(quantile(rho_gs_20_1_arr(rho_gs_20_1_arr < 1), [0.25, 0.75]));
rho_gs_20_100 = mean(quantile(rho_gs_20_100_arr(rho_gs_20_100_arr < 1), [0.25, 0.75]));
rho_gs_20_001 = mean(quantile(rho_gs_20_001_arr(rho_gs_20_001_arr < 1), [0.25, 0.75]));

rho_gs_40_1 = mean(quantile(rho_gs_40_1_arr(rho_gs_40_1_arr < 1), [0.25, 0.75]));
rho_gs_40_100 = mean(quantile(rho_gs_40_100_arr(rho_gs_40_100_arr < 1), [0.25, 0.75]));
rho_gs_40_001 = mean(quantile(rho_gs_40_001_arr(rho_gs_40_001_arr < 1), [0.25, 0.75]));

rho_sor_10_1 = mean(quantile(rho_sor_10_1_arr(rho_sor_10_1_arr < 1), [0.25, 0.75]));
rho_sor_10_100 = mean(quantile(rho_sor_10_100_arr(rho_sor_10_100_arr < 1), [0.25, 0.75]));
rho_sor_10_001 = mean(quantile(rho_sor_10_001_arr(rho_sor_10_001_arr < 1), [0.25, 0.75]));

rho_sor_20_1 = mean(quantile(rho_sor_20_1_arr(rho_sor_20_1_arr < 1), [0.25, 0.75]));
rho_sor_20_100 = mean(quantile(rho_sor_20_100_arr(rho_sor_20_100_arr < 1), [0.25, 0.75]));
rho_sor_20_001 = mean(quantile(rho_sor_20_001_arr(rho_sor_20_001_arr < 1), [0.25, 0.75]));

rho_sor_40_1 = mean(quantile(rho_sor_40_1_arr(rho_sor_40_1_arr < 1), [0.25, 0.75]));
rho_sor_40_100 = mean(quantile(rho_sor_40_100_arr(rho_sor_40_100_arr < 1), [0.25, 0.75]));
rho_sor_40_001 = mean(quantile(rho_sor_40_001_arr(rho_sor_40_001_arr < 1), [0.25, 0.75]));

disp([rho_jac_10_1, rho_jac_10_100, rho_jac_10_001]);
disp([rho_jac_20_1, rho_jac_20_100, rho_jac_20_001]);
disp([rho_jac_40_1, rho_jac_40_100, rho_jac_40_001]);

disp([rho_gs_10_1, rho_gs_10_100, rho_gs_10_001]);
disp([rho_gs_20_1, rho_gs_20_100, rho_gs_20_001]);
disp([rho_gs_40_1, rho_gs_40_100, rho_gs_40_001]);

disp([rho_sor_10_1, rho_sor_10_100, rho_sor_10_001]);
disp([rho_sor_20_1, rho_sor_20_100, rho_sor_20_001]);
disp([rho_sor_40_1, rho_sor_40_100, rho_sor_40_001]);
