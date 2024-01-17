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

