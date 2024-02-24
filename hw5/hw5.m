close all;
set(0,'DefaultFigureWindowStyle','docked')

pos = readmatrix("npeltr4.dat");
u_1 = readmatrix("output/u_1.dat");
u_3 = readmatrix("output/u_3.dat");
u_1_tri = readmatrix("output/u_1_tri.dat");
u_3_tri = readmatrix("output/u_3_tri.dat");
x = pos(:,2);
y = pos(:,3);

x_surf = reshape(x, [], 3);
y_surf = reshape(y, [], 3);
u_surf = reshape(u, [], 3);

element_list = readmatrix("epeltr4.dat");
element_list = element_list(:, 2:end-1);

tri_element_list = readmatrix("epeltr4_tri.dat");
tri_element_list = tri_element_list(:, 2:end-1);

fig = Plot2dTriMesh(x, y, tri_element_list, u_1);
title("Type I Boundary Condition");
saveas(fig, "figures/1.png");

fig = Plot2dTriMesh(x, y, tri_element_list, u_1_tri);
title("Type I Boundary Condition, Triangular Elements");
saveas(fig, "figures/1_tri.png");

fig = Plot2dTriMesh(x, y, tri_element_list, u_3);
title("Type III Boundary Condition");
saveas(fig, "figures/3.png");

fig = Plot2dTriMesh(x, y, tri_element_list, u_3_tri);
title("Type III Boundary Condition, Triangular Elements");
saveas(fig, "figures/3_tri.png");
