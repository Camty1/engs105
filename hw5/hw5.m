close all;
set(0,'DefaultFigureWindowStyle','docked')

pos = readmatrix("npeltr4.dat");
u_1 = readmatrix("output/u_1.dat");
u_3 = readmatrix("output/u_3.dat");
u_1_tri = readmatrix("output/u_1_tri.dat");
u_3_tri = readmatrix("output/u_3_tri.dat");
u_1_tri_dt = readmatrix("output/u_1_tri_dt.dat");
LHS_3 = readmatrix("output/LHS_3.dat");
LHS_3_tri = readmatrix("output/LHS_3_tri.dat");
RHS_3 = readmatrix("output/RHS_3.dat");
RHS_3_tri = readmatrix("output/RHS_3_tri.dat");

x = pos(:,2);
y = pos(:,3);

save = false;

element_list = readmatrix("epeltr4.dat");
element_list = element_list(:, 2:end-1);

tri_element_list = readmatrix("epeltr4_tri.dat");
tri_element_list = tri_element_list(:, 2:end-1);

fig = Plot2dTriMesh(x, y, tri_element_list, u_1);
title("Type I Boundary Condition");
if save
    saveas(fig, "figures/1.png");
end
fig = Plot2dTriMesh(x, y, tri_element_list, u_1_tri);
title("Type I Boundary Condition, Triangular Elements");
if save
saveas(fig, "figures/1_tri.png");
end

fig = Plot2dTriMesh(x, y, tri_element_list, u_3);
title("Type III Boundary Condition");
if save
saveas(fig, "figures/3.png");
end

fig = Plot2dTriMesh(x, y, tri_element_list, u_3_tri);
title("Type III Boundary Condition, Triangular Elements");
if save
saveas(fig, "figures/3_tri.png");
end
fig = Plot2dTriMesh(x, y, tri_element_list, u_1_tri_dt(end, :)');
title("Type I Boundary Condition, Triangular Elements, k = 500");
