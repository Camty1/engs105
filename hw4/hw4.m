%% ENGS 105 HW4
% Cameron Wolfe 2/16/2024

%% File Inputs
node_positions = readmatrix("hw4.nod", "FileType", "text", "Delimiter", " ");
node_positions = node_positions(:,2:3);

incidence_list = readmatrix("hw4.ele1", "FileType", "text", "Delimiter", " ");
incidence_list = incidence_list(:, 2:4);

u_1 = readmatrix("output/u_1.dat", "FileType", "text", "Delimiter", ",");
u_2 = readmatrix("output/u_2.dat", "FileType", "text", "Delimiter", ",");
u_3 = readmatrix("output/u_3.dat", "FileType", "text", "Delimiter", ",");
v_1 = readmatrix("output/v_1.dat", "FileType", "text");
v_2 = readmatrix("output/v_2.dat", "FileType", "text");
v_3 = readmatrix("output/v_3.dat", "FileType", "text");

%% Plotting
close all;
fig_1 = Plot2dTriMesh(node_positions(:,1), node_positions(:,2), incidence_list, u_1, v_1);
title("Scenario 1");
saveas(fig_1, "figures/1.png");

fig_2 = Plot2dTriMesh(node_positions(:,1), node_positions(:,2), incidence_list, u_2, v_2);
title("Scenario 2");
saveas(fig_2, "figures/2.png");

fig_3 = Plot2dTriMesh(node_positions(:,1), node_positions(:,2), incidence_list, u_3, v_3);
title("Scenario 3");
saveas(fig_3, "figures/3.png");
