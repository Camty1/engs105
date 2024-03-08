SAVE = true;

elem = readmatrix("epeltr4.dat");
elem_type = elem(:, 6);
elem = elem(:, 2:5);

node = readmatrix("npeltr4.dat");
node = node(:, 2:3);

LHS_lapack = readmatrix("output/LHS_lapack.dat");
LHS_unpacked = readmatrix("output/LHS_unpacked.dat");
RHS = readmatrix("output/RHS.dat");

u_lapack = readmatrix("output/u_lapack.dat");
u_ss_baseline = readmatrix("output/u_ss_baseline.dat");
u_ss_modified = readmatrix("output/u_ss_modified.dat");
u_transient = readmatrix("output/u_transient.dat")';

tumor_elems = find(elem_type == 8);
other_elems = find(elem_type ~= 8);
tumor_nodes = unique(elem(tumor_elems, :));
tumor_center = mean(node(tumor_nodes, :), 1);
other_nodes = unique(elem(other_elems, :));
boundary_nodes = intersect(tumor_nodes, other_nodes);
boundary_node_angles = zeros(length(boundary_nodes), 1);

for i=1:length(boundary_nodes)
	boundary_node_angles(i) = atan2(node(boundary_nodes(i), 2) - tumor_center(2), node(boundary_nodes(i), 1) - tumor_center(1));
end

[~, boundary_node_order] = sort(boundary_node_angles);
boundary_nodes = boundary_nodes(boundary_node_order);
boundary_nodes(end+1) = boundary_nodes(1);

baseline_tumor_temp = mean(u_ss_baseline(tumor_nodes));
modified_tumor_temp = mean(u_ss_modified(tumor_nodes));

co = colororder;

close all;
figs = {};

figs{end+1} = FEM_patch(elem, node, u_ss_baseline, "Baseline Temperature Contour Plot");
clim([min([u_ss_baseline; u_ss_modified]), max([u_ss_baseline; u_ss_modified])]);
hold on;
plot3(node(boundary_nodes, 1), node(boundary_nodes, 2), ones(size(boundary_nodes, 1), 1) * max(u_ss_baseline), "Color", co(2, :), "LineWidth", 2);
legend("Temperature", ['Tumor, T_{AVG} = ' num2str(baseline_tumor_temp) char(176) 'C']);
hold off;

figs{end+1} = FEM_patch(elem, node, u_ss_modified, "Increased Bloodflow Temperature Contour Plot");
clim([min([u_ss_baseline; u_ss_modified]), max([u_ss_baseline; u_ss_modified])]);
hold on;
plot3(node(boundary_nodes, 1), node(boundary_nodes, 2), ones(size(boundary_nodes, 1), 1) * max(u_ss_modified), "Color", co(2, :), "LineWidth", 2);
legend("Temperature", ['Tumor, T_{AVG} = ' num2str(modified_tumor_temp) char(176) 'C']);
hold off;

[bingo, bongo] = animated_FEM_patch(elem, node, u_transient, "Baseline Temperature Transient");

if SAVE
	for i=1:length(figs)
		saveas(figs{i}, regexprep(lower(['./figures/' figs{i}.Name '.png']), ' ', '_'));
	end
end
