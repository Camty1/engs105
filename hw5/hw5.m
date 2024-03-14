%% ENGS 105 HW 5
% Cameron Wolfe 3/12/2024

SAVE = false;
GIF = false;

%% Data loading
elem = readmatrix("epeltr4.dat");
elem_type = elem(:, 6);
elem = elem(:, 2:5);

node = readmatrix("npeltr4.dat");
node = node(:, 2:3);

u_lapack = readmatrix("output/u_lapack.dat");
u_ss_baseline = readmatrix("output/u_ss_baseline.dat");
u_ss_modified = readmatrix("output/u_ss_modified.dat");
u_trans_baseline = readmatrix("output/u_trans_baseline.dat")';
u_trans_modified = readmatrix("output/u_trans_modified.dat")';

tumor_elems = find(elem_type == 8);
bladder_elems = find(elem_type == 6);
other_elems = find(elem_type ~= 8);
tumor_nodes = unique(elem(tumor_elems, :));
bladder_nodes = unique(elem(bladder_elems, :));
tumor_center = mean(node(tumor_nodes, :), 1);
bladder_center = mean(node(bladder_nodes, :), 1);
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

%% Find closest, furthest, and middle nodes from tumor and bladder center
tumor_distances = sqrt(sum((node(tumor_nodes, :) - tumor_center).^2, 2));
[~, tumor_dist_idx] = sort(tumor_distances);
used_tumor_indices = tumor_nodes(tumor_dist_idx([1, floor(length(tumor_distances)/2), end]));

bladder_distances = sqrt(sum((node(bladder_nodes, :) - bladder_center).^2, 2));
[~, bladder_dist_idx] = sort(bladder_distances);
used_bladder_indices = bladder_nodes(bladder_dist_idx([1, floor(length(bladder_distances)/2), end]));

%% Transient comparison
timestep_names = {'quad', 'double', 'baseline', 'half', 'quarter',};
baseline_transients = cell(length(timestep_names), 1);
times = cell(size(baseline_transients));
timesteps = readmatrix("dts.dat");

for i=1:length(timestep_names)
	baseline_transients{i} = readmatrix(['output/u_trans_' timestep_names{i} '.dat'])';
	times{i} = 0:timesteps(i):50000;
end

interpolated_transients = zeros([size(baseline_transients{1}), 5]);
for i=1:length(timestep_names)
	interpolated_transients(:,:,i) = interp1(times{i}, baseline_transients{i}', times{1})';
end

% Calculate Relative Error
err = (interpolated_transients - interpolated_transients(:, :, 1));
err(:, 2:end, :) = err(:, 2:end, :) ./ interpolated_transients(:, 2:end, 1);

% Get average error
avg_err = reshape(mean(abs(err), [1,2]), [], 1);

%% Plotting
close all;
figs = {};
co = colororder;

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

for i=1:3
	figs{end+1} = figure;
	for j=1:5
		plot(times{j}, baseline_transients{j}(used_tumor_indices(i), :), 'DisplayName', ['Timestep = ', num2str(timesteps(j)), ' s']);
		if j == 1
			hold on;
		end
	end
	hold off;
	xlabel("Time (s)");
	ylabel(['Relative Temperature (' char(176) 'C)']);
	title(['Temperature of Node ' num2str(used_tumor_indices(i)) ' in Tumor for Different Timesteps']);
	legend('location', 'best');
end

for i=1:3
	figs{end+1} = figure;
	for j=1:5
		plot(times{j}, baseline_transients{j}(used_bladder_indices(i), :), 'DisplayName', ['Timestep = ', num2str(timesteps(j)), ' s']);
		if j == 1
			hold on;
		end
	end
	hold off;
	xlabel("Time (s)");
	ylabel(['Relative Temperature (' char(176) 'C)']);
	title(['Temperature of Node ' num2str(used_bladder_indices(i)) ' in Bladder for Different Timesteps']);
	legend('location', 'best');
end

if GIF
	animated_FEM_patch(elem, node, u_trans_baseline, "Baseline Temperature Transient");
	animated_FEM_patch(elem, node, u_trans_modified, "Increased Bloodflow Temperature Transient");
end

if SAVE
	for i=1:length(figs)
		saveas(figs{i}, regexprep(lower(['./figures/' figs{i}.Name '.png']), ' ', '_'));
	end
end
