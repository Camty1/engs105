%% ENGS 105 HW 6
% Cameron Wolfe 3/1/2024

%% File inputs
potential = readmatrix("output/u.dat");
potential_lapack = readmatrix("output/u_lapack.dat");
current = readmatrix("output/i.dat");
current_lapack = readmatrix("output/i_lapack.dat");

nodes = readmatrix("hw44_scaled.nod", 'FileType', 'Text');
elements = readmatrix("hw44.ele", 'FileType', 'Text');
bcs = readmatrix("hw44.dnd", 'FileType', 'Text');
element_centers = zeros(length(elements), 2);

source_node = nodes(503, 2:3);
nodes = nodes(1:502, 2:3);
elements = elements(:, 2:4);
bcs = bcs(:,1);

%% Covariance Calcs
LHS_packed = readmatrix("output/LHS.dat");
LHS = packed2unpacked_coeff(LHS_packed);
RHS = readmatrix("output/RHS.dat");

noise = zeros(length(RHS), 3);
point_source = readmatrix("point_source_local_coords.txt");

cov_mat_error = zeros(502, 502, 3);
cov_mat_potential = zeros(502, 502, 4);

for i=1:3
	node = elements(288,i);
	noise(node, 1) = RHS(i) * 0.5;
end

noise(492, 2) = abs(RHS(492) * 0.4);
noise(493, 2) = abs(RHS(493) * 0.4);

for i=1:length(bcs)
	bc = bcs(i);
	noise(bc, 3) = 0.5;
end

inverse_noise = LHS \ noise;
inverse_noise(:, 4) = sum(inverse_noise,2);
cov_diag = sqrt(inverse_noise .^ 2);

for l=1:length(elements)
	for i=1:3
		element_centers(l, 1) = element_centers(l, 1) + nodes(elements(l, i), 1) / 3;
		element_centers(l, 2) = element_centers(l, 2) + nodes(elements(l, i), 2) / 3;
	end
end

%% Plotting
close all;
figure(1);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', potential, 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
quiver(element_centers(:, 1), element_centers(:, 2), current(:, 1), current(:, 2), 'k', 'LineWidth', 1.5, 'DisplayName', 'Current');
hold off;
colorbar;
title("Potential and Current in Plate");
xlabel("x");
ylabel("y");
legend("Location", "northwest");
axis equal;

figure(2);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', potential_lapack, 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
quiver(element_centers(:, 1), element_centers(:, 2), current_lapack(:, 1), current_lapack(:, 2), 'k', 'LineWidth', 1.5, 'DisplayName', 'Current');
hold off;
colorbar;
title("Potential and Current in Plate");
xlabel("x");
ylabel("y");
legend("Location", "northwest");
axis equal;

figure(3);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', cov_diag(:,1), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
clim([0 max(cov_diag(:, 4))]);
axis equal;
title("Map of Error Covariance Diagonal Scenario 1");
xlabel("x");
ylabel("y");

figure(4);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', cov_diag(:,2), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Scenario 2");
xlabel("x");
ylabel("y");

figure(5);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', cov_diag(:,3), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
clim([0 max(cov_diag(:, 4))]);
axis equal;
title("Map of Error Covariance Diagonal Scenario 3");
xlabel("x");
ylabel("y");

figure(6);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', cov_diag(:,4), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
clim([0 max(cov_diag(:, 4))]);
axis equal;
title("Map of Error Covariance Diagonal Combined");
xlabel("x");
ylabel("y");
