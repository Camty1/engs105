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
LHS_inv = LHS^-1;
LHS_inv_trans = LHS_inv';
RHS = readmatrix("output/RHS.dat");

noise = zeros(length(RHS), 3);
point_source = readmatrix("point_source_local_coords.txt");
correlation_lengths = [0, 2, 5];

cov_mat_error = zeros(502, 502, 5);
cov_mat_potential = zeros(502, 502, 6);

for i=1:3
	node = elements(288,i);
	noise(node, 1) = RHS(node) * 0.5;
end

cov_mat_error(:, :, 1) = cov_mat_error(:,:,1) + diag(noise(:,1).^2);

noise(492, 2) = abs(RHS(492) * 0.4);
noise(493, 2) = abs(RHS(493) * 0.4);

cov_mat_error(:, :, 2) = cov_mat_error(:, :, 2) + diag(noise(:,2).^2);

for i=1:length(bcs)
	for j=i:length(bcs)
		for k=1:length(correlation_lengths)
			cov_mat_error(bcs(i), bcs(j), 2 + k) = context_length_erm(i, j, correlation_lengths(k), 0.5, bcs, nodes);
			cov_mat_error(bcs(j), bcs(i), 2 + k) = cov_mat_error(bcs(i), bcs(j), 2 + k);
		end
	end
end

cov_mat_error(:, :, 6) = cov_mat_error(:, :, 1) + cov_mat_error(:, :, 2) + cov_mat_error(:, :, 4);

for i=1:6
	cov_mat_potential(:, :, i) = LHS_inv * cov_mat_error(:, :, i) * LHS_inv_trans;
end

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
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', sqrt(diag(cov_mat_potential(:, :, 1))), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Scenario 1");
xlabel("x");
ylabel("y");

figure(4);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', sqrt(diag(cov_mat_potential(:, :, 2))), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Scenario 2");
xlabel("x");
ylabel("y");

figure(5);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', sqrt(diag(cov_mat_potential(:, :, 3))), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Scenario 3, l = 0");
xlabel("x");
ylabel("y");

figure(6);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', sqrt(diag(cov_mat_potential(:, :, 4))), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Scenario 3, l = 2");
xlabel("x");
ylabel("y");

figure(7);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', sqrt(diag(cov_mat_potential(:, :, 5))), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Scenario 3, l = 5");
xlabel("x");
ylabel("y");

figure(8);
patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', sqrt(diag(cov_mat_potential(:, :, 6))), 'FaceColor', 'interp', 'EdgeAlpha', 0.1, 'DisplayName', 'Potential');
hold on;
plot(source_node(1), source_node(2), 'ko', 'MarkerFaceColor', 'red', 'DisplayName', 'Current Source');
plot(nodes(492:493, 1), nodes(492:493, 2), "Color", 'red', 'LineWidth', 1.5, 'DisplayName', 'Current Sink');
colorbar;
axis equal;
title("Map of Error Covariance Diagonal Combined");
xlabel("x");
ylabel("y");

function correlation = context_length_erm(i, j, l, sigma, bcs, nodes)
delta = (nodes(bcs(j), :) - nodes(bcs(i), :)).^2;
dist = sqrt(sum(delta));
correlation = sigma^2;

if (l ~= 0)
	correlation = sigma^2 * (1 + dist / l) * exp(-dist / l);
end
end
