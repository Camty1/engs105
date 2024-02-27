potential = readmatrix("output/u.dat");
current = readmatrix("output/i.dat");

nodes = readmatrix("hw44_scaled.nod", 'FileType', 'Text');
elements = readmatrix("hw44.ele", 'FileType', 'Text');
element_centers = zeros(length(elements), 2);

source_node = nodes(503, 2:3);
nodes = nodes(1:502, 2:3);
elements = elements(:, 2:4);
u_elements = zeros(size(elements));

for l=1:length(elements)
	for i=1:3
		element_centers(l, 1) = element_centers(l, 1) + nodes(elements(l, i), 1) / 3;
		element_centers(l, 2) = element_centers(l, 2) + nodes(elements(l, i), 2) / 3;
	end
end
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


