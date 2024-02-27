nodes = readmatrix("hw44.nod", 'FileType', "Text");
elements = readmatrix("hw44.ele", 'FileType', "Text");

point_source_node = 503;
point_source_pos = nodes(point_source_node, 2:3);

x_elem = zeros(3, length(elements));
y_elem = zeros(3, length(elements));

% Get positions of element nodes
for i=1:length(elements)
	for j=1:3
		x_elem(j, i) = nodes(elements(i, j+1), 2);
		y_elem(j, i) = nodes(elements(i, j+1), 3);
	end
end

% Find distance to point source for all points on a given element
elem_distance = sqrt(sum((x_elem - point_source_pos(1)).^2, 1) + sum((y_elem - point_source_pos(2)).^2, 1));

% Find the minimum to find element that contains point source
[~, min_elem] = min(elem_distance);
co = colororder();
% Plot to make sure it looks good
plot(x_elem([1 2 3 1], min_elem), y_elem([1 2 3 1], min_elem))
hold on;
plot(point_source_pos(1), point_source_pos(2), 'o', 'MarkerFaceColor', co(2,:));
title("Element " + num2str(min_elem) + " and Point Source");
text(x_elem(1, min_elem), y_elem(1, min_elem), " 1")
text(x_elem(2, min_elem), y_elem(2, min_elem) + .002, " 2")
text(x_elem(3, min_elem), y_elem(3, min_elem), " 3")
xlim([min(x_elem(:, min_elem)) - 0.01, max(x_elem(:, min_elem)) + 0.01])
ylim([min(y_elem(:, min_elem)) - 0.01, max(y_elem(:, min_elem)) + 0.01])
xlabel("x");
ylabel("y");
hold off;

% Find local (area) coordinates
A = det([ones(3,1) x_elem(:, min_elem) y_elem(:, min_elem)]) / 2;
sub_A = zeros(3,1);

for i=1:3
	other_nodes = 1:3;
	other_nodes = other_nodes(other_nodes ~= i);
	sub_A(i) = abs(det([ones(3,1) [point_source_pos; [x_elem(other_nodes, min_elem) y_elem(other_nodes, min_elem)]]]) / 2);
end

phi = sub_A / A;

% Export area coords
writematrix(phi, "point_source_local_coords");
