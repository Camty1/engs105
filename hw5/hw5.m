close all;

pos = readmatrix("npeltr4.dat");
u = readmatrix("output/u_type1.dat");
x = pos(:,2);
y = pos(:,3);

x_surf = reshape(x, [], 3);
y_surf = reshape(y, [], 3);
u_surf = reshape(u, [], 3);

element_list = readmatrix("epeltr4.dat");

element_list = element_list(:, 2:end-1);

tri_element_list = [];

for i=1:length(element_list)
	if (element_list(i, 3) ~= element_list(i, 4))
		tri_element_list = [tri_element_list; element_list(i, [1 2 3]); element_list(i, [1 3 4])];
	else
		tri_element_list = [tri_element_list; element_list(i, 1:3)];
	end
end

fig = Plot2dTriMesh(x, y, tri_element_list, u);
title("Type I Boundary Condition Case");
saveas(fig, "figures/type1.png");
