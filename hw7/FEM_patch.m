function fig = FEM_patch(elements, nodes, u, plot_title, cmap)
fig = figure;
if nargin < 4
	plot_title = "Contour Plot";
end

if nargin < 5
	cmap = colormap(parula(10));
end
patch("Faces", elements, "Vertices", [nodes, u], "FaceVertexCData", u, 'FaceColor', 'interp', 'EdgeAlpha', 0.1);
colorbar;
colormap(cmap);
title(plot_title);
xlabel("x");
ylabel("y");
axis equal
fig.Name = plot_title;
end
