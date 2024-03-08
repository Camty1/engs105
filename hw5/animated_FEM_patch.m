function [fig, mov] = FEM_patch(elements, nodes, u, plot_title, cmap)
fig = figure;
if nargin < 4
	plot_title = "Contour Plot";
end

if nargin < 5
	cmap = colormap(parula(10));
end

T = (0:size(u, 2)-1) * 50;
p = patch("Faces", elements, "Vertices", [nodes, u(:, 1)], "FaceVertexCData", u(:, 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.1);
pos = get(fig, 'Position');
colorbar;
colormap(cmap);
title(plot_title);
xlabel("x");
ylabel("y");
axis equal
fig.Name = plot_title;
x_limits = get(gca, "XLim");
y_limits = get(gca, "YLim");
z_limits = [-16.909488231494944, 14.537130145485376];
xlim(x_limits);
ylim(y_limits);
zlim(z_limits);
clim(z_limits);

width = pos(3);
height = pos(4);
mov = zeros(height, width, 1, length(T), 'uint8');

f = getframe(fig);
[mov(:,:,1,1), map] = rgb2ind(f.cdata, 256, 'nodither');

for i=2:length(T)
    set(p, "Vertices", [nodes, u(:, i)], "FaceVertexCData", u(:, i));
    f = getframe(fig);
    mov(:, :, 1, i) = rgb2ind(f.cdata, map, 'nodither');
end

imwrite(mov, map, './figures/u_transient.gif', 'DelayTime', 0, 'LoopCount', inf);

end
