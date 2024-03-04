SAVE = true;
NUM_MEAS = 12;

u_mat = readmatrix("output/u_mat.dat");
R_plus_I = readmatrix("output/R_plus_I.dat");

cov_d = eye(NUM_MEAS) * 0.05^2;

R_inv = (R_plus_I)^-1;

cov_u = u_mat * R_inv * cov_d * R_inv' * u_mat';

elements = readmatrix("hw44.ele", "FileType", "text");
elements = elements(:, 2:4);
nodes = readmatrix("hw44.nod", "FileType", "text");
nodes = nodes(1:502, 2:3);

imprecision = sqrt(diag(cov_u));

sample_points = readmatrix("sample_points.dat");

close all;
fig = FEM_patch(elements, nodes, imprecision, "Imprecision Contour Plot");
hold on;
plot3(sample_points(:,1), sample_points(:,2), ones(NUM_MEAS, 1) * max(imprecision), 'r.');
hold off;
legend("Imprecision Contour Map", "Sample Points");

if SAVE
	saveas(fig, regexprep(lower(['figures/' fig.Name '.png']), ' ', '_'));
end
