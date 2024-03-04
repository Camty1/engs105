SAVE = true;
NUM_NODE = 502;

u_truth = readmatrix("output/truth.dat");
u_fit = readmatrix("output/u_fit.dat")';
b_truth = readmatrix("output/RHS.dat")';

elements = readmatrix("hw44.ele", "FileType", "text");
elements = elements(:, 2:4);
nodes = readmatrix("hw44.nod", "FileType", "text");
nodes = nodes(1:NUM_NODE, 2:3);
sample_points = readmatrix("sample_points.dat");
dnd_nodes = readmatrix("hw44.dnd", "FileType", "text");
dnd_nodes = dnd_nodes(:, 1);
internal_source_nodes = 1:NUM_NODE;
internal_source_nodes(dnd_nodes) = [];

K_packed = readmatrix("output/LHS.dat");
K = packed2unpacked_coeff(K_packed);

b_fit = K * u_fit;
source = b_fit;
source(dnd_nodes) = 0;
source_error = source - b_truth;

S = readmatrix("output/S.dat")';
d = readmatrix("output/d.dat");

model_data_misfit = d - S * u_fit;

error = u_fit - u_truth;

rms_bias = sqrt(mean(error.^2))
rms_mdm = sqrt(mean(model_data_misfit.^2))
rms_source = sqrt(mean((b_fit(internal_source_nodes) - b_truth(internal_source_nodes)).^2))
rms_dnd = sqrt(mean((b_fit(dnd_nodes) - b_truth(dnd_nodes)).^2))

figs = {};
close all;
figs{end+1} = FEM_patch(elements, nodes, u_truth, "Truth Contour Plot");
hold on;
plot3(sample_points(:,1), sample_points(:,2), ones(NUM_MEAS, 1) * max(u_truth), 'r.');
hold off;
legend("Truth Potential", "Sample Points");

figs{end+1} = FEM_patch(elements, nodes, u_fit, "Fit Contour Plot");
hold on;
plot3(sample_points(:,1), sample_points(:,2), ones(NUM_MEAS, 1) * max(u_fit), 'r.');
hold off;
legend("Fit Potential", "Sample Points");

figs{end+1} = FEM_patch(elements, nodes, source, "Source Distribution Contour Plot");
hold on;
plot3(sample_points(:,1), sample_points(:,2), ones(NUM_MEAS, 1) * max(source), 'r.');
hold off;
legend("Source Distribution", "Sample Points");

figs{end+1} = FEM_patch(elements, nodes, error, "Error Map Contour Plot");
hold on;
plot3(sample_points(:,1), sample_points(:,2), ones(NUM_MEAS, 1) * max(error), 'r.');
hold off;
legend("Error", "Sample Points");

figs{end+1} = FEM_patch(elements, nodes, source_error, "Source Error Contour Plot");
hold on;
plot3(sample_points(:,1), sample_points(:,2), ones(NUM_MEAS, 1) * max(source_error), 'r.');
hold off;
legend("Source Error", "Sample Points");

if SAVE
	for i=1:length(figs)
		saveas(figs{i}, regexprep(lower(['figures/' figs{i}.Name '.png']), ' ', '_'));
	end
end
