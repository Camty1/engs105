function covariance_matrices(sigma_model_data_misfit, sigma_source, sigma_ground, l_source, l_ground)

%% Manage inputs
if nargin == 0
	sigma_model_data_misfit = 0.05;
	sigma_source = 1.0;
	sigma_ground = 0.1;
	l_source = 0.3;
	l_ground = 0.4;
elseif nargin == 1
	sigma_source = 1.0;
	sigma_ground = 0.1;
	l_source = 0.3;
	l_ground = 0.4;
elseif nargin == 3
	l_source = 0.3;
	l_ground = 0.4;
end

%% Parameters
NUM_MEAS = 12;
NUM_NODE = 502;
NUM_GRND = 33;

%% File inputs
nodes = readmatrix("hw44.nod", "FileType", "text");
nodes = nodes(:, 2:3);

ground = readmatrix("hw44.dnd", "FileType", "text");
ground = ground(:,1);

%% Model data misfit covariance matrix
cov_delta = eye(NUM_MEAS) * sigma_model_data_misfit^2;

%% Source covariance matrix
cov_source = zeros(NUM_NODE);

for i=1:NUM_NODE
	for j=i:NUM_NODE
		cov_source(i,j) = context_length_erm(i, j, l_source, sigma_source, nodes);
		cov_source(j,i) = cov_source(i,j);
	end
end

for i=1:NUM_GRND
	cov_source(ground(i), :) = 0;
	cov_source(:, ground(i)) = 0;
end

%% Dirichlet BC (ground) covariance matrix
cov_ground = zeros(NUM_NODE);

for i=1:NUM_GRND
	for j=i:NUM_GRND
		cov_ground(ground(i), ground(j)) = context_length_erm(i, j, l_ground, sigma_ground, nodes, ground);
		cov_ground(ground(j), ground(i)) = cov_ground(ground(i), ground(j));
	end
end

%% RHS (b) covariance matrix
cov_b = cov_ground + cov_source;

%% File outputs
writematrix(cov_b, "output/cov_b.dat");
writematrix(cov_delta, "output/cov_delta.dat");
writematrix(cov_delta^-1, "output/cov_delta_inv.dat");

end
