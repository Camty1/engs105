num_measurements = 12;
num_nodes = 502;
num_grounds = 33;

nodes = readmatrix("hw44.nod", "FileType", "text");
nodes = nodes(:, 2:3);

ground = readmatrix("hw44.dnd", "FileType", "text");
ground = ground(:,1);

data_model_misfit = 0.05;

l_source = 0.3;
l_ground = 0.4;

sigma_source = 1.0;
sigma_ground = 0.1;

cov_delta = eye(num_measurements) * data_model_misfit^2;
cov_ground = zeros(num_nodes);
cov_source = zeros(num_nodes);

for i=1:num_nodes
	for j=i:num_nodes
		
		cov_source(i,j) = context_length_erm(i, j, l_source, sigma_source, nodes);
		cov_source(j,i) = cov_source(i,j);
		
	end
end

for i=1:num_grounds
	for j=i:num_grounds
		cov_ground(ground(i), ground(j)) = context_length_erm(i, j, l_ground, sigma_ground, nodes, ground);
		cov_ground(ground(j), ground(i)) = cov_ground(ground(i), ground(j));
	end
	
	cov_source(ground(i), :) = 0;
	cov_source(:, ground(i)) = 0;
end

cov_b = cov_ground + cov_source;

writematrix(cov_b, "output/cov_b.dat");
writematrix(cov_delta, "output/cov_delta.dat");
writematrix(cov_delta^-1, "output/cov_delta_inv.dat");

function correlation = context_length_erm(i, j, l, sigma, nodes, bcs)

if nargin == 5
	delta = (nodes(j, :) - nodes(i, :)).^2;
else
	delta = (nodes(bcs(j), :) - nodes(bcs(i), :)).^2;
end

dist = sqrt(sum(delta));
correlation = sigma^2;

if (l ~= 0)
	correlation = sigma^2 * (1 + dist / l) * exp(-dist / l);
end

end
