theta = [0 0 0.3 0.3 0.6];
r = [0.25 0.60 0.75 1.5 10.0];
N = 80;
D = 1;

r_min = 0.1;
r_max = 1.0;
phi_min = 0;
phi_max = pi / 2;

dr = (r_max - r_min) / (N - 1);
dphi = (phi_max - phi_min) / (N - 1);

set(0,'DefaultFigureWindowStyle','docked')
tracked_indices = randperm(6400, 5)

x = readmatrix_fortran('output/x080.dat');
y = readmatrix_fortran('output/y080.dat');
u_ss_1 = readmatrix('output/u_ss_1_080.dat');
u_ss_3 = readmatrix('output/u_ss_3_080.dat');

x_mesh = reshape(x, 80, 80);
y_mesh = reshape(y, 80, 80);
u_ss_1_mesh = reshape(u_ss_1, 80, 80);
u_ss_3_mesh = reshape(u_ss_3, 80, 80);

figure(1);
surf(x_mesh, y_mesh, u_ss_1_mesh, 'EdgeAlpha', 0.3);
title("Steady State for Type I BC");
xlabel("x");
ylabel("y");
zlabel("Potential")
view(15,30);
colorbar;
saveas(gcf, "figures/type_1_ss.png");

figure(2);
surf(x_mesh, y_mesh, u_ss_3_mesh, 'EdgeAlpha', 0.3);
title("Steady State for Type III BC");
xlabel("x");
ylabel("y");
zlabel("Potential")
view(15,30);
saveas(gcf, "figures/type_3_ss.png");

colorbar;

for i=1:5
	for type=1:2
		dt = r(i) * dr^2 / D;
		if type == 1
			u_transient = readmatrix(sprintf('output/u_transient_%d_%03d_%7.1E_%7.1E_%d.dat', 1, N, theta(i), r(i), D))';
		else
			u_transient = readmatrix(sprintf('output/u_transient_%d_%03d_%7.1E_%7.1E_%d.dat', 3, N, theta(i), r(i), D))';
		end
		u_tracked = u_transient(tracked_indices, :);
		
		figure(4*i+2*(type-1)-1);
		plot((0:1000)*dt, u_tracked);
		if type == 1
			title(sprintf("Type I Transient Response, theta = %.2f, r = %.2f", theta(i), r(i)));
		else
			title(sprintf("Type III Transient Response, theta = %.2f, r = %.2f", theta(i), r(i)));
		end
		xlabel("Time");
		ylabel("Potential");
		saveas(gcf, sprintf("figures/trans_%d_%d.png", i, type));
		
		figure(4*i+2*(type-1));
		if type == 1
			contourf(x_mesh, y_mesh, reshape(u_transient(:,end) - u_ss_1, 80, 80), 20, "EdgeAlpha", 0.3);
			title(sprintf("Type I Steady State Error Map, theta = %.2f, r = %.2f", theta(i), r(i)));
			colorbar;
			xlabel("x");
			ylabel("y");
			saveas(gcf, sprintf("figures/error_%d_%d.png", i, type));
		else
			contourf(x_mesh, y_mesh, reshape(u_transient(:,end) - u_ss_3, 80, 80), 20, "EdgeAlpha", 0.3);
			title(sprintf("Type III Steady State Error Map, theta = %.2f, r = %.2f", theta(i), r(i)));
			colorbar;
			xlabel("x");
			ylabel("y");
			saveas(gcf, sprintf("figures/error_%d_%d.png", i, type));
		end
		
	end
end

u_ss_1 = readmatrix(sprintf('output/u_ss_1_%03d.dat', N));
u_transient_1 = readmatrix(sprintf('output/u_transient_%d_%03d_%7.1E_%7.1E_%d.dat', 1, N, theta(5), r(5), D))';

u_ss_3 = readmatrix(sprintf('output/u_ss_3_%03d.dat', N));
u_transient_3 = readmatrix(sprintf('output/u_transient_%d_%03d_%7.1E_%7.1E_%d.dat', 3, N, theta(5), r(5), D))';

error_1 = zeros(size(u_transient_1));
error_3 = zeros(size(u_transient_3));
for i=1:size(u_transient_1, 2)
	error_1(:,i) = abs(u_transient_1(:,i) - u_ss_1);
	error_3(:,i) = abs(u_transient_3(:,i) - u_ss_3);
end

dt = r(5) * dr^2 / D;
error_1 = error_1(mod(1:N^2, N) ~= 0, 1:200);

ln_error_1 = log(error_1);
delta_ln_error_1 = error_1(:,2:end) - error_1(:,1:end-1);
tau_1 = abs(dt ./ ln_error_1);
tau_1 = mean(tau_1(tau_1 > 1e-4), "all");

error_3 = error_3(mod(1:N^2, N) ~= 0, 1:200);

ln_error_3 = log(error_3);
delta_ln_error_3 = error_3(:,2:end) - error_3(:,1:end-1);
tau_3 = abs(dt ./ ln_error_3);
tau_3 = mean(tau_3(tau_1 > 1e-4), "all");
