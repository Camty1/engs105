%% ENGS 105 HW 1
% Cameron Wolfe 1/16/24

%%  Setup
clear;

I0 = 1;
sigma = 1;
R = 1;
a = 0.1;
k = 3;

r_min = a;
r_max = R;

theta_min = 0;
theta_max = pi/2;

N_analytical = 50;

%% Analytical Solution
r_linspace = linspace(r_min, r_max, N_analytical);
theta_linspace = linspace(theta_min, theta_max, N_analytical);

[r_meshgrid, theta_meshgrid] = meshgrid(r_linspace, theta_linspace);

x = r_meshgrid .* cos(theta_meshgrid);
y = r_meshgrid .* sin(theta_meshgrid);

u_analytical = @(r,theta) -I0/sigma * R *((r.^k + a^(2*k)*r.^(-k))./(R^k + a^(2*k)*R^(-k))).*cos(k*theta);

grad_u_analytical_rhat = @(r,theta) I0/sigma*R*k*(r.^(k-1)-a.^(2*k)*r.^(-k-1))/(R^k+a^(2*k)*R^(-k)).*cos(k*theta);
grad_u_analytical_thetahat = @(r,theta) -I0/sigma*R*k*(r.^(k-1)-a^(2*k)*r.^(-k-1))/(R^k+a^(2*k)*R^(-k)).*sin(k*theta) ./ r;

grad_u_xhat = @(rhat, thetahat, theta) rhat .* cos(theta) - thetahat .* sin(theta);
grad_u_yhat = @(rhat, thetahat, theta) rhat .* sin(theta) + thetahat .* cos(theta);

u_meshgrid = u_analytical(r_meshgrid, theta_meshgrid);

grad_u_meshgrid_rhat = grad_u_analytical_rhat(r_meshgrid, theta_meshgrid);
grad_u_meshgrid_thetahat = grad_u_analytical_thetahat(r_meshgrid, theta_meshgrid);

grad_u_meshgrid_xhat = grad_u_xhat(grad_u_meshgrid_rhat, grad_u_meshgrid_thetahat, theta_meshgrid);
grad_u_meshgrid_yhat = grad_u_yhat(grad_u_meshgrid_rhat, grad_u_meshgrid_thetahat, theta_meshgrid);

%% Plotting Analytical Solution
colormap winter;
figure(1);
surf(x,y,u_meshgrid, 'EdgeAlpha',0.3);
view(15,30);
title(sprintf("Analytical Solution, N=%d", N_analytical));
xlabel("x");
ylabel('y');
zlabel("Potential");

figure(2);
quiver(x, y, grad_u_meshgrid_xhat, grad_u_meshgrid_yhat);
title(sprintf("Negative Gradient of Analytical Solution, N=%d", N_analytical));
xlabel("x");
ylabel("y");
xlim([-.06, 1.06]);
ylim([-.06, 1.06]);

%% Numerical Solution (coarse)
x_5 = readmatrix_fortran("output\x005.dat");
y_5 = readmatrix_fortran("output\y005.dat");
u_5 = readmatrix_fortran("output\u005.dat");
r_5 = readmatrix_fortran("output\r005.dat");
theta_5 = readmatrix_fortran("output\theta005.dat");

x_5_meshgrid = reshape(x_5, 5, 5);
y_5_meshgrid = reshape(y_5, 5, 5);
u_5_meshgrid = reshape(u_5, 5, 5);
r_5_meshgrid = reshape(r_5, 5, 5);
theta_5_meshgrid = reshape(theta_5, 5, 5);

figure(3);
surf(x_5_meshgrid, y_5_meshgrid, u_5_meshgrid);
view(15,30);
title("Numerical Solution, N=5");
zlabel("Potential");
xlabel("x");
ylabel("y");

error_5 = error_map(r_5_meshgrid, theta_5_meshgrid, u_5_meshgrid, u_analytical);

figure(4);
surf(x_5_meshgrid, y_5_meshgrid, error_5);
view(15, 30);
title("Error Map, N=5");
zlabel("Potential");
xlabel("x");
ylabel("y");

%% Error vs delta r
calculated_numbers = 5:5:100;
rms_errors = [];
delta_rs = [];

for n=calculated_numbers
    u = reshape(readmatrix_fortran(sprintf("output/u%03d.dat", n)), n, n);
    r = reshape(readmatrix_fortran(sprintf("output/r%03d.dat", n)), n, n);
    theta = reshape(readmatrix_fortran(sprintf("output/theta%03d.dat", n)), n, n);
    error = error_map(r, theta, u, u_analytical);
    delta_rs(end+1) = r(2) - r(1);
    rms_errors(end+1) = rms_error(error);
end

figure(5);
loglog(delta_rs, rms_errors);
title("RMS Error vs \Delta r");
xlabel("\Delta r");
ylabel("RMS Error");

%% Converged Numerical Solution
x_50 = readmatrix_fortran("output\x050.dat");
y_50 = readmatrix_fortran("output\y050.dat");
u_50 = readmatrix_fortran("output\u050.dat");
r_50 = readmatrix_fortran("output\r050.dat");
theta_50 = readmatrix_fortran("output\theta050.dat");

x_50_meshgrid = reshape(x_50, 50, 50);
y_50_meshgrid = reshape(y_50, 50, 50);
u_50_meshgrid = reshape(u_50, 50, 50);
r_50_meshgrid = reshape(r_50, 50, 50);
theta_50_meshgrid = reshape(theta_50, 50, 50);

colormap winter;
figure(6);
surf(x_50_meshgrid, y_50_meshgrid, u_50_meshgrid, 'EdgeAlpha', 0.3);
title("Converged Numerical Solution, N=50");
xlabel("x");
ylabel("y");
zlabel("Potential");
view(15,30);

% gradient calculation (current)
grad_r_50 = zeros(50);
grad_theta_50 = zeros(50);

for i=2:49
    for j=2:49
        grad_r_50(i, j) = -(u_50_meshgrid(i+1, j) - u_50_meshgrid(i-1, j)) / 2*(r_50(2) - r_50(1));
        grad_theta_50(i, j) = -1/r_50_meshgrid(i,j)*(u_50_meshgrid(i, j+1) - u_50_meshgrid(i, j-1)) / 2*(theta_50(51) - theta_50(1));
    end
end

grad_x_50 = grad_r_50 .* cos(theta_50_meshgrid) - grad_theta_50 .* sin(theta_50_meshgrid);
grad_y_50 = grad_r_50 .* sin(theta_50_meshgrid) + grad_theta_50 .* cos(theta_50_meshgrid);

figure(7);
quiver(x_50_meshgrid(2:49, 2:49), y_50_meshgrid(2:49, 2:49), grad_x_50(2:49, 2:49), grad_y_50(2:49, 2:49));
title("Negative Gradient of Numerical Solution, N=50");
xlabel("x");
ylabel("y");
