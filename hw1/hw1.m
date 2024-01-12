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
N_coarse = 5;
N_fine = 50;

r_linspace = linspace(r_min, r_max, N_analytical);
theta_linspace = linspace(theta_min, theta_max, N_analytical);

[r_meshgrid, theta_meshgrid] = meshgrid(r_linspace, theta_linspace);

x = r_meshgrid .* cos(theta_meshgrid);
y = r_meshgrid .* sin(theta_meshgrid);

u_analytical = @(r,theta) -I0/sigma * R *((r.^k + a^(2*k)*r.^(-k))./(R^k + a^(2*k)*R^(-k))).*cos(k*theta);

grad_u_analytical_rhat = @(r,theta) I0/sigma*R*k*(r.^(k-1)-a.^(2*k)*r.^(-k-1))/(R^k+a^(2*k)*R^(-k)).*cos(k*theta);
grad_u_analytical_thetahat = @(r,theta) -I0/sigma*R*k*(r^(k-1)-a^(2*k)*r.^(-k-1))/(R^k+a^(2*k)*R^(-k)).*sin(k*theta);

grad_u_xhat = @(rhat, thetahat, theta) rhat .* cos(theta) - thetahat .* sin(theta);
grad_u_yhat = @(rhat, thetahat, theta) rhat .* sin(theta) + thetahat .* cos(theta);

u_meshgrid = u_analytical(r_meshgrid, theta_meshgrid);

grad_u_meshgrid_rhat = grad_u_analytical_rhat(r_meshgrid, theta_meshgrid);
grad_u_meshgrid_thetahat = grad_u_analytical_thetahat(r_meshgrid, theta_meshgrid);

grad_u_meshgrid_xhat = grad_u_xhat(grad_u_meshgrid_rhat, grad_u_meshgrid_thetahat, theta_meshgrid);
grad_u_meshgrid_yhat = grad_u_yhat(grad_u_meshgrid_rhat, grad_u_meshgrid_thetahat, theta_meshgrid);

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