D = 0.5;
sigma = 0.1;
dx = 0.1;
dt_a = 0.005;
dt_b = 0.03;

sim_params = readmatrix("sim_params.cfg", "FileType", "text");
L = sim_params(1);
x_0 = sim_params(2);
timesteps = sim_params(3);

x = readmatrix("output/x.dat");
t_a = 0:dt_a:((timesteps)*dt_a);
t_b = 0:dt_b:((timesteps)*dt_b);

[x_mesh, t_mesh_a] = meshgrid(x,t_a);
[x_mesh, t_mesh_b] = meshgrid(x,t_b);

A = readmatrix("output/A_5.0E-01_5.0E-01_1.0E-01_5.0E-03.dat");
B = readmatrix("output/B_5.0E-01_5.0E-01_1.0E-01_5.0E-03.dat");

figure(1)
surf(x_mesh, t_mesh_a, A, "EdgeAlpha", 0.3);
xlabel("x");
ylabel("t");
zlabel("u");
title("Scheme A u(x,t), r = 1/4");
ylim([0 0.5]);

figure(2);
surf(x_mesh, t_mesh_b, B, "EdgeAlpha", 0.3);
xlabel("x");
ylabel("t");
zlabel("u");
title("Scheme B u(x,t), r = 1.5");

figure(3);
contourf(x_mesh, t_mesh_a, A-B, 40, "EdgeAlpha", 0);
xlabel("x");
ylabel("t");
title("Difference Between Schemes A and B, r = 1/4");
colorbar;

figure(4);
tiledlayout(2,1);
nexttile;
plot(x, A(21:5:41, :));
title("Scheme A, t = 20\Delta t, 25\Delta t, ..., 50\Delta t");
xlabel("x");
ylabel("u");
xlim([4.8 5.2])
ylim([0.7 0.85])
nexttile;
plot(x, B(21:5:41, :));
title("Scheme B, t = 20\Delta t, 25\Delta t, ..., 40\Delta t");
xlim([4.8 5.2])
ylim([0.7 0.85])
xlabel("x");
ylabel("u");



