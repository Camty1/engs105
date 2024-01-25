x = readmatrix_fortran('output/x080.dat');
y = readmatrix_fortran('output/y080.dat');

u_ss = readmatrix_fortran('output/u_ss_1_080.dat');
u_transient_0_025 = readmatrix('output/u_transient_1_0.0E+00_2.5E-01.dat');

x_mesh = reshape(x, 80, 80);
y_mesh = reshape(y, 80, 80);
u_mesh = reshape(u_ss, 80, 80);
u_transient_mesh = reshape(u_transient_0_025, 80, 80, []);

figure(1)
contourf(x_mesh, y_mesh, u_mesh, 20, 'EdgeAlpha', 0.3);

figure(2)
surf(x_mesh, y_mesh, u_transient_mesh(:,:,1), 'EdgeAlpha', 0.3);


