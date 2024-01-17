function error = error_map(r_mesh, theta_mesh, u_mesh, u_func)
    u_analytical = u_func(r_mesh, theta_mesh);

    error = u_mesh - u_analytical;
end