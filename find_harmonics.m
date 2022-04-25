%% find_harmonics
function[theta, theta_root, theta_all] = find_harmonics(beta, beta_dot,...
    beta_star, theta, theta_dot, theta_star, trim_inflow, u_p, u_t, psi);
theta_old = 0;
psi_n = 1:361;

for n = 1:5
theta_knot = trim_inflow(5, :);
theta_nc = trim_inflow(8, :);
theta_ns = trim_inflow(6, :);
theta_sum(psi_n) = (theta_nc.*cos(n.*psi_n))...
    + (theta_ns.*sin(n.*psi_n))...
    + theta_old;
theta_root(psi_n) = theta_knot + theta_sum(psi_n);
theta_old = theta_sum(psi_n);
theta_all(n, psi_n) = theta_root(psi_n);
end 
theta = theta_root;
end