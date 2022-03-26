%% find_bemt function
% r: radial station
% psi: azimuth angle
% v_inf: free stream velocity
% v_tip: tip velocity
% alpha: angle of attack

function[u_t, u_p_linear, u_p_MS, T_linear_bar, Q_linear_bar, T_MS_bar,...
    Q_MS_bar, P_linear_bar, P_MS_bar, dL_MS, alpha_linear,alpha_MS, dFx_MS, dFz_MS] = ...
    find_bemt(r, psi, mu, beta, beta_star,v_inf, lambda_TPP_UH60_FF,...
    Cl_alpha,lambda_MS_total, theta, rho, c, R, cd0, N_psi, dr, Nb, omega, v_tip)
% mu is a scalar
% lambda_TPP_UH_60_FF is a scalar
% want lambda_1_i_MS to be length(r) x length(psi)
%% constants
xpr = r.*R;
%% defining functions
u_t_fn = @(r,psi) v_tip*(r + mu.*sin(psi));
u_p_linear_fn = @(r, psi, beta, beta_star) v_tip*(lambda_TPP_UH60_FF + ...
    r.*beta_star + mu.*beta.*cos(psi));
u_p_MS_fn = @(r, psi, beta, beta_star) v_tip*(lambda_MS_total + ...
    r.*beta_star + mu.*beta.*cos(psi));
phi_linear_fn = @(u_t_linear, u_p_linear) u_p_linear./u_t_linear;
phi_MS_fn = @(u_p_MS, u_t_MS) u_p_MS./u_t_MS;
alpha_linear_fn = @(u_t_linear, u_p_linear) theta - u_p_linear./u_t_linear;
alpha_MS_fn = @(u_p_MS, u_t_MS) theta - u_p_MS./u_t_MS;
Cl_linear_fn = @(alpha_linear) Cl_alpha.*alpha_linear;
Cl_MS_fn = @(alpha_MS) Cl_alpha.*alpha_MS;
dL_linear_fn = @(u_t, Cl_linear) 0.5.*rho.*(u_t.^2).*c.*dr.*Cl_linear;
dL_MS_fn = @(u_t, Cl_MS) 0.5.*rho.*(u_t.^2).*c.*dr.*Cl_MS;
dD_fn = @(u_t) 0.5.*rho.*(u_t.^2).*c.*dr.*cd0;        % is it cd0 NOT SURE
dFz_linear_fn = @(dL_linear, phi_linear, dD) dL_linear.*cos(phi_linear) - dD.*sin(phi_linear);
dFz_MS_fn = @(dL_MS, phi_MS, dD) dL_MS.*cos(phi_MS) - dD.*sin(phi_MS);
dFx_linear_fn = @(dL_linear, phi_linear, dD) dL_linear.*sin(phi_linear) + dD.*cos(phi_linear);
dFx_MS_fn = @(dL_MS, phi_MS, dD) dL_MS.*sin(phi_MS) + dD.*cos(phi_MS);

T_linear_psi_fn = @(psi, dFz_linear) sum(dFz_linear, 2);
T_linear_bar_fn = @(T_linear_psi) (Nb/N_psi)*sum(T_linear_psi);
T_MS_psi_fn = @(psi, dFz_MS) sum(dFz_MS);
T_MS_bar_fn = @(T_MS_psi) (Nb/N_psi)*sum(T_MS_psi);


Q_linear_psi_fn = @(r, psi, dFx_linear) sum(R*r.*dFx_linear, 1);
Q_linear_bar_fn = @(Q_linear_psi) (Nb/N_psi)*sum(Q_linear_psi);
Q_MS_psi_fn = @(r, psi, dFx_MS_fn) (sum(R*r.*dFx_MS_fn));
Q_MS_bar_fn = @(Q_MS_psi) (Nb/N_psi)*sum(Q_MS_psi);

P_linear_bar_fn = @(Q_linear_bar) Q_linear_bar*omega;
P_MS_bar_fn = @(Q_MS_bar) abs(Q_MS_bar)*omega;

%% 
[r_bemt_grid, beta_grid] = meshgrid(r, beta);
[r_bemt_grid, beta_star_grid] = meshgrid(r, beta_star);
[r_bemt_grid, psi_bemt_grid] = meshgrid(r, psi);

u_t = u_t_fn(r_bemt_grid,psi_bemt_grid);
u_p_linear = u_p_linear_fn(r_bemt_grid,psi_bemt_grid,beta_grid, beta_star_grid);
u_p_MS = u_p_MS_fn(r_bemt_grid,psi_bemt_grid,beta_grid, beta_star_grid);
phi_linear = phi_linear_fn(u_t, u_p_linear);
phi_MS = phi_MS_fn(u_p_MS, u_t);
alpha_linear = alpha_linear_fn(u_t, u_p_linear);
alpha_MS = alpha_MS_fn(u_p_MS, u_t);
Cl_linear = Cl_linear_fn(alpha_linear);
Cl_MS = Cl_MS_fn(alpha_MS);

dL_linear = dL_linear_fn(u_t, Cl_linear); 
dL_MS = dL_MS_fn(u_t, Cl_MS);
dD = dD_fn(u_t);
dFz_linear = dFz_linear_fn(dL_linear, phi_linear,dD);
dFz_MS = dFz_MS_fn(dL_MS, phi_MS, dD);
dFx_linear = dFx_linear_fn(dL_linear, phi_linear, dD); 
dFx_MS = dFx_MS_fn(dL_MS, phi_MS, dD); 
T_linear_psi = T_linear_psi_fn(psi, dFz_linear);
T_linear_bar = T_linear_bar_fn(T_linear_psi);
T_MS_psi = T_linear_psi_fn(psi, dFz_linear);
T_MS_bar = T_linear_bar_fn(T_MS_psi);

Q_linear_psi = Q_linear_psi_fn(r, psi, dFx_linear);
Q_linear_bar = Q_linear_bar_fn(Q_linear_psi);
Q_MS_psi = Q_MS_psi_fn(r, psi, dFx_MS);
Q_MS_bar = Q_MS_bar_fn(Q_MS_psi);

P_linear_bar = P_linear_bar_fn(Q_linear_bar);
P_MS_bar = P_MS_bar_fn(Q_MS_bar);





