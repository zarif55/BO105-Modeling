%% find_coleman_inflow function
% r: radial station
% psi: azimuth angle
% v_inf: free stream velocity
% v_tip: tip velocity
% alpha: angle of attack
function[coleman_data] = find_coleman_inflow(r, psi, v_inf, v_tip, alpha)
%% constants
X(1) = 0;
epsilon_cole = 1;
lambda_i_cole = 0;
%% defining functions
kx_fn = @(X) tan(X/2);
lambda_io = @(lambda_i_TPP_UH60_FF) lambda_i_TPP_UH60_FF;
lambda_i_cole_fn = @(r_grid, psi_grid, lambda_io, kx) lambda_io.*(1+(kx.*r.*cos(psi))); %%% CHECK EQUATION
mu_z_fn = @(v_inf) v_inf.*sin(alpha)/v_tip;                         %%% LAMBDA_IO should be varied same as v_inf because same amount of indices
mu_x_fn = @(v_inf) v_inf.*cos(alpha)/v_tip;                         %%% BUT psi and lambda_io vary independently
X_fn = @(lambda_i_cole, mu_x, mu_z) atan(mu_x./(lambda_i_cole-mu_z));

%% Coleman Model
[r_grid, psi_grid] = meshgrid(r, psi);
while epsilon_cole > 0.004
    lambda_i_cole_old = lambda_i_cole;
    kx = kx_fn(X);
    lambda_i_cole = lambda_i_cole_fn(r_grid, psi_grid, lambda_io, kx);
    mu_z = mu_z_fn(v_inf);
    mu_x = mu_x_fn(v_inf);
    X = X_fn(lambda_i_cole, mu_x, mu_z);
    epsilon_cole = abs(lambda_i_cole - lambda_i_cole_old);
end



coleman_data = lambda_i_cole;
end