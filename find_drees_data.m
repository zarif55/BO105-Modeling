%% find_drees_data function
% r: radial station
% psi: azimuth angle
% v_inf: free stream velocity
% v_tip: tip velocity
% alpha: angle of attack
function[drees_data] = find_drees_data(r, psi, mu, v_inf, v_tip, alpha)
%% constants
lambda_io = lambda_TPP_UH60_FF;
%% Drees Model   
kx_dree_fn = @(mu, ) (4/3)*((1-cos(X)-(1.8*mu.^2))/(sin(psi)));
ky_dree_fn = @(mu) -2*mu;
lambda_i_dree = @(lambda_io, kx_dree) lambda_io*(1+(kx_dree.*r.*cos(psi)));


for m = 1:223
for n = 1:100
    kx_dree = @()
    ky_dree = -2*(mu(m));
    lambda_i_dree = lambda_io(n)*(1+(kx*r*cos(psi(n))));
end
end

