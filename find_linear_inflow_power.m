%% find_linear_inflow_power function
% sigma: blade solidity
% f:
% A: rotor disk area
% L_tr: moment arm tail rotor
% R: total radius
function[linear_inflow_power] = find_linear_inflow_power(mu, CT_constant_T, sigma, f, A,k)
%% constants
lambda_climb_FF = 0;
cd0 = 0.001;
%% defining functions
% c(r) is constant because assume rectangular blade
% blade elements are all assumed as the same with the cd0 being for airfoil SC1095
CP_L_H_60_lambda = @(CT_constant_T) (k.*(CT_constant_T.^(3/2))/sqrt(2)) + (sigma.*cd0./8);
CP_L_FF_60_lambda = @(mu, CT_constant_T)...
    ((k*CT_constant_T.^2)./(2.*mu)) + ((sigma.*cd0./8)*(1+(4.6*mu.^2)))...
    + ((1./2).*(f./A).*(mu.^3)) + (lambda_climb_FF.*CT_constant_T);   %%% POWER EQUATION IS RIGHT BUT GRAPH IS WRONG
%YF_W = @(CP_UH60_FF) CP_UH60_FF.*R./CT_constant_T./L_tr;

%% Finding Linear Coefficient of Power Required CP_UH60_FF based on lambda_TPP_UH60_FF
CP_UH60_H = CP_L_H_60_lambda(CT_constant_T(1));
CP_UH60_FF = CP_L_FF_60_lambda(mu(2:223), CT_constant_T(2:223));
CP_UH60_FF = [CP_UH60_H CP_UH60_FF];
%YF_W_FF = YF_W(CP_UH60_FF);


linear_inflow_power = CP_UH60_FF; 

end