function [Cl, L, phi_wp_jones] = find_beddoes_mdpt(A, v, rho, alpha, theta)
% approximate solutions from W.P. Jones
%% constants
A1 = 0.165;
b1 = 0.0455;
A2 = 0.335;
b2 = 0.3;

%% zeros
phi_wp_jones = zeros(1, 100);

X = 0;
Y = 0;
ds = 1;
for s = 2:length(alpha)
X_old = X;
Y_old = Y;
d_alpha = alpha(s) - alpha(s-1);
X = (X_old*exp(-b1*ds)) + (A1*d_alpha*exp(-b1*ds));
Y = (Y_old*exp(-b2*ds)) + (A2*d_alpha*exp(-b2*ds));
alpha_eq = alpha - X - Y;
cl_eq = 2*pi*alpha_eq;
epsilon_mdpt = 2 - ((b1*ds*exp(-b1*ds))/(1-exp(-b1*ds)))...
- ((b2*ds*exp(-b2*ds))/(1-exp(-b2*ds)));
phi_wp_jones(s) = 1 - (A1*exp(-b1*s)) - (A2*exp(-b2*s));
end
Cl = cl_eq(1, :);
L = 0.5.*rho.*(v.^2).*A.*Cl;

