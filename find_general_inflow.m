%% find_general_inflow
function[R,A,A_blade,f,v_inf_mph,v_inf,c,sigma,omega,v_tip,alpha_d,alpha,...
    mu,T,CT_constant_T,lambda_hover,psi,e,r,Cl_alpha,nu,dr,...
    linear_inflow_power,MS_inflow_power,w] = find_general_inflow(Nb, R_i,...
    a, omega_i,A_i,A_blade_i,f_i,T_i,c_i,rho)
R = R_i*0.3048;
A = A_i*0.092903;
A_blade = A_blade_i*0.092903;
f = f_i*0.092903;
v_inf_mph = 0:1:222;
v_inf = v_inf_mph*0.44704;
c = c_i*0.0254;
sigma = (Nb*c)/(pi*R);
omega = (omega_i/60)*2*pi;
v_tip = omega*R;
alpha_d = -2;                                                          %degrees
alpha = -2*pi/180 * ones(1, 223);                                      %radians
mu = v_inf.*cos(alpha)/(omega*R);                                      % small angle approximation
T = T_i*4.4482216153;
CT_constant_T = (T/(rho*A*v_tip^2)) * ones(1, 223);
lambda_hover = sqrt(CT_constant_T/2);
psi = (0:1:360)*pi/180;
e = 0.17;
r = linspace(e, 1, 100);
Cl_alpha = a;
nu = sqrt(1-r.^2);
dr = (1-(e))/99;
linear_inflow_power = 0.01;
MS_inflow_power = 0.01;
w = 0.3;