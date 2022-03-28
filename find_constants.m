%% constants
function [Nb, N_psi, R_i, a, omega_i,A_i,A_blade_i,cd0,f_i,T_i,c_i,rho,k,...
    L_tr,x_cg,y_cg,h,vb, lock, theta_tw] = find_constants
Nb = 4;
N_psi = 360;
R_i = 26.85;                                                           %ft
a = 2*pi;                                                              %1/rad
omega_i = 258;                                                              %rpm
A_i = 2261.5;                                                          %ft^2
A_blade_i = 186.8;                                                     %ft^2
cd0 = 0.01;
f_i = 35.04;                                                           %ft
T_i = 20000;                                                           %lbs
c_i = 20.76;                                                           %in
rho = 1.225;                                                           %kg/m^3
k = 1.15;
L_tr = 45;  
x_cg = 0;
y_cg = 0;
h = 5;  
vb = 1;
lock = 8;
theta_tw = 0;