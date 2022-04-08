clc
clear
close all
%% Velocity Profile x - x^3;
%% constants
x_se = 0.655;
V_inf = 1; %m/s
nu = 1.48*10^-5; %m^2/s
mu = 1.81*10^-5;
L = 1;
h = 0.1;
dx = 0.0001; %m
dy = 0.0001; %m
Nx = (L/dx) + 1;
Ny = (h/dy) + 1;
N_th = (0.1/dx)+1;
%% arrays
x_pts = 0:1:Nx;
y_pts = 0:1:Ny;
x = x_pts.*dx;
y = y_pts.*dy;
U_inf = x - x.^3 + 0.001;
U_inf_spacer = ones(1, length(x_pts));
%% problem 1
[Re_x_blasius, delta_blasius, delta_star_blasius, theta_blasius,...
u_blasius, v_blasius, alpha_b, beta_b, cf_x_blasius]...
= find_curve(dx, dy, x,y, nu, U_inf_spacer, N_th, Ny);
[Re_x_blasius2, delta_blasius2, delta_star_blasius2, theta_blasius2,...
u2_blasius, v2_blasius, alpha_b2, beta_b2, cf2_x_blasius]...
= find_curve(dx, dy, x,y, nu, U_inf_spacer, 10, Ny);
[Re_x_blasius3, delta_blasius3, delta_star_blasius3, theta_blasius3,...
u3_blasius, v3_blasius, alpha_b3, beta_b3, cf3_x_blasius]...
= find_curve(dx, dy, x,y, nu, U_inf_spacer, 10000, Ny);
Re_a = U_inf*L/nu;
delta_a = 4.9*L/sqrt(Re_x_blasius3(10000));
delta_star_a = 1.72*L/sqrt(Re_x_blasius3(10000));
theta_a = 0.664/2*L/sqrt(Re_x_blasius3(10000));
% plot
figure(1)
hold on, grid on
title("Parallel Velocity U_inf = 1")
xlim([0, 1.25]), ylim([0,27])
plot(u_blasius(:)/U_inf_spacer(1), y/delta_blasius3(10000))
plot(u2_blasius(:)/U_inf_spacer(1), y/delta_blasius3(10000))
plot(u3_blasius(:)/U_inf_spacer(1), y/delta_blasius3(10000))
legend("0.1 m","0.00009 m","0.9999 m")
xlabel("u/U_inf"); ylabel("Eta")
figure(2)
hold on, grid on
title("Perpendicular Velocity U_inf = 1")
plot(v_blasius(:)/U_inf_spacer(1), y/delta_blasius3(10000))
plot(v2_blasius(:)/U_inf_spacer(1), y/delta_blasius3(10000))
plot(v3_blasius(:)/U_inf_spacer(1), y/delta_blasius3(10000))
legend("0.1 m","0.00009 m","0.9999 m")
xlabel("v/U_inf"), ylabel("Eta")
figure(3)
hold on, grid on
title("Freestream Velocity")
plot(x, U_inf(:), '-')
xlabel("x"), ylabel("U_inf")
% display
disp("delta star blasius = "); disp(delta_star_blasius(1000))
disp("delta star blasius analytic = "); disp(delta_star_a)
disp("delta blasius = "); disp(delta_blasius(1000))
disp("delta blasius analytic = "); disp(delta_a)
disp("theta blasius = "); disp(theta_blasius(1000))
disp("theta blasius analytic = "); disp(theta_a)
%% problem 2
[Re1_x, delta1, delta1_star, theta1, u, v1, alpha1, beta1, cf1_x]...
= find_curve(dx, dy, x,y, nu, U_inf, 100, Ny);
[Re2_x, delta2, delta2_star, theta2, u2, v2, alpha2, beta2, cf2_x]...
= find_curve(dx, dy, x,y, nu, U_inf, 1000, Ny);
[Re3_x, delta3, delta_star3, theta3, u3, v3, alpha3, beta3, cf3_x]...
= find_curve(dx, dy, x,y, nu, U_inf, 10000, Ny);
[Re_x, delta, delta_star, theta, u1, v, alpha, beta, cf_x]...
= find_curve(dx, dy, x,y, nu, U_inf, Nx, Ny);
disp("stable")
figure(4)
hold on, grid on
title("Parallel Velocity U_inf = x - x^3 + 0.001")
plot(u1(:)/U_inf(1000), y/delta3(10000))
plot(u2(:)/U_inf(10000), y/delta3(10000))
plot(u3(:)/U_inf(10000), y/delta3(10000))
legend("0.0099 m","0.0999 m","0.9999 m")
figure(5)
hold on, grid on
title("Perpendicular Velocity U_inf = x - x^3 + 0.001")
plot(v1(:)/U_inf(10000), y/delta3(10000))
plot(v2(:)/U_inf(10000), y/delta3(10000))
plot(v3(:)/U_inf(10000), y/delta3(10000))
legend("0.0099 m","0.0999 m","0.9999 m")
figure(6)
loglog(x, cf_x.*sqrt(Re_x))
roots = x((cf_x<75)&(cf_x>72));
%% functions
function[Re_x, delta, delta_star, theta, u, v, alpha, beta, cf_x] = ...
    find_curve(dx,dy, x,y, nu, U_inf, Nx, Ny);
%zeros
u = ones(1, length(y));
v = zeros(1, length(y));
alpha = zeros(1,length(y));
beta = zeros(1,length(y));
Re_x = zeros(1, length(x));
delta = zeros(1, length(x));
delta_star = zeros(length(y));
theta = zeros(length(y));
cf_x = ones(1, length(x));
u_new = zeros(1,length(y));
for m=1:Nx
u(1) = 0;
u(end) = U_inf(m);
v(1) = 0;
for n=2:Ny
alpha(n) = (nu*dx)/((u(n)*(dy^2)));
beta(n) = (dx*v(n))/((2*u(n)*(dy)));
if alpha(n)<1/2 && beta(n)<=alpha(n)
u_new(n) = ((alpha(n) - beta(n)).*(u(n+1))) ...
+ ((1-(2.*alpha(n))).*u(n)) +...
((alpha(n)+beta(n)).*(u(n-1)))...
+( (((U_inf(m+1)).^2)-(U_inf(m)).^2) / (2.*u(n)));
v(n) = v(n-1)...
-((dy./(2.*dx)).*(u_new(n)...
-u(n)+u_new(n-1)-u(n-1) ));
Re_x(m+1) = U_inf(m+1)*x(m+1)/nu;
delta(m+1) = x(m+1)/sqrt(Re_x(m+1));
delta_star(n) = sum(y(n)...
*(1-(u_new(n)./U_inf(m+1))));
theta(n) = sum(y(n)*(u_new(n)./U_inf(m+1))...
*(1-(u_new(n)./U_inf(m+1))));
end
du(n) = u(n+1)-u(n-1);
u(n) = u_new(n);
end
cf_x(m+1) = (nu*(du(2)/dy))/(0.5*(U_inf(m+1)^2));
end
cf_x(1) = 0;
delta(1) = 0;
end