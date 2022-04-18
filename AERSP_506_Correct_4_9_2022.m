clc
clear
close

Nb  = 4;
R = 20; %ft
c  = 1; %ft
v_tip = 650; %ft/sec
W = 20000; %lbs
h  = 5; %ft
a  = 5.7; %per
x_cg = 0; %ft
vb = 1.04;
lock = 10;
cd0 = 0.012;
L_tr = 25; %ft ;
f = 15; %ft^2
y_cg = 0;
k = 1.15;
rho = 0.002378;


%% General
A = pi*R^2;
omega = v_tip/R;
psi = 0:pi/180:2*pi;
sigma = Nb*c/pi/R;
mu_FF = 200/v_tip ;
mu_HtFF = 0:0.1:0.5;
f_HtFF = [f*0.5, f, f*1.5];

[plot_data_hover] = trim_506(R, v_tip, W, h, a, x_cg, vb, lock, cd0,...
    L_tr, f, y_cg, k, rho, omega, psi, 0, A, sigma);
[plot_data_FF] = trim_506(R, v_tip, W, h, a, x_cg, vb, lock, cd0,...
    L_tr, f, y_cg, k, rho, omega, psi, mu_FF, A, sigma);

for m = 1:6
    [plot_data_HtFF_1] = trim_506(R, v_tip, W, h, a, x_cg, vb, lock, cd0,...
        L_tr, f_HtFF(1), y_cg, k, rho, omega, psi, mu_HtFF(2), A, sigma);
end

for m = 1:6
    [plot_data_HtFF_2] = trim_506(R, v_tip, W, h, a, x_cg, vb, lock, cd0,...
        L_tr, f_HtFF(2), y_cg, k, rho, omega, psi, mu_HtFF(m), A, sigma);
end

for m = 1:6
    [plot_data_HtFF_3] = trim_506(R, v_tip, W, h, a, x_cg, vb, lock, cd0,...
        L_tr, f_HtFF(3), y_cg, k, rho, omega, psi, mu_HtFF(m), A, sigma);
end






function[plot_data] = trim_506(R, v_tip, W, h, a, x_cg, vb, lock, cd0,...
    L_tr, f, y_cg, k, rho, omega, psi, mu, A, sigma)
M_yF = 0;
M_xF = 0;
CT = W/(rho*(pi*R^2)*(v_tip)^2);
plot_data = zeros(1, 8);
er = 0.005;
lambda(1) = k*sqrt(CT/2);
er_lambda = 1;
while er_lambda > er
    lambda_old = lambda;
    lambda = (CT/(2*sqrt((mu^2)+...
        (lambda)^2))) + ((1/2)*(f/A)*(mu^3/CT));
    er_lambda = lambda_old - lambda;
end

CP = k*lambda*CT + sigma*cd0/8;
YF_W = CP*R/CT/L_tr;
CH_TPP = 0;
CY_TPP = 0;
er_beta_1c = 1;
er_beta_1s = 1;

beta_1c = 1;
beta_1s = 1;

while [er_beta_1c, er_beta_1s] > er%change to while loop, errors in a vector
    beta_1c_old = beta_1c;
    beta_1s_old = beta_1s;

    beta_1c = (-x_cg/h + M_yF / (h * W) + (CH_TPP / CT)) / (1 + ((vb^2 - 1)...
        / lock) / ((h * 2 * CT) / (R * sigma * a)));
    beta_1s = (y_cg / h - M_xF / (h * W) + (CY_TPP / CT)) / (1 + ((vb^2 - 1)...
        / lock) / (h * 2 * CT / (R * sigma * a)));
    alpha_shaft = (x_cg / h - M_yF / (h .* W) + ((vb^2 - 1) / lock)...
        .*(CH_TPP / CT) / (h * 2 * CT / (R * sigma * a))) / (1 + ((vb * vb - 1)...
        / lock) / (h .* 2 .* CT / (omega .* sigma .* a)))...
        + 0.5 .* f .* mu .* mu / (A .* CT);
    phi_shaft = (y_cg / h - M_xF / (h .* W) - ((vb^2 - 1) / lock)...
        .* (CH_TPP / CT) / (h .* 2 .* CT / (sigma .* R .* a))) / (1 + (vb^2 - 1) .* R .* sigma .* a / (lock .* 2 .* h .* CT)) - YF_W;
    theta_knot = ((6 .* CT / (sigma .* a)) .* (1 + 3 .* mu.^2 / 2)...
        + lambda.* 1.5 .* (1 - 0.5 .* mu.^2) + 3 * mu .* beta_1s .* (vb^2 - 1) / lock) / (1 - mu.^2 + 0.25 * 9 * (mu.^4));
    theta_1s = - beta_1c + (1 ./ (1 + 1.5 .* mu.^2)) .* (8 .* (vb^2 - 1)...
        .* beta_1s / lock - (8 .* mu / 3) .* (theta_knot - 3 .* lambda .* 0.25));
    beta_knot = (lock / (vb^2)) .* (theta_knot .* (1 + mu.^2) / 8 + mu.* (beta_1c + theta_1s) / 6 - lambda / 6);
    theta_1c = beta_1s + (1 / (1 + 1.5.* mu.^2)) .* ((4 .* mu .* beta_knot ./ 3) + (8 ./ lock) .* (vb^2 - 1) .* beta_1c);
    CH_TPP = sigma*a*0.5*(1/2.*mu.* lambda.*theta_knot- ...
        (beta_knot.*theta_1c/6)...
        +0.25.*theta_1s.*lambda...
        +0.25.*mu.*beta_knot.^2) ...
        + sigma.*cd0.*mu./ 4;
    CY_TPP = - sigma .* a .* 0.5 .* (3 .* mu .* beta_knot .* theta_knot ./ 4 + 0.25 .* theta_1c .* lambda...
        + beta_knot .* theta_1s .* (1 + 3 .* mu.^2) / 6 - 3 .* mu .* beta_knot .* lambda ./ 2);
    er_beta_1c = beta_1c - beta_1c_old/beta_1c;
    er_beta_1s = beta_1s - beta_1s_old/beta_1s;
    plot_data(1) = beta_1c;
    plot_data(2) = beta_1s;
    plot_data(3) = theta_1c;
    plot_data(4) = theta_1s;
    plot_data(5) = alpha_shaft;
    plot_data(6) = phi_shaft;
    plot_data(7) = beta_knot;
    plot_data(8) = theta_knot;
end
end


