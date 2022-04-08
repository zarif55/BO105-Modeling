%% find_prescribed_wake
function[xv, yv, zv, mu_z, mu_x, E, r_wake] =...
    find_prescribed_wake(lambda_TPP_linear_it, mu, alpha, v_tip,...
    r, psi, psi_w, psi_b)

%% constants
lambda_io = lambda_TPP_linear_it;
X = atan(-lambda_io./mu);
E = X/2;    
mu_x = mu.*cos(alpha)./v_tip;
mu_z = mu.*sin(alpha)./v_tip;
rv = 1;

%% defining functions
zv1_fn = @(lambda_io, mu_z, psi_b, psi_w, mu_x, E, y0, x0) mu_z*psi_w - ...
    lambda_io.*psi_w.*(1+(E.*(x0+(0.5.*mu_x.*psi_w)-(abs(y0).^3))));
zv2_fn = @(lambda_io, mu_z, psi_b, psi_w, mu_x, E, y0, x0) mu_z.*psi_w - ...
    ((2.*lambda_io.*psi_w).*(1-(E.*(abs(y0).^3))));
zv3_fn = @(lambda_io, mu_z, psi_b, psi_w, mu_x, E, y0, x0) mu_z.*psi_w - ...
    ((2.*lambda_io.*x0./mu_x).*(1-(E.*(abs(y0).^3))));
%%
%%zv = %% initialize
zv1_count = 0;
zv2_count = 0;
zv3_count = 0;
xv = zeros(length(psi_w), length(psi_b));
zv = zeros(length(psi_w), length(psi_b));
for n = 1:length(psi_w)
    for m = 1:length(psi_b)
        x0(n, m) = rv.*cos(psi_b(m)-psi_w(n));
        y0(n, m) = rv.*sin(psi_b(m)-psi_w(n));
        yv(n, m) = y0(n, m);
        xv(n, m) = x0(n, m) + mu_x*psi_w(n);
        if xv(n, m) < -cos(psi_b(m)-psi_w(n))
            zv(n, m) = zv1_fn(lambda_io, mu_z, psi_b(m), psi_w(n), mu_x, E, yv(n, m), xv(n, m));
            zv1_count = zv1_count + 1;
        elseif xv(n, m) > cos(psi_b(m)-psi_w(n))
            zv(n, m) = zv2_fn(lambda_io, mu_z, psi_b(m), psi_w(n), mu_x, E, yv(n, m), xv(n, m));
            zv2_count = zv2_count + 1;
        else
            zv(n, m) = zv3_fn(lambda_io, mu_z, psi_b(m), psi_w(n), mu_x, E, yv(n, m), xv(n, m));
            zv3_count = zv3_count + 1;
        end
    end
end
zv(1, :) = 0;
% disp(size(xv))
% disp(size(yv))
% disp(size(zv))
% disp(zv1_count)
% disp(zv2_count)
% disp(zv3_count)
r_wake = cat(3, xv, yv, zv);


