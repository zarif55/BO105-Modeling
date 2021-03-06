%% find_acoustics
function [P_prime, P_prime_inter, t_real, P_prime_old]= find_acoustics(R, r, rotor, c0, rho, u_t, L_5, dr)
P_prime_old = [];
t_old = [];
for tau = 0:0.1:5
x_obs = rotor;
y_pos = 4;
r_ac = x_obs - y_pos; %%should be a unit vector
r_mag = r_ac./(sqrt(sum(((r_ac).^2), 4)));
t = tau + norm(squeeze(r_mag(1, 1, 1, :)))/c0;
dtau = 0.1;
v_surface = u_t;
M = v_surface(1 ,:)/c0;
Mr = dot(M, squeeze(r_mag(1, 1, :, 3))');
L = L_5./(R.*dr);
L_m = dot(L, M); %% units
L_r = dot(L, squeeze(r_mag(1, 1, :, 3))'); %% make sure
L_dot_r = L_r./(dtau);
M_r_dot = Mr./(dtau);
dS = dr*R;
doppler = 1/(abs(1-Mr));
P_prime = (1/(4*pi))*...
    ((dS./c0)*((rho*(L_dot_r)) ./ ((squeeze(r_mag(1, 1, :, 3))').*(abs(1-Mr)).^2))...
    +((dS)*((L_r-L_m) ./ (((squeeze(r_mag(1, 1, :, 3))').^2).*(abs(1-Mr)).^2)))...
    +((dS./c0)*((L_r*((squeeze(r_mag(1, 1, :, 3))'.*M_r_dot)+...
    (c0.*Mr)-(c0*M.^2)))) ./ ((squeeze(r_mag(1, 1, :, 3))'.^2).*(abs(1-Mr)).^3)));
P_prime_old = [P_prime_old mean(P_prime)];
t_old = [t_old t];
end
t_real = 0:0.1:5;
 %= interp1(P_prime)
% %%choose t evenly spaced define seperately 
%P_prime_sum = 
% %%% 10 Pa on the order of 1 or less Pa 
P_prime_inter = interp1(t_old, P_prime_old, t_real);
%%P_prime_sum = sum(P_prime_inter(1,:));
end