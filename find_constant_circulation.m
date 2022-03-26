%% find_constant_circulation

function[constant_circulation,constant_circulation_bar] = find_constant_circulation(rho, omega, R, Nb, T_MS_bar, v_tip, A)
%% constants
T = T_MS_bar;
CT = T/(0.5*rho*A*(v_tip)^2);

%% defining functions
constant_circulation_fn = @(T)(2*T) / (Nb*rho*omega*R^2);
constant_circulation_bar_fn = @(CT) 2*pi*CT/Nb;


%% 
constant_circulation = constant_circulation_fn(T);
constant_circulation_bar = constant_circulation_bar_fn(CT);