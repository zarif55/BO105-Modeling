%% find_duhmel_integral
%solve duhmels integral using trapezoidal method
function [U] = find_duh(Thrust)
%% constants
m = 1000;
zeta = 0.05;
f = 1.5;
wn = 2*pi*f;
wd = wn*sqrt(1-zeta^2);


t = 0:0.1:100;
freq = 1;
wf = 2*pi*freq;
Thurst_duhmel = Thrust*sin(wf*t);

%% initialize
U = zeros(length(t));
A = 0; 
B = 0;
A_integral(1) = 0;

for n = 2:length(t)
A_old = A_integral(n-1);
A_integral(n)= (exp(zeta*wn*t(n)))*Thurst_duhmel*cos(wd*t(n));
Area_A = trapz(A_old(n), A_integral(n));
A = (1/(wd*m))*(A(n-1) + Area_A(n));

B_old = B_integral(n-1);
B_integral= (exp(zeta*wn*t(n)))*Thurst_duhmel*sin(wd*t(n));
Area_B = trapz(B_old(n), B_integral(n));
B = (1/(wd*m))*(B(n-1) + Area_B(n));

U = (A(n)*(exp(-zeta(wn)*t(n)))*sin(wd*t(n)))...
    - B(n)*(exp(-zeta*wn*t(n)))*(cos(wd*t(n)));
end
end