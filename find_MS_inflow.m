%% find_MS_inflow

function[lambda_i_MS] = find_MS_inflow(mu, psi, r, nu, alpha, CT_constant_T)
%% constants
c1 = zeros(1, 100);
epsilon = 1;
n = 1;

cn_even = @(n, nu, alpha)((-1)^((n-2)/2))*(3/4)* ((n+nu)/((n^2)-1)) ...
    * (((1-nu)/(1+nu))^(n/2)) * ...
    (((1+sin(alpha))/(1-sin(alpha)) )^(n/2));
for p = psi
    for s = r
        for m = 1:1:223
            for cn = 1:20
                if cn < 3
                    c1(1) = (3/4)*nu;
                    c1(2) = -(3*pi/16)*(sqrt(1-nu^2))*((1+sin(alpha))/(1-sin(alpha)))^(1/2);
                    c1(cn+1) = cn_even(n,nu,alpha);

                elseif (mod(cn,2)==0)
                    c1(cn+1) = cn_even(n,nu,alpha);

                elseif mod(cn,2)==1
                    c1(cn+1) = 0;
                end
            end
            while epsilon > 0.004
                lambda_i_MS(n+1) = ((2*CT_constant_T(1))/mu(m)) *( c1(1) + sum(c1(2:length(c1)).*cos(n*phi))) ;
                epsilon = abs(lambda_i_MS(n+1)-lambda_i_MS(n));
                n = n+1;
            end
        end
    end
end
end

% plot() = ;
