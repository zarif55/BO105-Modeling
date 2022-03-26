%% find_MS_1_inflow

function[cn_1] = find_MS_1_inflow(mu, psi, nu, alpha, CT_constant_T)
%% constants
%epsilon = 1;

cn_even = @(n, nu, alpha)((-1)^((n-2)/2))*(3/4)* ((n+nu)/((n^2)-1)) ...
    * (((1-nu)/(1+nu))^(n/2)) * ...
    (((1+sin(alpha))/(1-sin(alpha)) )^(n/2));
cn_1 = zeros(length(nu), 6);
% lambda_1_i_MS = zeros(length(psi), length(nu));
% for i = 1:length(psi)
    for j = 1:length(nu)
        n=5;
%         epsilon = 1;
        %         while epsilon > 0.004
        c1 = zeros(1, n+1);
        for cn = 0:n
            if cn == 0
                c1(1) = (3/4).*nu(j);
            elseif cn == 1
                c1(2) = -(3*pi/16).*(sqrt(1-nu(j).^2))*((1+sin(alpha))/(1-sin(alpha))).^(1/2);

            elseif (mod(cn,2)==0)
                c1(cn+1) = cn_even(n,nu(j), alpha);

            elseif mod(cn,2)==1
                c1(cn+1) = 0;
            end
        end
        %             lambda_1_i_MS_old = lambda_1_i_MS(i,j);
        %             lambda_1_i_MS(i,j) = ((2*CT_constant_T(1))/mu) *( c1(1) + sum(c1(2:length(c1)).*cos((1:n).*psi(j)))) ;
        %             epsilon = abs(lambda_1_i_MS(i,j)-lambda_1_i_MS_old);
%         n = n+1;
    cn_1(j,:) = c1;
    end
    %     end
% end

end