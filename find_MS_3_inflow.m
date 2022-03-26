%% find_MS_3_inflow

function[cn_3] = find_MS_3_inflow(mu, psi, nu, alpha, CT_constant_T)
%% constants

%epsilon = 1;

cn_even = @(n, nu, alpha)((-1)^((n-2)/2)).*(15/8).*  ...
    ((((n+nu)/((n^2)-1))* (((9*(nu.^2) + ...
    (n^2) - 6))/((n^2)-9)))+(3*nu/((n^2)-9))).*...
    (((1-nu)/(1+nu))^(n/2)) * ...
    (((1+sin(alpha))/(1-sin(alpha)) )^(n/2));
lambda_3_i_MS = zeros(length(psi), length(nu));
cn_3 = zeros(length(nu), 6);
% for i = 1:length(psi)
for j = 1:length(nu)
    n=5;
    %         epsilon = 1;
    %         while epsilon > 0.004
    c1 = zeros(1, n+1);
    for cn = 0:n
        if cn == 0
            c1(1) = (15/8)*nu(j).*(1-(nu(j)).^2);
        elseif cn == 1
            c1(2) = -(15*pi/256)*(5-(9*(nu(j)).^2)).*(sqrt(1-(nu(j)).^2)).*(((1+sin(alpha))/(1-sin(alpha))).^0.5);

        elseif cn ==  3
            c1(4) = (45*pi/256)*((1-((nu(j)).^2)).^(3/2)).*(((1+sin(alpha))/(1-sin(alpha))).^1.5);

        elseif (mod(cn,2)==0)
            c1(cn+1) = cn_even(n,nu(j), alpha);

        elseif mod(cn,2)==1
            c1(cn+1) = 0;
        end
        %             end
        %             lambda_3_i_MS_old = lambda_3_i_MS(i,j);
        %             lambda_3_i_MS(i,j) = ((2.*CT_constant_T(1))/mu).*( c1(1) + sum(c1(2:length(c1)).*cos((1:n).*psi(j)))) ;
        %             epsilon = abs(lambda_3_i_MS(i,j)-lambda_3_i_MS_old);
        %             n = n+1;
        %         end
        cn_3(j,:) = c1;
    end
end

end