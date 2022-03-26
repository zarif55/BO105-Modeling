%% find_linear_inflow function
% mu: advance ratio
% alpha: angle of attack
% CT_constant_T: thrust coefficient
function[lambda_TPP_UH60_FF,lambda_i_TPP_UH60_FF, lambda_output] = find_linear_inflow(mu, alpha, CT_constant_T)
%% defining constants
%alpha = -2 * ones(1, 223);

%% defining functions

lambda_H_60 = @(CT_constant_T) sqrt(CT_constant_T./2);
lambda_i_FF_60 = @(lambda_FF_60_old) CT_constant_T./ ...
    (2*sqrt((((mu.*sin(alpha)-lambda_FF_60_old)).^2)...
    +(mu.*cos(alpha)).^2));
lambda_FF_60 = @(lambda_i_FF_60) (-mu.*sin(alpha) + lambda_i_FF_60);

%% define lambda_TPP_UH60_FF
% assuming CT_constant_T is 1 x N
% lambda_TPP_UH60_FF would be 1 x N
lambda_TPP_UH60_FF = lambda_H_60(CT_constant_T);
lambda_output = lambda_TPP_UH60_FF;
%% Finding Linear inflow of lambda_TPP_UH60_FF
% Using 100 iterations to find 223 different values of mu
% mu should be 1 x 223
% lambda_TPP_UH60_FF should be 1 x 223 as well
for n=1:100
    lambda_TPP_UH60_FF_old = lambda_TPP_UH60_FF; %lambda_TPP_UH60_FF_old is 1 x N
    lambda_i_TPP_UH60_FF= lambda_i_FF_60(lambda_TPP_UH60_FF_old);
    lambda_TPP_UH60_FF = lambda_FF_60(lambda_i_TPP_UH60_FF);
end




