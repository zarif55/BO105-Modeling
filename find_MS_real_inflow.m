%% find_MS_real_inflow
function[w,lambda_MS_total] = find_MS_real_inflow(cn_1, cn_3, w, psi,nu, CT_constant_T, mu)

%% Adding padding of zeros
length_cn1 = length(cn_1);
length_cn3 = length(cn_3);
padding = [];
if length_cn1 > length_cn3
    padding = zeros(1, length_cn1 - length_cn3);
    cn_3 = [cn_3 padding];
elseif length_cn3 > length_cn1
    padding = zeros(1, length_cn3 - length_cn1);
    cn_1 = [cn1 padding];
end

%% defining functions
cn_fn = @(cn_1, cn_3, range_w) (range_w.*cn_1) + ((1-range_w).*cn_3);
cn_total_fn = @(cn) sum(cn);
lambda_MS_total = zeros(length(psi), length(nu));
w_fn = @(min_i) 0.3 + min_i.*(0.4/100);
epsilon = 1;

while epsilon > 0.004
    lambda_MS_old = lambda_MS_total;
    for range_w = w:0.0001:0.7
        cn = cn_fn(cn_1, cn_3, range_w);
        cn_total = cn_total_fn(cn)/length(range_w);
        [val, min_i] = min(cn_total);
        w = w_fn(min_i);
    end
    for i = 1:length(psi)
        lambda_MS_total(i,:) = ((2*CT_constant_T(1))./mu).*( cn(:,1) + ...
            sum(cn(:,2:size(cn, 2))...
            .*cos((1:(size(cn, 2)-1))...
            .*psi(i)), 2)) ;
    end
    epsilon = abs(lambda_MS_total-lambda_MS_old);
end
end