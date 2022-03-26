%% find_biot_savart
function[V, lambda_wake, r1, r2] = find_biot_savart(circulation_v,...
    v_tip, r_wake, psi_w, psi_b) %% first problem is how r wake is set up

%% constants
% rv = 1;
% lv = 30*pi*rv;
% %dlv = rv*(pi/12)/(2*pi);
r_biot = @(r_wake) r_wake;

%% defining functions
r1_fn = @(r_wake) r_wake;
r2_fn = @(r_wake) r_wake;
lv_fn = @(r1, r2) r2 - r1;
h_fn = @(r_biot, lv) abs(cross(r_biot, lv));
V_fn = @(circulation_v, r1, r2, h, lv) (circulation_v/(4*pi*h))*(dot(lv, r1-r2).*cross(r1, r2));

%%
for n = 1:length(psi_b)
r1(n) = r1_fn(r_wake(n));
r2(n) = r2_fn(r_wake(n+1));
lv(n) = lv_fn(r1(n), r2(n));
h(n) = h_fn(r_biot(n), lv(n));
V(n) = V_fn(circulation_v, r1(n), r2(n), h(n));
lambda_wake(n) = V(n)/v_tip;
end