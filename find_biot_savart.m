%% find_biot_savart
function[V, lambda_wake, r1, r2] = find_biot_savart(circulation_v,...
    v_tip, r_wake, psi_w, psi_b, r_rotor, r) %% first problem is how r wake is set up

%% constants

%% zeros
lambda_wake = zeros(length(r),length(psi_b),3);
%% defining functions
lv_fn = @(r1, r2) r2 - r1;
h_fn = @(r1, lv) abs(cross(r1, lv));
V_fn = @(circulation_v, r1, r2, h, lv) ...
    (circulation_v/(4.*pi.*h)).*...
    (dot(lv, r1-r2).*cross(r1, r2));

%%
for PB = 1:length(psi_b)
    for R_R = 1:length(r)
        V = zeros(length(psi_w),3);
        for PW = 1:length(psi_w)-1
            r1 = r_wake(PW,PB,:) - (r_rotor(R_R,PB,:));
            r2 = r_wake(PW+1,PB,:) - (r_rotor(R_R,PB,:));
            lv = lv_fn(r1, r2);
            h = h_fn(r1, lv);
            V(PW, :) = V_fn(circulation_v, r1, r2, h, lv);
        end
        V_total = sum(V);
        lambda_wake(R_R, PB, :) = V_total./v_tip;
    end
end