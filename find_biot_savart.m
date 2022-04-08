%% find_biot_savart
function[V, lambda_wake] = find_biot_savart(circulation_v,...
    v_tip, r_wake, psi_w, psi_b, r_rotor, r) 

%% constants

%% zeros
lambda_wake = zeros(length(r),length(psi_b),3);
%% defining functions
%%
for PB = 1:length(psi_b)
    for R_R = 1:length(r)
        V = zeros(length(psi_w),3);
        for PW = 1:length(psi_w)-1
            r1 = r_wake(PW,PB,:) - (r_rotor(R_R,PB,:));
            r2 = r_wake(PW+1,PB,:) - (r_rotor(R_R,PB,:));
            lv = r2 - r1;
%             disp("display 3")
            h = abs(cross(r1, lv));
%             disp("display 4")
            V(PW, :) = (circulation_v/(4.*pi.*h)).*...
                (dot(lv, r1-r2).*cross(r1, r2));
        end
        V_total = sum(V);
        lambda_wake(R_R, PB, :) = V_total./v_tip;
    end
end