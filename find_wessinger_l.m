%% find_wessinger_l

%% near wake 
function[control_pt] = find_wessinger_l(c, xv, yv, zv, r_wake) 
del_psi = pi/12;
del_psi2 = pi/6;
psi_b = 0:del_psi:90*pi;
psi_w = 0:del_psi2:180*pi;
control_pt_x_fn = @(xv) xv + c/2*(cos(psi_w))
control_pt_y_fn = @(yv) yv + c/2*(sin(psi_w))

control_pt_x = control_pt_x_fn(xv)
control_pt_y = control_pt_y_fn(yv)
control_pt = [control_pt_x, control_pt_y, zv]


