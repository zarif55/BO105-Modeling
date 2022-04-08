function [r_rotor] = find_rotor(psi_b, r, c)

x_rotor_c_4 = r'*cos(psi_b);
y_rotor_c_4 = r'*sin(psi_b);
x_rotor = x_rotor_c_4 - cos(psi_b).*sqrt(x_rotor_c_4.^2+((c/2)^2));
y_rotor = y_rotor_c_4 - sin(psi_b).*sqrt(y_rotor_c_4.^2+((c/2)^2));
z_rotor= zeros(size(x_rotor));
r_rotor = cat(3, x_rotor, y_rotor, z_rotor);