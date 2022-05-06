%% find_rotor
function[rotor, NACA23012, r_real, X, Y, Z] = find_rotor_23012(R, r)
import_23012 = readtable("C:\Users\Zarif Rahman"+...
    "\OneDrive\Desktop\Inflow Modeling\NACA23012.txt");
r_real = r.*ones(61, 100);
NACA23012 = table2array(import_23012);
x = NACA23012(:, 1);
y = NACA23012(:, 2);
[X, Y, Z] = meshgrid(NACA23012(:, 1), NACA23012(:, 2), r);
rotor = cat(4, X, Y, Z);
end