clear, clc

u = 3986000;

r0 = [6320 0 7750];
v0 = [0 11 0];
mag_r0 = norm(r0);
mag_v0 = norm(v0);

h = cross(r0, v0);
mag_h = norm(h);

e = cross(v0, h)/u - r0/mag_r0;
mag_e = norm(e);

energy = mag_v0^2 / 2 - u / mag_r0;
a = -u / (2 * energy);

p = e / mag_e;
w = h / mag_h;
q = cross(w, p);

C_PF_ijk = [p; q; w];

r0_PF = C_PF_ijk * r0.';

theta0 = atan2(r0_PF(2), r0_PF(1));
theta0_deg = theta0 * 180 / pi;

E0 = 2 * atan((sqrt((1 - mag_e) / (1 + mag_e)))*tan(theta0 / 2));
t0 = (E0 - mag_e*sin(E0)) ./ (sqrt(u ./ a^3));

syms E
eqn = sqrt(u / a^3) * 600 == E - mag_e * sin(E);
solution = solve(eqn, E)

theta = 2*atan(sqrt((1 + mag_e) / (1 - mag_e))*tan(2.848/2));

deltaTheta = (theta0 - theta) * 180 / pi;