%% Problem 6.13

clear, clc

u = 398600;

rp1 = 7000;
e1 = 0.3;

rp2 = 32000;
e2 = 0.5;

h1 = sqrt(rp1 * u * (1 + e1));

vP1 = u / h1 * (1 + e1);

ra2 = (-rp2 - e2*rp2) / (e2 - 1);

vPT1 = sqrt((2 * u * ra2) / (rp1 * (rp1 + ra2)));

VP = vPT1 - vP1;

h2 = sqrt(rp2 * u * (1 + e2));

vA2 = u / h2 * (1 - e2);

vAT2 = sqrt((2 * u * rp1) / (ra2 * (rp1 + ra2)));

VA = vA2 - vAT2;

delta_V = VP + VA

a = (rp1 + ra2) / 2;

TOF = pi * a^(3/2) / sqrt(u) / 3600


ra1 = (-rp1 - e1 * rp1) / (e1 - 1);

eT = (ra1 - rp2) / (ra1 + rp2);
hT = sqrt((ra1 + rp2) / 2 * u * (1 - eT^2));

vPT2 = hT / ra1;

vA1 = h1 / ra1;

vP2 = h2 / rp2;

vAT2 = hT / rp2;

delta_V2 = vPT2 - vA1 + vP2 - vAT2

TOF2 = pi * ((ra1 + rp2)/2)^(3/2) / sqrt(u) / 3600