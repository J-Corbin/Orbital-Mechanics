%% Problem 6.14

clear, clc

u = 398600;

rp1 = 6674;
ra1 = 6680;

rp2 = 6637;
ra2 = 6669;

rp3 = 6637;
ra3 = 6637;

rp4 = 6633;
ra4 = 6572;

a1 = (rp1 + ra1) / 2;
a2 = (rp2 + ra2) / 2;
a3 = (rp3 + ra3) / 2;
a4 = (rp4 + ra4) / 2;

e1 = (ra1 - rp1) / (ra1 + rp1);
e2 = (ra2 - rp2) / (ra2 + rp2);
e3 = (ra3 - rp3) / (ra3 + rp3);
e4 = (ra4 - rp4) / (ra4 + rp4);

h1 = sqrt(a1 * u * (1 - e1^2));
h2 = sqrt(a2 * u * (1 - e2^2));
h3 = sqrt(a3 * u * (1 - e3^2));
h4 = sqrt(a4 * u * (1 - e4^2));


%%
vP1 = h1 / rp1;
vA2 = h2 / ra2;
vA1 = h1 / ra1;
vP2 = h2 / rp2;

vPT1 = sqrt((2 * u * ra2) / (rp1 * (rp1 + ra2)));
vAT2 = sqrt((2 * u * rp1) / (ra2 * (rp1 + ra2)));

vAT1 = sqrt((2 * u * ra1) / (rp2 * (ra1 + rp2)));
vPT2 = sqrt((2 * u * rp2) / (ra1 * (ra1 + rp2)));

dV1_1 = vPT1 - vP1 + vA2 - vAT2;
dV1_2 = vAT1 - vA1 + vP2 - vPT2;

if (dV1_2 / dV1_1 > 1)
   delta_V1 = dV1_1; 
else
    delta_V1 = dV1_2;
end

%%
vP2 = h2 / rp2;
vA3 = h3 / ra3;
vA2 = h2 / ra2;
vP3 = h3 / rp3;

vPT2 = sqrt((2 * u * ra3) / (rp2 * (rp2 + ra3)));
vAT3 = sqrt((2 * u * rp2) / (ra3 * (rp2 + ra3)));

vAT2 = sqrt((2 * u * ra2) / (rp3 * (ra2 + rp3)));
vPT3 = sqrt((2 * u * rp3) / (ra2 * (ra2 + rp3)));

dV2_1 = vPT2 - vP2 + vA3 - vAT3;
dV2_2 = vAT2 - vA2 + vP3 - vPT3;

if (dV2_2 / dV2_2 > 1)
   delta_V2 = dV2_1; 
else
    delta_V2 = dV2_2;
end

%%
vP3 = h3 / rp3;
vA4 = h4 / ra4;
vA3 = h3 / ra3;
vP4 = h4 / rp4;

vPT3 = sqrt((2 * u * ra4) / (rp3 * (rp3 + ra4)));
vAT4 = sqrt((2 * u * rp3) / (ra4 * (rp3 + ra4)));

vAT3 = sqrt((2 * u * ra3) / (rp4 * (ra3 + rp4)));
vPT4 = sqrt((2 * u * rp4) / (ra3 * (ra3 + rp4)));

dV3_1 = vPT3 - vP3 + vA4 - vAT4;
dV3_2 = vAT3 - vA3 + vP4 - vPT4;

if (dV3_2 / dV3_2 > 1)
   delta_V3 = dV3_1; 
else
    delta_V3 = dV3_2;
end

deltaV = delta_V1 + delta_V2 + delta_V3