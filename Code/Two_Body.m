clear, clc
%% Constants
u_earth = 398600;
r_earth = 6371;

%% Initial Conditions
r0 = [8182 -6865.9 18281];
v0 = [.2 5 7];
mag_r0 = norm(r0);
mag_v0 = norm(v0);

%% Orbit Calculations
    %Momentum
    h = cross(r0, v0);
    mag_h = norm(h);

    %Eccentricity
    e = (cross(v0, h) / u_earth) - (r0 / mag_r0);
    mag_e = norm(e);

    %Perifocal calculations
    w = h / mag_h;
    p = e / mag_e;
    q = cross(w, p);

    %Direction Cosine Matric (PF to ijk)
    C_pf = [p; q; w];

    %radius in PF coords
    r0_pf = C_pf * r0.';
    r0_pf = r0_pf.';

    %Initial actual anomally 
    theta0 = atan2(r0_pf(2), r0_pf(1)) * 180 / pi;

    %Energy
    E = (mag_v0^2 / 2) - ( u_earth / mag_r0);

    %Semi-Major Axis
    a = -(u_earth / (2 * E));

    %Perigee and Apogee Radii
    r_perigee = a * (1 - mag_e) * p;
    r_p = norm(r_perigee);
    r_apogee = a * (1 + mag_e) * (-p);
    r_a = norm(r_apogee);

%% Polar Plot
    %Orbit Calculations from Actual Anomollay
    theta = 0:.01:2*pi;
    radius = mag_h^2 ./ (u_earth * (1 + mag_e * cos(theta)));
    
    %Setting up and plotting orbit
    figure(1)
    polarplot(theta, radius, '--k')
    hold on
    r_1 = ones(1, length(theta)) * r_earth;
    polarplot(theta, r_1, 'b')
    legend('Orbit', 'Earth')

%% 3D Plot
   %Cartesian Coordinates of the Orbit
    x = radius.*cos(theta);
    y = radius.*sin(theta);
    z = y * 0;
    
    %Setting up and plotting orbit
    figure(2)
    plot3(x, y, z)
    grid on
    hold on
    plotEarth