function [sat_orbit, position] = satellite(e, a, Omega, omega, i, theta, t)

    global u
    
    h = sqrt(u * a * (1-e^2));
    
    %Computing direction cosine matrix from perifocal to N for unknown object
    C_N_PF_obj = [-sin(Omega)*cos(i)*sin(omega)+cos(Omega)*cos(omega), -sin(Omega)*cos(i)*cos(omega)-cos(Omega)*sin(omega), sin(Omega)*sin(i);...
        cos(Omega)*cos(i)*sin(omega)+sin(Omega)*cos(omega), cos(Omega)*cos(i)*cos(omega) - sin(Omega)*sin(omega), sin(i)*cos(omega);...
        sin(i)*sin(omega), sin(i)*cos(omega), cos(i)];

    %Calculating initThial eccentric anomaly from given theta initial
    E0 = 2 * atan2((tan(theta / 2) / sqrt((1 + e)/(1 - e))), 1);
    
    %Calculating time since pergogee t0
    t0 = (E0 - e*sin(E0)) / sqrt(u / a^3);
    
    %Calculating the mean anomalies at times
    M1 = sqrt(u / a^3) * (t(2) - t0);
    M2 = sqrt(u / a^3) * (t(3) - t0);
    
    %Calculating eccentric anomaly at times
    E1 = kepler_E(e, M1);
    E2 = kepler_E(e, M2);
    
    %Calculating actual anomaly values at times
    theta1 = 2 * atan2(sqrt((1 + e)/(1 - e)) * tan(E1 / 2), 1);
    theta2 = 2 * atan2(sqrt((1 + e)/(1 - e)) * tan(E2 / 2), 1);
    
    %For loop to calculate X,Y,Z positions of the estimated orbit
    counter = 1;
    coordinate = zeros(3, 629);
    for angle = 0:.01:(2*pi)
        radius_polar = [h^2/(u * (1 + e * cos(angle)));0;0];
        C_PF_polar_obj = [cos(angle), sin(angle), 0;...
        -sin(angle), cos(angle), 0;...
        0, 0, 1].';
       radius_pf = (C_PF_polar_obj * radius_polar); 
       radius = (C_N_PF_obj * radius_pf).';
       for z = 1:3
           coordinate(z, counter) = radius(z);
       end
       counter = counter + 1;
    end
    
    %Determing position vector at the three theta values
    radius_polar = [h^2/(u * (1 + e * cos(theta)));0;0];
        C_PF_polar_obj = [cos(theta), sin(theta), 0;...
        -sin(theta), cos(theta), 0;...
        0, 0, 1].';
    radius_pf = (C_PF_polar_obj * radius_polar); 
    radius = (C_N_PF_obj * radius_pf).';
    R(1, :) = radius;
    
    radius_polar = [h^2/(u * (1 + e * cos(theta1)));0;0];
        C_PF_polar_obj = [cos(theta1), sin(theta1), 0;...
        -sin(theta1), cos(theta1), 0;...
        0, 0, 1].';
    radius_pf = (C_PF_polar_obj * radius_polar); 
    radius = (C_N_PF_obj * radius_pf).';
    R(2, :) = radius;
    
    radius_polar = [h^2/(u * (1 + e * cos(theta2)));0;0];
        C_PF_polar_obj = [cos(theta2), sin(theta2), 0;...
        -sin(theta2), cos(theta2), 0;...
        0, 0, 1].';
    radius_pf = (C_PF_polar_obj * radius_polar); 
    radius = (C_N_PF_obj * radius_pf).';
    R(3, :) = radius;
    
    sat_orbit = coordinate; 
    position = R;
end