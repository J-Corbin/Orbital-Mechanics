function [E0, C, mag_e, R] = orbit_parameters(r, v, t)
    global u
    
    h = cross(r, v);
    e = cross(v, h)/u - r/norm(r);
    mag_e = norm(e);
    
    term1 = sqrt((1 + norm(e))/(1 - norm(e)));
    
    p = e/norm(e);
    w = h/norm(h);
    q = cross(w, p);
    
    C = [p; q; w];
    
    r_pf = C * r.';
    r_pf = r_pf.';
    theta0 = atan2(r_pf(2), r_pf(1));
    
    term2 = tan(theta0 / 2);
    
    E0 = 2 * atan2(term2, term1);
    
    r_perigee = norm(h)^2 / (u * (1 + mag_e));
    
    a = r_perigee / (1 - mag_e);
    
    t0 = (E0 - mag_e*sin(E0))/sqrt(u/a^3);
    
    t1 = t0 + t(2);
    t2 = t1 + t(3);
    
    M1 = t1 * sqrt(u / a^3);
    M2 = t2 * sqrt(u / a^3);
    
    E1 = kepler_E(mag_e, M1);
    E2 = kepler_E(mag_e, M2);
    
    theta1 = 2* atan2((sqrt((1 + mag_e)/(1 - mag_e))) * tan(E1 / 2), 1);
    theta2 = 2* atan2((sqrt((1 + mag_e)/(1 - mag_e))) * tan(E2 / 2), 1);
    
    r1 = [norm(h)^2 / (u * (1 + mag_e * cos(theta1))), 0, 0];
    r2 = [norm(h)^2 / (u * (1 + mag_e * cos(theta2))), 0, 0];
    
    C1 = [cos(theta1), -sin(theta1), 0;...
        sin(theta1), cos(theta1), 0;...
        0, 0, 1]; %Polar to perifocal for theta1
    C2 = [cos(theta2), -sin(theta2), 0;...
        sin(theta2), cos(theta2), 0;...
        0, 0, 1]; %Polar to perifocal for theta2
    
    r1_pf = (C1 * r1.');
    r2_pf = (C2 * r2.');
    
    R1 = (C * r1_pf);
    R2 = (C * r2_pf);
    
    R = [r; R1.'; R2.'];
end