function [mag_e, a, theta, Omega, i, omega] = orbit_elements(r2, v2)
    global u

    mag_r2 = norm(r2);
    mag_v2 = norm(v2);

    %Momentum
    h = cross(r2, v2);
    mag_h = norm(h);

    %Eccentricity
    e = (cross(v2, h) / u) - (r2 / mag_r2);
    mag_e = norm(e);

    %Perifocal calculations
    w = h / mag_h;
    p = e / mag_e;
    q = cross(w, p);

    %Direction Cosine Matric (PF to ijk)
    C_pf = [p; q; w];

    %radius in PF coords
    %r2_pf = C_pf * r2.';
    %r2_pf = r2_pf.';

    %Initial actual anomally 
    %theta0 = atan2(r2_pf(2), r2_pf(1)) * 180 / pi;

    %Energy
    E = (mag_v2^2 / 2) - (u / mag_r2);

    %Semi-Major Axis
    a = -(u / (2 * E));

    %Perigee and Apogee Radii
    %r_perigee = a * (1 - mag_e) * p;
    %r_p = norm(r_perigee);
    %r_apogee = a * (1 + mag_e) * (-p);
    %r_a = norm(r_apogee);
    
    %Calculating inclination
    i = acos(h(3)/mag_h) * 180 / pi;

    %Right ascension of the ascending node
    N = cross([0; 0; 1], h);
    mag_N = norm(N);
    if N(2) > 0
        Omega = acos(N(1)/mag_N) * 180 / pi;
    else
        Omega = 360 - acos(N(1)/mag_N) * 180 / pi;
    end

    %Argument of Perigee
    if e(3) > 0
        omega = acos(dot(N/mag_N, e/mag_e)) * 180 / pi;
    else
        omega = 360 - acos(dot(N/mag_N, e/mag_e)) * 180 / pi;
    end

    %True anomaly
    if dot(r2, v2) >= 0
        theta = acos(dot(e/mag_e, r2/mag_r2)) * 180 / pi;
    else
        theta = 360 - acos(dot(e/mag_e, r2/mag_r2)) * 180 / pi;
    end
end