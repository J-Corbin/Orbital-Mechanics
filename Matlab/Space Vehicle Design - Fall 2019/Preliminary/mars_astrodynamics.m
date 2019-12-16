clear, clc, close all

%% Variable Declaration
% Gravitational Parameters
mu_mars = 42828;
mu_sun = 1.32712e11;

% Unit conversion
deg_rad = pi / 180; %converting degree to rad
day2sec = 24*3600;
sec2day = 1/day2sec;

%% Intercept Trajectory
% Plot the orbits of Earth and Mars, along with the intercept trajectory
%Planned dates: launch: Sep 10, 2022
%               arrive: Mar 28, 2023
launch_day = 'Sep 10, 2022 12:00:00.000000';
arrival_day = 'Mar 28, 2023 12:00:00.000000';
[r, v, c3, v_inf, tof] = Trajectory('n', launch_day, arrival_day, 'n');
tof = tof * sec2day;
fprintf('------------------------------------------------------------------')
fprintf('\n The transfer to Mars will take approximately %g days.\n', tof)

% Calculate the orbital elements of the intercept trajectory
[mag_e, a, ~, Omega, i, omega] = orbit_elements(r, v, mu_sun);
fprintf('------------------------------------------------------------------')
fprintf('\n The orbital elements of the transfer trajectory are:')
fprintf('\n Eccentricity =                                       %g', mag_e)
fprintf('\n Semi-Major Axis [km] =                              %g', a)
fprintf('\n Right of Ascension of the Ascending Node [degrees] = %g', Omega)
fprintf('\n Inclination [degrees] =                              %g', i)
fprintf('\n Argument of Perogee [degrees] =                      %g', omega)
fprintf('\n')
fprintf('------------------------------------------------------------------')
fprintf('\n')


%% Orbital Elements
% Calculate the orbital elements 
r_perogee = 3396; %Perogee = radius of Mars for landing
e = 1 + r_perogee * v_inf^2 / (mu_mars^2);
h = r_perogee * sqrt(v_inf^2 + 2 * mu_mars / r_perogee);