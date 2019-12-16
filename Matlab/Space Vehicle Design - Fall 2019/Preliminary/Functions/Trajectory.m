function [r, v, c3, v_inf_arrival, tof] = Trajectory(porkchop, launch_day, arrival_day, dark)
cspice_furnsh('naif0012.tls');
cspice_furnsh('de421.bsp');
cspice_furnsh('pck00010.tpc');

%% Generate Pork Chop Plots
% Generates porkchop plot if user requests
if porkchop == 'y'
    Porkchop(launch_day, arrival_day)
else
    fprintf('\n No porkchop plots will be displayed.\n')
    fprintf('------------------------------------------------------------------')
end
%% Constants Declaration
% Gravitational parameters
mu_sun = 1.32712e11;

% Values for plots
r_earth = .5e7;
r_sun = 1.5e7; 
r_mars = r_earth;
yellow = [1 1 .1];
blue = [0 .5 1];
red = [.8 0 0];
if dark == 'y'
    intercept_color = 'w'
else
    intercept_color = 'k'
end

%% State Vector Calculation
depart_time = cspice_str2et(launch_day);
arrival_time = cspice_str2et(arrival_day);
tof = arrival_time - depart_time;

% Produce Earth state vectors for departure time
[e_state_depart, ~] = cspice_spkezr('EARTH', depart_time, 'J2000', 'LT', 'SUN');
% Produce Mars state vectors for departure time
[m_state_depart, ~] = cspice_spkezr('MARS', depart_time, 'J2000', 'LT', 'SUN');

% Produce Earth state vectors for arrival time
[e_state_arrival, ~] = cspice_spkezr('EARTH', arrival_time, 'J2000', 'LT', 'SUN');
% Produce Mars state vectors for arrival time
[m_state_arrival, ~] = cspice_spkezr('MARS', arrival_time, 'J2000', 'LT', 'SUN');

%% Solving Lambert's Problem for Intercept Trajectory
%[v1, v2] = lambert(e_state_depart(1:3), m_state_arrival(1:3), tof, mu_sun, 'pro');
[v1, v2] = glambert(mu_sun, [e_state_depart(1:3);e_state_depart(4:6)],...
    [m_state_arrival(1:3);m_state_arrival(4:6)], tof, 0);

dv_depart = v1 - e_state_depart(4:6);
c3 = dv_depart' * dv_depart;

dv_arrival = v2 - m_state_arrival(4:6);
v_inf_arrival = norm(dv_arrival);

r = m_state_arrival(1:3);
v = v2;

%% Plotting the Intercept Trajectory
figure
l_angle = 0;
fudge_factor = .06;
u_angle = fudge_factor + l_angle + acos(dot(e_state_depart(1:3), m_state_arrival(1:3))...
    / (norm(m_state_arrival(1:3) * norm(e_state_depart(1:3)))));
plot_orbit(m_state_arrival(1:3), v2, mu_sun, l_angle, u_angle, intercept_color)

%% Plotting Earth's Orbit
l_angle = 0;
u_angle = 2 * pi;
plot_orbit(e_state_depart(1:3), e_state_depart(4:6), mu_sun, l_angle, u_angle, '--b')

%% Plotting Mar's Orbit
l_angle = 0;
u_angle = 2 * pi;
plot_orbit(m_state_depart(1:3), m_state_depart(4:6), mu_sun, l_angle, u_angle, '--r')

%% Plotting Various Items
% Plotting vector to positions and text
plot3([0, e_state_depart(1)], [0, e_state_depart(2)], [0, e_state_depart(3)], '-.b')

plot3([0, m_state_depart(1)], [0, m_state_depart(2)], [0, m_state_depart(3)], '-.r')

plot3([0, e_state_arrival(1)], [0, e_state_arrival(2)], [0, e_state_arrival(3)], 'b')

plot3([0, m_state_arrival(1)], [0, m_state_arrival(2)], [0, m_state_arrival(3)], 'r')


% Plot the sun
plot_planet(r_sun, yellow, [0, 0, 0])

% Plot the earth at departure
plot_planet(r_earth, blue, e_state_depart(1:3))

% Plot the earth at arrival
plot_planet(r_earth, blue, e_state_arrival(1:3))

% Plot Mars at departure
plot_planet(r_mars, red, m_state_depart(1:3))

% Plot Mars at arrival
plot_planet(r_mars, red, m_state_arrival(1:3))

% Random graph adjustments
grid on
axis equal
xlabel('km')
ylabel('km')
zlabel('km')
view(3)
title('Orbits and Intercept Trajectory to Mars')
legend('Intercept Trajectory', 'Earth Orbit', 'Mars Orbit')

fprintf('\n The C3 value for this launch is: %g km^2/s^2', c3)
fprintf('\n The v_inf value at Mars is: %g km/s', v_inf_arrival)
fprintf('\n')

% Changing plot colors to DARKMODE if desired
if dark == 'y'
    set(gca, 'color', [0 0 0])
    set(gca, 'xcolor', 'W')
    set(gca, 'ycolor', 'W')
    set(gca, 'zcolor', 'W')
    text(e_state_depart(1)*1.2, e_state_depart(2)*1.2, e_state_depart(3)*1.2, 'Departure', 'Color', 'W')
    text(m_state_depart(1)*1.2, m_state_depart(2)*1.2, m_state_depart(3)*1.2, 'Departure', 'Color', 'W')
    text(e_state_arrival(1)*1.2, e_state_arrival(2)*1.2, e_state_arrival(3)*1.2, 'Arrival', 'Color', 'W')
    text(m_state_arrival(1)*1.2, m_state_arrival(2)*1.2, m_state_arrival(3)*1.2, 'Arrival', 'Color', 'W')
else
    text(e_state_depart(1)*1.2, e_state_depart(2)*1.2, e_state_depart(3)*1.2, 'Departure')
    text(m_state_depart(1)*1.2, m_state_depart(2)*1.2, m_state_depart(3)*1.2, 'Departure')
    text(e_state_arrival(1)*1.2, e_state_arrival(2)*1.2, e_state_arrival(3)*1.2, 'Arrival')
    text(m_state_arrival(1)*1.2, m_state_arrival(2)*1.2, m_state_arrival(3)*1.2, 'Arrival')
end