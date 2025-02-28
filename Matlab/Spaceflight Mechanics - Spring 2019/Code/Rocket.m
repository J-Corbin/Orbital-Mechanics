%% Rocket Problem
% Spaceflight Mechanics - Spring 2019
% Johnathan Corbin

function Rocket

%% Clearing console and variables.
clear, clc

%% Declaring constant values given from problem statement.
Isp = 353; %LOx - Kerosene propellant, [s]
epsilon = 0.01; %Rocket's structural ratio
mass_L = 1000; %Payload mass, [kg]
mass_P = 10000; %Propellant mass, [kg]
diameter = 1.65; %Rocket diameter, [m]
t_burn = 240; %Total burn time, [s]
r_earth = 6378 * 1000; %Radius of the Earth, [m]
t0 = 0; %Initial time, [s]
pitch = 500; %Height at which pitchover begins, [m]
u = 398600 * 1000^3; %Standard gravitational parameter for Earth, [m^3/s^2]
s_g = 9.81; %Standard gravity, [m/s^2]
rho0 = 1.225; %Sea level air density, [kg/m^3]
hscale = 7500; %Density scale height, [m]

%% Intermediate calculations of various values.
A = pi * (diameter / 2)^2; %Frontal area of the rocket, [m^2]
tspan = [t0, t_burn]; %Time range for the integrator, [s]
mass_structure = epsilon * mass_P / (1 - epsilon); %Structural mass of the rocket, [kg]
m_dot = mass_P / t_burn; %Mass rate of ejecta, [kg/s]
thrust = m_dot * Isp * s_g; %Thrust generated by ejecta, [N]

%% No atmosphere, no pitchover angle.
psi0 = 0; %Pitchover angle, [rad]
Cd = 0; %Drag coefficient
v0 = 0; %Initial velocity, [m/s]
x0 = 0; %Initial downrange distance, [m]
h0 = 0; %Initial height, [m]
theta0 = 0; %Initial angular displacement, [rad]
phi0 = 0; %Initial angle between launch site and tangential direction, [rad]

f0 = [phi0, v0, h0, x0, theta0, psi0]; % Initial conditions vector

%Integrating with no pitch until pitch height
Opt1 = odeset('Events', @Begin_Pitch);
[time, f] = ode45(@rocketMan, tspan, f0);
v = f(:, 2); %Velocity, [m/s]
h = f(:, 3); %Height, [m]
x = f(:, 4); %Downrange, [m]

figure(1)
subplot(2, 2, [1 2])
plot(x / 1000, h / 1000, '--k')
hold on

%{
%New initial conditions for second integration
phi0 = phi(end);
v0 = v(end); %Initial velocity, [m/s]
x0 = x(end); %Initial downrange distance, [m]
h0 = h(end); %Initial height, [m]
theta0 = theta(end);
psi0 = 0;
tspan = [time(end), t_burn];
f0 = [phi0, v0, h0, x0, theta0, psi0];

%Secondary integration with pitchover
Opt2 = odeset('Events', @Crashed);
[~, f] = ode45(@rocketMan, tspan, f0, Opt2);
v = f(:, 2); %Velocity, [m/s]
h = f(:, 3); %Height, [m]
x = f(:, 4); %Downrange, [m]

plot(x / 1000, h / 1000, '--b')
grid on
title('No Atmosphere, No Pitch')
xlabel('Downrange [km]')
ylabel('Height [km]')
%}

%% Functions
function dydt = rocketMan(t, y)
    dydt = zeros(size(y));
    
    phi = y(1);
    v = y(2);
    h = y(3);
    x = y(4);
    theta = y(5);
    psi = y(6);
    
    if t < t_burn
       mass = mass_structure + mass_L + mass_P - m_dot * t;
       T = thrust;
    else
        T = 0;
        mass = mass_structure + mass_L;
    end
    
    %Intermediate calculations
    gravity = -(mass * u) / (r_earth + h)^2;
    density = rho0 * exp(-h / hscale);
    drag = -0.5 * density * Cd * A * v^2;
    
    %Differential equations
    phi_dot = (u / (r_earth + h)^2) * sin(psi) / v;
    v_dot = (T + drag + gravity * cos(psi)) / mass;
    h_dot = v * cos(psi);
    x_dot = v * sin(psi);
    theta_dot = (v * sin(psi)) / (r_earth + h);
    psi_dot = phi_dot - theta_dot;

    dydt(1) = phi_dot;
    dydt(2) = v_dot;
    dydt(3) = h_dot;
    dydt(4) = x_dot;
    dydt(5) = theta_dot;
    dydt(6) = psi_dot;
end

function [value,isterminal,direction] = Begin_Pitch(~,y)
%Event funtion to stop integration when rocket reaches 500[m]
if y(4) < pitch
 value = 1; %Keep going
else %If not
 value = 0; %Then stop
end
isterminal = 1; %Terminate integration when condtion met
direction = 0; %Direction doesn't matter
end

function [value,isterminal,direction] = Crashed(~,y)
%Event funtion to stop integration when rocket reaches 500[m]
if y(4) > 0
 value = 1; %Keep going
else %If not
 value = 0; %Then stop
end
isterminal = 1; %Terminate integration when condtion met
direction = 0; %Direction doesn't matter
end

end