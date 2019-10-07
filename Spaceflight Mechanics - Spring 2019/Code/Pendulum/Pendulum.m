%% Pendulum Coding Problem 
% Johnathan Corbin - Spaceflight Mechanics - Spring 2019
%%
clear, clc

global L m g rho drag_coefficient interval area

%Declaring constants
L = 2; %length of rod [meters]
m = .04; %mass of sphere [kg]
radius = .05; %radius of sphere [m]
g = 9.81; %acceleration due to gravity [m/s^2]
rho = 1.225; %density of air [kg/m^3]
drag_coefficient = .42; %coefficient of drag
interval = [0,10]; %time interval [s]

area = pi * radius^2; %frontal area of sphere [m^2]

%% 75 Degree Displacement - Numerical with Drag
initial_conditions = [75 * pi / 180; 0]; %initial conditions [deg deg/s]

fname = 'PendulumFuncDrag'; %Setting name of the .m file for the state 
%variable time derivatives 

[time, y] = ode45(fname, interval, initial_conditions);

figure(1)
plot(time, y(:,1) * 180 / pi, '--b');
ylabel('Displacement (deg)')
xlabel('Time (sec)')
hold on
grid on

%% 75 Degree Displacement - Numerical without Drag
fname = 'PendulumFunc';

[time, y] = ode45(fname, interval, initial_conditions);

plot(time, y(:,1) * 180 / pi, 'k');

%% 75 Degree Displacement - Analytical 
time = 0:.01:10; %time interval (sec)
theta0 = 75 * pi / 180;

theta = theta0 * cos(sqrt(g / L) * time);

plot(time, theta * 180 / pi, '-.r')
legend('Numerical with Drag', 'Numerical without Drag', 'Analytical',...
    'Location', 'bestoutside')
title('75 Degree Initial Displacement')

%% 170 Degree Displacement - Numerical with Drag 
initial_conditions = [170 * pi / 180; 0]; %initial conditions [deg deg/s]

fname = 'PendulumFuncDrag'; %Setting name of the .m file for the state 
%variable time derivatives 

[time, y] = ode45(fname, interval, initial_conditions);

figure(2)
plot(time, y(:,1) * 180 / pi, '--b');
ylabel('Displacement (deg)')
xlabel('Time (sec)')
hold on
grid on

%% 170 Degree Displacement - Numerical without Drag
fname = 'PendulumFunc';

[time, y] = ode45(fname, interval, initial_conditions);

plot(time, y(:,1) * 180 / pi, 'k');

%% 170 Degree Displacement - Analytical
time = 0:.01:10; %time interval (sec)
theta0 = 170 * pi / 180;

theta = theta0 * cos(sqrt(g / L) * time);

plot(time, theta * 180 / pi, '-.r')
legend('Numerical with Drag', 'Numerical without Drag', 'Analytical',...
    'Location', 'bestoutside')
title('170 Degree Initial Displacement')