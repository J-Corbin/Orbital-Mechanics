% This script solves the equations of motion associated with a Maass-spring-damper system
% discussed in an associated handout (which include a cubic stiffness term and hydrodynamic
% damping.

% The problem is solve the equations associated with the linearized model, using the closed
% closed form solution for that idealized model. Next the full nonlinear equations are integrated
% numerical using "ode45"
clear, clc
clear global variable %clear the global variables so that you are sure of what your are starting with
global m k1 k2 D1 D2 % declare what quantities are global so that they
 % can be easily used by the the dydt script

% Specify System parameters
m = 1; %[kg]
k1 = 100; %[n/m]
k2 = 1.0; %[N/m^3]
D1 = 0.5; %[(N*s)/m]
D2 = 0.01; % %[(N*s^2)/m^2]

% Specify Inital Conditions
delta0 = 10; % [m]
ddelta0 = 0; % [m/s]
% Specify duration of simulation
Tfinal = 5.0; % [s]

%% Numerical Solution
tspan = [0,Tfinal]; %specify the time interval over which the simulation is
 % run
y0 = [delta0; ddelta0]; % Specify intial value for first order state variables
fname = 'MassSpringDamperFunc'; % specify name of .m file containing the
 % script which determines the state variable time derivatives

% Run ODE45 to integrate these equations and return the result
 [t2,y] = ode45(fname,tspan,y0);

%% Closed Form Solution (for linearized System)
t = linspace(0,Tfinal,length(y));

omega = sqrt(k1/m); % equation (14) from Handout
zeta = D1/(2*sqrt(k1*m)); % equation (13) from handout

A=delta0;
B= (ddelta0+zeta*omega*A)/(omega*sqrt(1-zeta^2)); % eqns (16) from handout

% Determine delta as a function of t
omegad = omega*sqrt(1-zeta^2); % Damped natural circular frequency
delta = exp(-zeta*omega*t).*(A*cos(omegad*t) + B*sin(omegad*t));

% Plot results for both the Anlytic and Numeric Solutions on single plot
 plot(t,delta, 'r', t2,y(:,1), 'b')
 legend('Analytic Solution','Numerical Solution')
 title ('Comparison of Analytic and Numerical Results')
 xlabel('Time [sec]')
 ylabel('Dispacement [m]')