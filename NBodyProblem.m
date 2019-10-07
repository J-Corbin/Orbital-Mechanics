%% N-Body Problem
% Johnathan Corbin

%% Clearing console and variables.
clear, clc

%% Declaring parameters of every body
%body number = [x, y, z, Vx, Vy, Vz, mass]
body1 = [-10, 10, 0, 0, 0, 0, 10];
body2 = [10, 10, 0, 0, 0, 0, 10];
body3 = [10, -10, 0, 0, 0, 0, 10];
body4 = [-10, -10, 0, 0, 0, 0, 10];

bodyMatrix = [body1, body2, body3, body4];

global G numBodies
G = 6.6743 * 10^-11;
numBodies = 4;

% Time to run the equation solver.
tspan = [0, 1000];

%% Intermediate calculations
r = zeros([numBodies-1 1]);

for index = 1:(numBodies - 1)
        r(index) = sqrt((bodyMatrix(1 + (7 * (numBodies - 1))) - bodyMatrix(1 + (7 * (index - 1))))^2 + ...
        (bodyMatrix(2 + (7 * (numBodies - 1))) - bodyMatrix(2 + (7 * (index - 1))))^2 + ...
        (bodyMatrix(3 + (7 * (numBodies - 1))) - bodyMatrix(3 + (7 * (index - 1))))^2);
end

%% Solving equations of motion
bodyMatrix = [body1, body2, body3]; %Initial conditions vector
[time, results] = ode45(@nBody, tspan, bodyMatrix);

%% Rates function called by ODE45.
function dydt = nBody(t, params)
    global G numBodies
    
    %dydt = zeros([size(params)-numBodies 1]);
    dydt = zeros([6 1]);
    
    for index = 1:numBodies
        vx_dot = vx_dot - (G*params(7*index)*(params(1 + (7 * (numBodies - 1))) - params(1 + (7 * (index - 1))) / r(index));
        vy_dot = vy_dot - (G*params(7*index)*(params(2 + (7 * (numBodies - 1))) - params(2 + (7 * (index - 1)))) / r(index));
        vz_dot = vz_dot - (G*params(7*index)*(params(3 + (7 * (numBodies - 1))) - params(3 + (7 * (index - 1)))) / r(index));
    end
    
    
    start = 4;
    for index = 1:numBodies
        for velocity = 1:9
           dydt(velocity) = params(start + (7 * (index - 1)));
           start = start + 1;
        end
        for acceleration = 10:18
           dydt(acceleration) =  
        end
        if start == 6
           start = 4; 
        end
    end
    
    
    %{
    dydt = zeros(size(y));
    
    v = y(1);
    theta = y(2);
    x = y(3);
    h = y(4);
    
    theta_dot = 0;
    v_dot = (T + drag + g) / mass;
    x_dot = 0;
    h_dot = v;
    
    dydt(1) = v_dot;
    dydt(2) = theta_dot;
    dydt(3) = x_dot;
    dydt(4) = h_dot;
    %}
end
