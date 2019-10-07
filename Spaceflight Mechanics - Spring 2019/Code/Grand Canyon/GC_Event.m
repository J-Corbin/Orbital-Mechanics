function [value,isterminal,direction] = GC_Event(t,y)
% The second arguement of this function is the system state vector being
% integrated by ODE45
global d
%Event funtion to stop integration when rock hit -1500[m]
%if y(3) > -1500 %If rock is higher than -1500[m]
if y(3) > -d %If rock is higher than -1500[m]
 value = 1; %Keep going
else %If not
 value = 0; %Then stop
end
isterminal = 1; %Terminate integration when condtion met
direction = 0; %Direction doesn't matter
end