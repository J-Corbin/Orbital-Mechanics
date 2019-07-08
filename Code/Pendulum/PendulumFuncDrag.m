 function dydt = PendulumFuncDrag(t, y)
% Program determine the first order time derivatives of variable vector y = [y1,y2] for
% the current values of y and time t.
global L m g rho drag_coefficient area
 % can be easily used by the the dydt script

 dydt = zeros(size(y)); % initialize to zero and makes column vector

 dydt(1) = y(2);

 dydt(2) = (-sign(y(2)) * rho * L * drag_coefficient * y(2)^2 * area) / (2 * m) - (g * sin(y(1))) / L;
