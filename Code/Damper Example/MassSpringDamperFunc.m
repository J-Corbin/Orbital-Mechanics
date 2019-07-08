 function dydt =MassSpringDamperFunc(t, y)
% Program determine the first order time derivatives of variable vector y = [y1,y2] for
% the current values of y and time t.
% IMPORTANT NOTE!!!
% This script will only work if the parameters m, k1, k2, D1, and D2 have been declared
% as glogal variables within the calling script
global m k1 k2 D1 D2 % declare what quantities are global so that they
 % can be easily used by the the dydt script
dydt = zeros(size(y)); % initialize to zero and makes column vector
dydt(1)=y(2); % Equation (19) from in handout
dydt(2) = -((D1*y(2)+D2*sign(y(2))*y(2)^2)+(k1*y(1)+k2*y(1)^3))/m; % Equation (20) in handout
%EOF MassSpringDamperDyDt