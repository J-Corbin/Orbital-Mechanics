function dydt = Canyon(t, y)

global A CD rOP m rho Omega d u 

dydt = zeros(size(y));

F_g = -(u*m/(rOP + y(3))^2)*[0; 0; 1];

F_D = -(1/2) * rho * CD * A * sqrt(y(4)^2 + y(5)^2 + y(6)^2) * [y(4); y(5); y(6)];

rOQ = (rOP) * [0; 0; 1] + [y(1); y(2); y(3)];

acceleration = (F_g + F_D - m * (cross(Omega, cross(Omega, rOQ)) + 2*cross(Omega,[y(4);y(5);y(6)])))/m;

dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);
dydt(4) = acceleration(1);
dydt(5) = acceleration(2);
dydt(6) = acceleration(3);