function [r, v] = improvement(r2, v2, P1, P2, P3, tau, f1, f3, g1, g3)
global u
P1_old = P1;  P2_old = P2;  P3_old = P3;
diff1    = 1;     diff2    = 1;     diff3    = 1;
n    = 0;
nmax = 1000;
tol  = 1.e-8;

%...Iterative improvement loop:
while ((diff1 > tol) & (diff2 > tol) & (diff3 > tol)) & (n < nmax)
    n = n+1;

%...Compute quantities required by universal kepler's equation:
    ro  = norm(r2);
    vo  = norm(v2);
    vro = dot(v2,r2)/ro;
    alpha   = 2/ro - vo^2/u;

%...Solve universal Kepler's equation at times tau1 and tau3 for
%   universal anomalies x1 and x3:
    x1 = kepler_U(tau(1), ro, vro, alpha);
    x3 = kepler_U(tau(3), ro, vro, alpha);

%...Calculate the Lagrange f and g coefficients at times tau1
%   and tau3:
    [ff1, gg1] = f_and_g(x1, tau(1), ro, alpha);
    [ff3, gg3] = f_and_g(x3, tau(3), ro, alpha);

%...Update the f and g functions at times tau1 and tau3 by 
%   averaging old and new:
    f1    = (f1 + ff1)/2;
    f3    = (f3 + ff3)/2;
    g1    = (g1 + gg1)/2;
    g3    = (g3 + gg3)/2;

%...Equations 5.96 and 5.97:
    c1    =  g3/(f1*g3 - f3*g1);
    c3    = -g1/(f1*g3 - f3*g1);

%...Equations 5.109a, 5.110a and 5.111a:
    P1  = 1/D0*(      -D(1,1) + 1/c1*D(2,1) - c3/c1*D(3,1));
    P2  = 1/D0*(   -c1*D(1,2) +      D(2,2) -    c3*D(3,2));
    P3  = 1/D0*(-c1/c3*D(1,3) + 1/c3*D(2,3) -       D(3,3));

%...Equations 5.86:
    r1 = R(1,:) + P1*p(1,:);
    r2 = R(2,:) + P2*p(2,:);
    r3 = R(3,:) + P3*p(3,:);

%...Equation 5.118:
    v2    = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);

%...Calculate differences upon which to base convergence:
    diff1 = abs(P1 - P1_old);
    diff2 = abs(P2 - P2_old);
    diff3 = abs(P3 - P3_old);

%...Update the slant ranges:
    P1_old = P1;  P2_old = P2;  P3_old = P3;
end
%...End iterative improvement loop

fprintf('\n( **Number of Gauss improvement iterations = %g)\n\n',n)

if n >= nmax
	fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax);
end
end