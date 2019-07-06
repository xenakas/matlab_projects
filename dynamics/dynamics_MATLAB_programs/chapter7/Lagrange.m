function [L] = Lagrange(E,Q,q,t)

% dT/d(dq)
Tdq = deriv(E, diff(q,t));

% d(dT/d(dq))/dt
Tt = diff(Tdq, t);

% dT/dq
Tq = deriv(E, q);

% left hand side of Lagrange's eom
LHS = Tt - Tq;
 
% Lagrange's equation of motion
L = LHS-Q;
