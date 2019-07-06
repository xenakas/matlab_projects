% example 4.4

clear all; clc; close 

syms m g mu alpha h x y vB d real

ls = {alpha, h, mu, g};
ln = {pi/6, 2, 0.2, 9.8};

% (a)
G = [0 -m*g 0];

Ff = mu*...
[-m*g*(cos(alpha))^2 m*g*cos(alpha)*sin(alpha) 0];

UAB = int(G(1)+Ff(1),x,0,h/tan(alpha))...
    +int(G(2)+Ff(2),y,h,0);

U_AB = simplify(UAB);

fprintf('work done by the particle \n')
fprintf('U_AB = %s  \n',char(U_AB))
fprintf('\n\n')

TB = (m*vB^2)/2;
TA = 0;
DeltaTAB = TB-TA;
fprintf('Delta(T_AB) = %s \n',char(DeltaTAB));
fprintf('\n\n');

vBs = solve(U_AB-DeltaTAB,vB);
fprintf('vB =  \n')
pretty(vBs(2))
fprintf('\n\n')

vBn = subs(vBs(2), ls, ln);
fprintf('vB = %6.3f (m/s)\n\n',vBn)

% (b)
Ffx=subs(Ff(1),alpha, pi);
Ffy=subs(Ff(2),alpha, pi);

U_BC=int(G(1)+Ffx,x,h/tan(alpha),h/tan(alpha)+d);
fprintf('U_BC = %s  \n',char(U_BC));
fprintf('\n\n');

TC = 0;
DeltaTBC=subs(TC-TB,'vB',vBs(2));
fprintf('Delta T_BC = %s  \n',char(DeltaTBC))
fprintf('\n\n')

d=solve(U_BC-DeltaTBC,'d');
fprintf('d = %s  \n',char(d))
fprintf('\n\n')

dn = subs(d, ls, ln);
fprintf('d = %6.3f (m)\n', dn)

% end of program