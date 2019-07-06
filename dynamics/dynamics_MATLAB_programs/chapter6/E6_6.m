% example 6.6
clear all; clc; close all 

syms b h d e 
rO = [0,0,0];
rC = [0,b/2,d+h/2];
rP = [0,0,d+h+e];

syms t
q = sym('q(t)');
theta = [0,0,q];
omega = diff(theta,t);
alpha = diff(omega,t);

syms FOx FOy FOz FPx FPy FPz
% reaction force at O
FO = [FOx,FOy,FOz];
% reaction force at P
FP = [FPx,FPy,0];

% gravitational force at C
syms m g mE
G = [0,0,-m*g];

ICyy=m*h^2/12;
ICzz=m*b^2/12;
ICxx=ICyy+ICzz;
ICyz=0;
ICxy=0;
ICxz=0;

Iyy=ICyy+m*(h/2+d)^2;
Iyy=simplify(expand(Iyy));
Izz=ICzz++m*(b/2)^2;
Izz=simplify(expand(Izz));
Ixx=ICxx+m*((d+h/2)^2+(b/2)^2);
Ixx=simplify(expand(Ixx));
Iyz=ICyz+m*(d+h/2)*(b/2);
Iyz=simplify(expand(Iyz));
Ixy=0;
Ixz=0;

syms omega0
list= {b, h, d, e, m, g  , omega0, mE};
listn={1, 2, 1, 1, 1, 9.8, 2, 1}; 

Ixxn=subs(Ixx, list, listn);
Iyyn=subs(Iyy, list, listn);
Izzn=subs(Izz, list, listn);
Iyzn=subs(Iyz, list, listn);

fprintf('Ixx=%s=%g (kg m^2) \n',char(Ixx),Ixxn)
fprintf('Iyy=%s=%g (kg m^2) \n',char(Iyy),Iyyn)
fprintf('Izz=%s=%g (kg m^2) \n',char(Izz),Izzn)
fprintf('Iyz=%s=%g (kg m^2) \n',char(Iyz),Iyzn)

I=[
     Ixx, -Ixy, -Ixz;
    -Ixy,  Iyy, -Iyz;
    -Ixz, -Iyz,  Izz
  ];
pretty(I)
In=subs(I, list, listn)

aO=[0,0,0];
aC=aO+cross(alpha,rC)+cross(omega,cross(omega,rC));

% Newton eom
% m aC = G+FO+FP 
eqF = m*aC-(G+FO+FP);
% Euler eom 
eqM = -alpha*I-cross(omega,omega*I)-...
     (cross(rC,G)+cross(rP,FP));

% => 
eqFx = eqF(1); % (1)
eqFy = eqF(2); % (2)
eqFz = eqF(3); % (3) 
eqMx = eqM(1); % (4)
eqMy = eqM(2); % (5)
eqMz = eqM(3); % (6) =>  

pretty(eqFx)
pretty(eqFy)
pretty(eqFz)

pretty(eqMx)
pretty(eqMy)
pretty(eqMz)

% from last equation
qlist={diff('q(t)',t,2),diff('q(t)',t)};
qlistn={0,'omega0'};

eq1 = subs(eqFx,qlist,qlistn);
eq2 = subs(eqFy,qlist,qlistn);
eq3 = subs(eqFz,qlist,qlistn);
eq4 = subs(eqMx,qlist,qlistn);
eq5 = subs(eqMy,qlist,qlistn);
sol=...
solve(eq1,eq2,eq3,eq4,eq5,...
      'FOx,FOy,FOz,FPx,FPy');
FOxs=sol.FOx;
FOys=sol.FOy;
FOzs=sol.FOz;
FPxs=sol.FPx;
FPys=sol.FPy;

FOxn=subs(FOxs, list, listn);
FOyn=subs(FOys, list, listn);
FOzn=subs(FOzs, list, listn);
FPxn=subs(FPxs, list, listn);
FPyn=subs(FPys, list, listn);

fprintf('\n\n')
fprintf('FOx = \n')
pretty(simplify(FOxs))
fprintf('\n')
fprintf('FOy = \n')
pretty(simplify(FOys))
fprintf...
('FOy = %g (N)\n',FOyn)
fprintf('\n')
fprintf('FOz = \n')
pretty(simplify(FOzs))
fprintf...
('FOz = %g (N)\n',FOzn)
fprintf('\n')
fprintf('FPx = \n')
pretty(simplify(FPxs))
fprintf('\n')
fprintf('FPy = \n')
pretty(simplify(FPys))
fprintf...
('FPy = %g (N)\n',FPyn)
fprintf('\n\n')


syms xE yE zE 
% position vector for the 
% point mass mE
rE = [xE,yE,zE];

% COM of the system has to be on 
% z-axis
% =>
% m xC + mE xE = 0
% m yC + mE yE = 0

eqEx = m*rC(1) + mE*xE;
eqEy = m*rC(2) + mE*yE;

pretty(eqEx)
pretty(eqEy)

% z-axis has to be principal axis
% mE*yE*zE+Iyz = 0 
eqEm = mE*yE*zE+Iyz;

solE = solve(eqEx,eqEy,eqEm,'xE,yE,zE');

xEs=solE.xE;
yEs=solE.yE;
zEs=solE.zE;

xEn=subs(xEs, list, listn);
yEn=subs(yEs, list, listn);
zEn=subs(zEs, list, listn);

fprintf('xE=%s=%g (m) \n',char(xEs),xEn)
fprintf('yE=%s=%g (m) \n',char(yEs),yEn)
fprintf('zE=%s=%g (m) \n',char(zEs),zEn)

% end of program