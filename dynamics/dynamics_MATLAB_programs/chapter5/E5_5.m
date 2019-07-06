% example 5.5

clear all; clc; close all

% Input data (m)
AB = 0.15; AC = 0.30; CD = 0.20; 	
DE = 0.14; CF = 0.50; a  = 0.25; 
phi = 200*pi/180; % rad
% Position of joint A
xA=0; yA=0; rA=[xA yA 0]; 
% Position of joint C
xC=-AC; yC=0; rC=[xC yC 0];
% Position of joint B
xB=AB*cos(phi); yB=AB*sin(phi); rB=[xB yB 0]; 
% Link 3
phi3=atan((yB-yC)/(xB-xC)); phi2=phi3;
% Position of joint F
xF=xC+CF*cos(phi2); yF=yC+CF*sin(phi2); rF=[xF,yF,0];
% Position of joint D
xD=xC-CD*cos(phi3); yD=yC-CD*sin(phi3); rD=[xD yD 0]; 
% Position of joint E
xE=-AC-a; yE=yD+sqrt(DE^2-(xE-xD)^2); rE=[xE yE 0];
% Link 4
phi4=atan((yE-yD)/(xE-xD));

% Graphic of the mechanism
la=-xE+0.05;
plot(...
[xA,xB],[yA,yB],'r-o',...
[xD,xF],[yD,yF],'k-o',...
[xC,xC],[yC,yC],'k-o',...
[xD,xE],[yD,yE],'b-o',...
[xE,xE],[-la,la],'r--')
grid
text(xA,yA,'   A')
text(xB,yB,'   B')
text(xC,yC,'   C')
text(xD,yD,'  D')
text(xE,yE,'  E')
text(xF,yF,'  F')
xlim([-la la]);
ylim([-la la]);

fprintf('Results \n\n') 
fprintf('phi = phi1 = %g (degrees)\n',phi*180/pi) 
fprintf('rA = [%g,%g,%g] (m)\n', rA) 
fprintf('rB = [%6.3f,%6.3f,%g] (m)\n', rB) 
fprintf('rC = [%6.3f,%6.3f,%g] (m)\n', rC) 
fprintf('rD = [%6.3f,%6.3f,%g] (m)\n', rD) 
fprintf('rE = [%6.3f,%6.3f,%g] (m)\n', rE) 
fprintf('rF = [%6.3f,%6.3f,%g] (m)\n', rF) 
fprintf...
('phi2 = phi3 = %6.3f (degrees) \n',phi2*180/pi) 
fprintf...
('phi4 = %6.3f (degrees) \n',phi4*180/pi) 

fprintf('\n') 
fprintf('Velocity and acceleration analysis \n\n') 

n = 700.; % rpm 
omega1 = [ 0 0 pi*n/30 ]; 
alpha1 = [0 0 0 ]; 
fprintf...
('omega1 = [%g,%g,%6.3f] (rad/s)\n', omega1) 
fprintf...
('alpha1 = [%g,%g,%g] (rad/s^2)\n', alpha1)  
fprintf('\n') 
vA = [0 0 0]; aA = [0 0 0];
vB1 = vA + cross(omega1,rB);  
aB1 = aA + cross(alpha1,rB) - ...
    dot(omega1,omega1)*rB;
vB2 = vB1;
aB2 = aB1;
fprintf...
('vB1 = vB2 = [%6.3f,%6.3f,%g] (m/s)\n',vB1) 
fprintf('|vB1| = %6.3f (m/s)\n', norm(vB1))
fprintf...
('aB1 = aB2 = [%6.3f,%6.3f,%g] (m/s^2)\n',aB1) 
fprintf('|aB1| = %6.3f (m/s^2)\n', norm(aB1))
fprintf('\n')

aB1n = - dot(omega1,omega1)*rB;
aB1t = cross(alpha1,rB);
fprintf('aB1n = [%6.3f,%6.3f,%d] (m/s^2)\n', aB1n)
fprintf('aB1t = [%g,%d,%g] (m/s^2)\n', aB1t)
fprintf('\n');

omega3z=sym('omega3z','real'); 
vB32=sym('vB32','real');
omega3 = [ 0 0 omega3z ]; 
% omega3z unknown (to be calculated)
% vB32 unknown (to be calculated)
vC = [0 0 0 ]; % C is fixed
% vB3 = vC + omega3 x rCB   
% (B3 & C are points on link 3)
vB3 = vC + cross(omega3,rB-rC); 
% point B2 is on link 2 and point B3 is on link 3
% vB3 = vB2 + vB3B2
% between the links 2 and 3 there is a  
% translational joint B_T
% vB3B2 is the relative velocity of B3 wrt 2
% vB3B2 is parallel to the sliding direcion BC
% vB3B2 is written as a vector
vB3B2 = vB32*[cos(phi2) sin(phi2) 0]; 
% vB3 = vB2 + vB3B2
eqvB = vB3 - vB2 - vB3B2; 
% vectorial equation
% the component of the vectorial equation on x-axis
eqvBx = eqvB(1); 
% the component of the vectorial equation on y-axis
eqvBy = eqvB(2); 
% two equations eqvBx and eqvBy with two unknowns 
% solve for omega3z and vB32 
solvB = solve(eqvBx,eqvBy); 
omega3zs=eval(solvB.omega3z); 
vB32s=eval(solvB.vB32);

omega3v = [0 0 omega3zs]; 
omega2v = omega3v;
vB32v = vB32s*[cos(phi2) sin(phi2) 0];

digits 3
% print the equations for calculating 
% omega3 and vB32
fprintf...
('vB3 = vC + omega3 x rCB = vB2 + vB3B2 => \n') 
qvBx=vpa(eqvBx);
fprintf('x-axis: %s = 0 \n', char(qvBx)) 
qvBy=vpa(eqvBy);
fprintf('y-axis: %s = 0 \n', char(qvBy)) 
fprintf('=>\n') 
fprintf('omega3z = %6.3f (rad/s)\n', omega3zs) 
fprintf('vB32 = %6.3f (m/s)\n', vB32s) 
fprintf('\n') 
fprintf...
('omega2=omega3 = [%g,%g,%6.3f](rad/s)\n', omega3v) 
fprintf('vB3B2 = [%6.3f,%6.3f,%d] (m/s)\n\n', vB32v)  

% Coriolis acceleration
aB3B2cor = 2*cross(omega3v,vB32v); 

alpha3z=sym('alpha3z','real'); 
aB32=sym('aB32','real'); % aB32 unknown

alpha3 = [ 0 0 alpha3z ]; % alpha3z unknown
aC = [0 0 0 ]; % C is fixed
% aB3 acceleration of B3
aB3 = aC + cross(alpha3, rB-rC) -...
    dot(omega3v, omega3v)*(rB-rC); 
% aB3B2 relative velocity of B3 wrt 2
% aB3B2 parallel to the sliding direcion BC

aB3B2 = aB32*[ cos(phi2) sin(phi2) 0]; 
% aB3 = aB2 + aB3B2 + aB3B2cor
eqaB = aB3 - aB2 - aB3B2 - aB3B2cor; 
% vectorial equation
eqaBx = eqaB(1); % equation component on x-axis
eqaBy = eqaB(2); % equation component on y-axis
solaB = solve(eqaBx,eqaBy);
alpha3zs=eval(solaB.alpha3z);
aB32s=eval(solaB.aB32);
alpha3v = [0 0 alpha3zs];
alpha2v = alpha3v;
aB32v = aB32s*[cos(phi2) sin(phi2) 0];

% print the equations for calculating 
% alpha3 and aB32
fprintf...
('aB32cor = [%6.3f,%6.3f,%g](m/s^2)\n', aB3B2cor)  
fprintf('\n') 
fprintf('aB3=aC+alpha3xrCB-(omega3.omega3)rCB\n') 
fprintf('aB3=aB2+aB3B2+aB3B2cor =>\n') 
qaBx=vpa(eqaBx);
fprintf('x-axis:\n') 
fprintf('%s = 0\n',char(qaBx)) 
qaBy=vpa(eqaBy);
fprintf('y-axis:\n') 
fprintf('%s = 0\n',char(qaBy))  
fprintf('=>\n');
fprintf('alpha3z = %6.3f (rad/s^2)\n', alpha3zs) 
fprintf('aB32 = %6.3f (m/s^2)\n', aB32s) 
fprintf('\n') 
fprintf...
('alpha2=alpha3=[%g,%g,%6.3f](rad/s^2)\n',alpha3v) 
fprintf...
('aB3B2=[%6.3f,%6.3f,%d] (m/s^2)\n\n',aB32v) 

% vD3 velocity of D3
% D3 & C points on link 3
vD3 = vC + cross(omega3v,rD-rC);
vD4 = vD3;
fprintf...
('vD3 = vD4 = [%6.3f,%6.3f,%g] (m/s)\n', vD3) 
fprintf('|vD3| = %6.3f (m/s)\n', norm(vD3))
% aD3 acceleration of D3
aD3 = aC + cross(alpha3v, rD-rC)-...
    dot(omega3v, omega3v)*(rD-rC);  
aD4 = aD3;
fprintf...
('aD3 = aD4 = [%6.3f,%6.3f,%g] (m/s^2)\n', aD3) 
fprintf('|aD3| = %6.3f (m/s^2)\n', norm(aD3))
fprintf('\n')

vDC = cross(omega3v,rD-rC);
aDC=cross(alpha3v,rD-rC)-...
    dot(omega3v,omega3v)*(rD-rC);
aDCn=-dot(omega3v,omega3v)*(rD-rC);
aDCt=cross(alpha3v,rD-rC);

fprintf...
('vDC = [%6.3f,%6.3f,%g] (m/s)\n', vDC)
fprintf...
('aDC = [%6.3f,%6.3f,%g] (m/s^2)\n', aDC)
fprintf...
('aDCn = [%6.3f,%6.3f,%d] (m/s^2)\n', aDCn)
fprintf...
('|aDCn| = %6.3f (m/s^2)\n', norm(aDCn))
fprintf...
('aDCt = [%6.3f,%6.3f,%g] (m/s^2)\n', aDCt)
fprintf...
('|aDCt| = %6.3f (m/s^2)\n', norm(aDCt))
fprintf('\n')

% F & C points on link 3
vF = vC + cross(omega3v,rF-rC);
fprintf('vF = [%6.3f,%6.3f,%g] (m/s)\n', vF) 
% aD3 acceleration of D3
aF = aC + cross(alpha3v, rF-rC)-...
    dot(omega3v, omega3v)*(rF-rC);  
fprintf('aF = [%6.3f,%6.3f,%g] (m/s^2)\n', aF) 
fprintf('|vF| = %6.3f (m/s)\n', norm(vF))
fprintf('|aF| = %6.3f (m/s^2)\n', norm(aF))
fprintf('\n')

omega4z = sym('omega4z','real'); 
% omega4z unknown
vEy = sym('vEy','real'); 
% vEy unknown
omega4 = [ 0 0 omega4z ];
vE5 = [ 0 vEy 0];
% vEy velocity of E
vE4 = vD4 + cross(omega4,rE-rD); 
% vE parallel to the sliding direcion y
% vE5 = vE4
eqvE = vE5-vE4; 
% vectorial equation
eqvEx = eqvE(1); % component on x-axis
eqvEy = eqvE(2); % component on y-axis
solvE = solve(eqvEx,eqvEy);
omega4zs=eval(solvE.omega4z);
vEys=eval(solvE.vEy);
omega4v = [0 0 omega4zs];
vE = [0 vEys 0];

% print the equations for calculating 
% omega4 and vE
fprintf...
('vE = vD + omega4 x (rE-rD) => \n') 
qvEx=vpa(eqvEx);
fprintf('x-axis: %s = 0 \n', char(qvEx)) 
qvEy=vpa(eqvEy);
fprintf('y-axis: %s = 0 \n', char(qvEy)) 
fprintf('=>\n') 
fprintf('omega4z = %6.3f (rad/s)\n', omega4zs) 
fprintf('vEy = %6.3f (m/s)\n', vEys) 
fprintf('\n') 
fprintf...
('omega4 = [%g,%g,%6.3f] (rad/s)\n', omega4v) 
fprintf('vE = [%g,%6.3f,%g] (m/s)\n', vE) 
fprintf('|vE| = %6.3f (m/s)\n', norm(vE))

alpha4z = sym('alpha4z','real');
aEy = sym('aEy','real');
alpha4 = [ 0 0 alpha4z ]; % alpha5z unknown
% aE5 acceleration of E5
aE4=aD4+cross(alpha4,rE-rD)-...
    dot(omega4v,omega4v)*(rE-rD); 
% aE parallel to the sliding direcion y
aE5 = [0 aEy 0]; 
eqaE = aE5 - aE4; 
% vectorial equation
eqaEx = eqaE(1); % component on x-axis
eqaEy = eqaE(2); % component on y-axis
solaE = solve(eqaEx,eqaEy);
alpha4zs = eval(solaE.alpha4z);
aEys = eval(solaE.aEy);
alpha4v = [0 0 alpha4zs];
aEv = [0 aEys 0]; 

% print the equations for calculating 
% alpha4 and aE
fprintf('\n') 
fprintf('aE4=aD+alpha4xrDE-(omega4.omega4)rDE\n') 
fprintf('aE5=aE4 =>\n') 
qaEx=vpa(eqaEx);
fprintf('x-axis: %s = 0 \n', char(qaEx)) 
qaEy=vpa(eqaEy);
fprintf('y-axis: %s = 0 \n', char(qaEy)) 
fprintf('=>\n') 
fprintf('alpha4z = %6.3f (rad/s^2)\n', alpha4zs) 
fprintf('aEy = %6.3f (m/s^2)\n', aEys) 
fprintf('\n') 
fprintf...
('alpha4=[%g,%g,%6.3f] (rad/s^2)\n', alpha4v) 
fprintf('aE = [%g,%6.3f,%g] (m/s^2)\n', aEv) 
fprintf('|aE| = %6.3f (m/s^2)\n', norm(aEv))
fprintf('\n') 

vED = cross(omega4v,rE-rD);
aED=cross(alpha4v,rE-rD)-...
    dot(omega4v,omega4v)*(rE-rD);
aEDn=-dot(omega4v,omega4v)*(rE-rD);
aEDt=cross(alpha4v,rE-rD);

fprintf...
('vED = [%6.3f,%6.3f,%g] (m/s)\n', vED)
fprintf...
('aED = [%6.3f,%6.3f,%g] (m/s^2)\n', aED)
fprintf...
('aEDn = [%6.3f,%6.3f,%d] (m/s^2)\n', aEDn)
fprintf...
('|aEDn| = %6.3f (m/s^2)\n', norm(aEDn))
fprintf...
('aEDt = [%6.3f,%6.3f,%g] (m/s^2)\n', aEDt)
fprintf...
('|aEDt| = %6.3f (m/s^2)\n', norm(aEDt))
fprintf('\n')

% end of program


% 
% Results 
% 
% phi = phi1 = 200 (degrees)
% rA = [0,0,0] (m)
% rB = [-0.141,-0.051,0] (m)
% rC = [-0.300, 0.000,0] (m)
% rD = [-0.490, 0.061,0] (m)
% rE = [-0.550, 0.188,0] (m)
% rF = [ 0.176,-0.153,0] (m)
% phi2 = phi3 = -17.878 (degrees) 
% phi4 = -64.778 (degrees) 
% 
% Velocity and acceleration analysis 
% 
% omega1 = [0,0,73.304] (rad/s)
% alpha1 = [0,0,0] (rad/s^2)
% 
% vB1 = vB2 = [ 3.761,-10.332,0] (m/s)
% |vB1| = 10.996 (m/s)
% aB1 = aB2 = [757.409,275.674,0] (m/s^2)
% |aB1| = 806.018 (m/s^2)
% 
% aB1n = [757.409,275.674,0] (m/s^2)
% aB1t = [0,0,0] (m/s^2)
% 
% vB3 = vC + omega3 x rCB = vB2 + vB3B2 => 
% x-axis: 0.0513*omega3z - 0.952*vB32 - 3.76 = 0 
% y-axis: 0.159*omega3z + 0.307*vB32 + 10.3 = 0 
% =>
% omega3z = -51.934 (rad/s)
% vB32 = -6.751 (m/s)
% 
% omega2=omega3 = [0,0,-51.934](rad/s)
% vB3B2 = [-6.425, 2.073,0] (m/s)
% 
% aB32cor = [215.270,667.364,0](m/s^2)
% 
% aB3=aC+alpha3xrCB-(omega3.omega3)rCB
% aB3=aB2+aB3B2+aB3B2cor =>
% x-axis:
% 0.0513*alpha3z - 0.952*aB32 - 1400.0 = 0
% y-axis:
% 0.307*aB32 + 0.159*alpha3z - 805.0 = 0
% =>
% alpha3z = 7157.347 (rad/s^2)
% aB32 = -1086.945 (m/s^2)
% 
% alpha2=alpha3=[0,0,7157.347](rad/s^2)
% aB3B2=[-1034.459,333.682,0] (m/s^2)
% 
% vD3 = vD4 = [ 3.189, 9.885,0] (m/s)
% |vD3| = 10.387 (m/s)
% aD3 = aD4 = [73.937,-1527.948,0] (m/s^2)
% |aD3| = 1529.736 (m/s^2)
% 
% vDC = [ 3.189, 9.885,0] (m/s)
% aDC = [73.937,-1527.948,0] (m/s^2)
% aDCn = [513.385,-165.601,0] (m/s^2)
% |aDCn| = 539.433 (m/s^2)
% aDCt = [-439.448,-1362.347,0] (m/s^2)
% |aDCt| = 1431.469 (m/s^2)
% 
% vF = [-7.972,-24.713,0] (m/s)
% aF = [-184.842,3819.871,0] (m/s^2)
% |vF| = 25.967 (m/s)
% |aF| = 3824.340 (m/s^2)
% 
% vE = vD + omega4 x (rE-rD) => 
% x-axis: 0.127*omega4z - 3.19 = 0 
% y-axis: 0.0597*omega4z + vEy - 9.89 = 0 
% =>
% omega4z = 25.176 (rad/s)
% vEy =  8.383 (m/s)
% 
% omega4 = [0,0,25.176] (rad/s)
% vE = [0, 8.383,0] (m/s)
% |vE| =  8.383 (m/s)
% 
% aE4=aD+alpha4xrDE-(omega4.omega4)rDE
% aE5=aE4 =>
% x-axis: 0.127*alpha4z - 112.0 = 0 
% y-axis: aEy + 0.0597*alpha4z + 1611.0 = 0 
% =>
% alpha4z = 882.339 (rad/s^2)
% aEy = -1660.866 (m/s^2)
% 
% alpha4=[0,0,882.339] (rad/s^2)
% aE = [0,-1660.866,0] (m/s^2)
% |aE| = 1660.866 (m/s^2)
% 
% vED = [-3.189,-1.502,0] (m/s)
% aED = [-73.937,-132.917,0] (m/s^2)
% aEDn = [37.814,-80.279,0] (m/s^2)
% |aEDn| = 88.739 (m/s^2)
% aEDt = [-111.751,-52.638,0] (m/s^2)
% |aEDt| = 123.527 (m/s^2)
