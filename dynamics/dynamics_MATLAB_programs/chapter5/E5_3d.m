
% velocity and acceleration analysis
% R-RTR-RRT
% contour method

clear all; clc; close all

AB = 0.15;  % m 	
AC = 0.30;	% m 	
CD = 0.20;	% m		
DE = 0.14;	% m	
CF = 0.50;  % m
a  = 0.25;	% m
phi1 = 200*pi/180;	% rad

% joint A
xA = 0; yA = 0; rA = [xA,yA,0];
% joint C
xC = -AC; yC = 0; rC = [xC yC 0];
% joint B
xB=AB*cos(phi1); yB=AB*sin(phi1); 
rB = [xB yB 0];  
% link 3
phi3 = atan((yB-yC)/(xB-xC));
phi2 = phi3;
% joint F
xF = xC+CF*cos(phi2);
yF = yC+CF*sin(phi2);
rF = [xF,yF,0];
% joint D
xD = xC - CD*cos(phi3);
yD = yC - CD*sin(phi3);
rD = [xD yD 0]; 
% joint E
phi5 = pi/2;
xE = -AC-a;
yEs = sym('yEs','real');
eqnE=(xE-xD)^2+(yEs-yD)^2-DE^2';
solE = solve(eqnE,'yEs');
yE1 = eval(solE(1));
yE2 = eval(solE(2));
if yE1 > yC yE = yE1;
else yE = yE2;
end
rE = [xE yE 0];
% link 4
phi4 = atan((yE-yD)/(xE-xD));
% plot
plot(...
[xA,xB],[yA,yB],'r-o',...
[xD,xF],[yD,yF],'b-o',...
[xC,xC],[yC,yC],'k-o',...
[xD,xE],[yD,yE],'g-o')
axis square
xlabel('x (m)')
ylabel('y (m)')
title('positions for \theta = 200 (deg)')
text(xA,yA,'  A')
text(xB,yB,'  B')
text(xC,yC,'  C')
text(xD,yD,'  D')
text(xE,yE,'  E')
text(xF,yF,'  F')
grid

fprintf('\n') 
fprintf('Velocity and acceleration analysis\n\n')  
fprintf('Contour method \n') 
fprintf('\n') 

n = 700.; % rpm
omega1 = [ 0 0 pi*n/30 ]; alpha1 = [0 0 0 ]; 
fprintf...
('omega1 = [%g,%g,%6.3f] (rad/s)\n',omega1) 
fprintf...
('alpha1 = [%g,%g,%g] (rad/s^2)\n',alpha1)   
fprintf('\n') 

vA = [0 0 0]; aA = [0 0 0];

vB1 = vA + cross(omega1,rB); vB2 = vB1; 
aB1 = aA + cross(alpha1,rB) - ...
    dot(omega1,omega1)*rB;

vB2 = vB1;
aB2 = aB1;
fprintf...
('vB1=vB2=[%6.3f,%6.3f,%g] (m/s)\n',vB1) 
fprintf...
('aB1=aB2=[%6.3f,%6.3f,%g] (m/s^2)\n',aB1) 
fprintf('\n')

fprintf('\n') 
fprintf('Contour I: 0-1-2-3-0 \n\n')  
fprintf('Relative velocities \n') 
fprintf('\n') 

omega10 = omega1; 
omega21v = [ 0 0 sym('omega21z','real')]; 
omega03v = [ 0 0 sym('omega03z','real')]; 
v32v = ...
sym('vB32','real')*[cos(phi2) sin(phi2) 0]; 

eqIomega = omega10 + omega21v + omega03v; 
eqIvz=eqIomega(3);
eqIv = ...
 cross(rB,omega21v)+cross(rC,omega03v)+v32v;
eqIvx=eqIv(1);
eqIvy=eqIv(2);

digits 3
Ivz=vpa(eqIvz);
fprintf('%s = 0 \n', char(Ivz)) 
Ivx=vpa(eqIvx);
fprintf('%s = 0 \n', char(Ivx)) 
Ivy=vpa(eqIvy);
fprintf('%s = 0 \n', char(Ivy)) 

solIv=solve(eqIvz,eqIvx,eqIvy);
omega21 = [ 0 0 eval(solIv.omega21z) ]; 
omega03 = [ 0 0 eval(solIv.omega03z) ]; 
vB3B2 = ...
 eval(solIv.vB32)*[cos(phi2) sin(phi2) 0]; 

fprintf...
('omega21=[%g,%g,%6.3f] (rad/s)\n',omega21) 
fprintf...
('omega03=[%g,%g,%6.3f] (rad/s)\n',omega03) 
fprintf...
('vB32 = %6.3f (m/s)\n', eval(solIv.vB32)) 
fprintf...
('vB3B2=[%6.3f,%6.3f,%d] (m/s)\n', vB3B2)  
fprintf('\n') 

fprintf('Absolute velocities \n\n') 
omega30 = - omega03;
omega20 = omega30;
vC = [0 0 0 ]; 
vD3 = vC + cross(omega30,rD-rC);
fprintf...
('omega20=omega30=[%d,%d,%6.3f](rad/s)\n',omega30) 
fprintf...
('vD3=vD4 =[%6.3f,%6.3f,%d] (m/s)\n', vD3) 

fprintf('\n') 
fprintf('Relative accelerations \n') 
fprintf('\n') 

alpha10 = alpha1; 
alpha21v = [ 0 0 sym('alpha21z','real') ]; 
alpha03v = [ 0 0 sym('alpha03z','real') ]; 
a32v=sym('aB32','real')*[cos(phi2) sin(phi2) 0]; 

eqIalpha = alpha10 + alpha21v + alpha03v; 
eqIaz=eqIalpha(3);
eqIa=cross(rB,alpha21v)+cross(rC,alpha03v)+...
 a32v+2*cross(omega20,vB3B2)-...
 dot(omega1,omega1)*rB-dot(omega30,omega30)*(rC-rB);
eqIax=eqIa(1);
eqIay=eqIa(2);

Iaz=vpa(eqIaz);
fprintf('%s=0 \n',char(Iaz)) 
Iax=vpa(eqIax);
fprintf('%s=0 \n',char(Iax)) 
Iay=vpa(eqIay);
fprintf('%s=0 \n',char(Iay))  

solIa=solve(eqIaz,eqIax,eqIay);
alpha21 = [ 0 0 eval(solIa.alpha21z) ]; 
alpha03 = [ 0 0 eval(solIa.alpha03z) ]; 
aB3B2=eval(solIa.aB32)*[cos(phi2) sin(phi2) 0]; 

fprintf...
('alpha21=[%d,%d,%6.3f] (rad/s^2)\n',alpha21) 
fprintf...
('alpha03=[%d,%d,%6.3f] (rad/s^2)\n',alpha03) 
fprintf('aB32=%6.3f (m/s^2)\n',eval(solIa.aB32)) 
fprintf('aB3B2=[%6.3f,%6.3f,%d] (m/s^2)\n',aB3B2)  
fprintf('\n') 

fprintf('Absolute accelerations \n\n') 
alpha30 = - alpha03;
alpha20 = alpha30;
aC = [0 0 0 ]; 
aD3=aC+cross(alpha30,rD-rC)-...
    dot(omega20,omega20)*(rD-rC);
fprintf...
('alpha20=alpha30=[%d,%d,%6.3f] (rad/s^2)\n',alpha30) 
fprintf('aD3=aD4 = [%6.3f,%6.3f,%d] (m/s^2)\n',aD3) 

fprintf('\n') 
fprintf('Contour II: 0-3-4-5-0 \n\n')
fprintf('Relative velocities \n') 
fprintf('\n') 

omega43v = [ 0 0 sym('omega43z','real') ]; 
omega54v = [ 0 0 sym('omega54z','real') ]; 
v05v = [0 sym('v05','real') 0]; 

eqIIomega = omega30 + omega43v + omega54v; 
eqIIvz=eqIIomega(3);
eqIIv=cross(rC,omega30)+cross(rD,omega43v)+...
cross(rE,omega54v)+v05v;
eqIIvx=eqIIv(1);
eqIIvy=eqIIv(2);

IIvz=vpa(eqIIvz); fprintf('%s = 0 \n', char(IIvz)) 
IIvx=vpa(eqIIvx); fprintf('%s = 0 \n', char(IIvx)) 
IIvy=vpa(eqIIvy); fprintf('%s = 0 \n', char(IIvy)) 

solIIv=solve(eqIIvz,eqIIvx,eqIIvy);
omega43 = [ 0 0 eval(solIIv.omega43z) ]; 
omega54 = [ 0 0 eval(solIIv.omega54z) ]; 
vE05 = [0 eval(solIIv.v05) 0]; 

fprintf...
('omega43=[%d,%d,%6.3f] (rad/s)\n',omega43)
fprintf...
('omega54=[%d,%d,%6.3f] (rad/s)\n',omega54)
fprintf('vE05=%6.3f (m/s)\n',eval(solIIv.v05))
fprintf('vE0E5 =[%d,%6.3f,%d] (m/s)\n',vE05) 
fprintf('\n')

fprintf('Absolute velocities \n\n')
omega40=omega30+omega43;
fprintf...
('omega40=[%d,%d,%6.3f] (rad/s)\n',omega40)
fprintf('vE5=[%d,%6.3f,%d] (m/s)\n',-vE05) 

fprintf('\n') 
fprintf('Relative accelerations \n\n')
 
alpha43v = [ 0 0 sym('alpha43z','real') ]; 
alpha54v = [ 0 0 sym('alpha54z','real') ]; 
a05v=[0,sym('a05','real'),0]; 

eqIIalpha = alpha30 + alpha43v + alpha54v; 
eqIIaz=eqIIalpha(3);
eqIIa=cross(rC,alpha30)+cross(rD,alpha43v)+...
cross(rE,alpha54v)+a05v-...
dot(omega30,omega30)*(rD-rC)-...
dot(omega40,omega40)*(rE-rD);
eqIIax=eqIIa(1);
eqIIay=eqIIa(2);

IIaz=vpa(eqIIaz); fprintf('%s=0 \n',char(IIaz)) 
IIax=vpa(eqIIax); fprintf('%s=0 \n',char(IIax))
IIay=vpa(eqIIay); fprintf('%s=0 \n',char(IIay)) 

solIIa=solve(eqIIaz,eqIIax,eqIIay);
alpha43 = [0 0 eval(solIIa.alpha43z)]; 
alpha54 = [0 0 eval(solIIa.alpha54z)]; 
aE0E5=[0,eval(solIIa.a05),0]; 

fprintf...
('alpha43=[%d,%d,%6.3f] (rad/s^2)\n',alpha43) 
fprintf...
('alpha54=[%d,%d,%6.3f] (rad/s^2)\n',alpha54) 
fprintf('aE05=%6.3f (m/s^2)\n',eval(solIIa.a05)) 
fprintf('aE0E5=[%d,%6.3f,%d] (m/s^2)\n',aE0E5)  
fprintf('\n') 

fprintf('Absolute accelerations \n\n')
alpha40 = alpha30+alpha43;

fprintf...
('alpha40=[%d,%d,%6.3f] (rad/s^2)\n',alpha40)
fprintf('aE5=[%d,%6.3f,%d] (m/s^2)\n',-aE0E5)  


% end of program