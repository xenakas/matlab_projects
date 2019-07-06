% example 5.3a

clear all; clc; close all

% Input data 
AB = 0.15;  % m 	
AC = 0.30;	% m 	
CD = 0.20;	% m		
DE = 0.14;	% m	
CF = 0.50;  % m
a  = 0.25;	% m
phi = 60*pi/180; % rad

% Position of joint A
xA=0; yA=0; rA=[xA yA 0]; 

% Position of joint B
xB=AB*cos(phi); yB=AB*sin(phi);
rB=[xB yB 0]; 

% Position of joint C
xC=-AC; yC=0; rC=[xC yC 0];

xE=-AC-a;

% Position of joint D 
% Distance formula: CD=constant
xDsol = sym('xDsol','real');
yDsol = sym('yDsol','real');
eqD1=(xC-xDsol)^2+(yC-yDsol)^2-CD^2 ;
% Distance formula: CD=constant
eqD2=(yC-yB)/(xC-xB)-(yDsol-yC)/(xDsol-xC);

% Simultaneously solve above equations
solD= solve(eqD1, eqD2, 'xDsol, yDsol');
% Two solutions for xD - vector form
xDpos= eval(solD.xDsol);
% Two solutions for yD - vector form
yDpos = eval(solD.yDsol);
% Separate the solutions in scalar form
xD1 = xDpos(1); 
xD2 = xDpos(2); 
yD1 = yDpos(1); 
yD2 = yDpos(2); 

% Select the correct position for D 
% for the given input angle
if(phi>=0 && phi<=pi/2)||(phi >= pi/2 && phi<=pi)
if yD1 <= yC xD=xD1; yD=yD1; else xD=xD2; yD=yD2;end
 else
if yD1 >= yC xD=xD1; yD=yD1; else xD=xD2; yD=yD2;end
end

rD = [xD yD 0]; % Position vector of D 

% Position of joint E 
% Distance formula: DE=constant
yEsol = sym('yEsol','real');
eqE1 = (xD-xE)^2 + (yD-yEsol)^2 - DE^2;
solE = solve(eqE1,'yEsol');
yEpositions=eval(solE);
yE1 = yEpositions(1); 
yE2 = yEpositions(2); 

if yE1 > yD 
    yE=yE1;
else
    yE=yE2;
end
rE = [xE yE 0]; % Position vector of E 

% Angles of the links with the horizontal
phi3 = atan((yC-yB)/(xC-xB));
phi2 = phi3;
phi4 = atan((yE-yD)/(xE-xD));

% Position of joint F 
CF = AB+AC+0.1;
xF = xC+CF*cos(phi2);
yF = yC+CF*sin(phi2);
rF = [xF,yF,0];

fprintf('Results \n\n')
fprintf('rA = [%6.3g, %6.3g, %g] (m)\n',rA)
fprintf('rB = [%6.3g, %6.3g, %g] (m)\n',rB)
fprintf('rC = [%6.3g, %6.3g, %g] (m)\n',rC)
fprintf('rD = [%6.3g, %6.3g, %g] (m)\n',rD)
fprintf('rE = [%6.3g, %6.3g, %g] (m)\n',rE)
fprintf('rF = [%6.3g, %6.3g, %g] (m)\n',rF)
fprintf...
    ('phi3 = %6.3g (degrees) \n',phi3*180/pi)
fprintf...
    ('phi4 = %6.3g (degrees) \n',phi4*180/pi)

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
% x limits to the specified values
xlim([-la la]);
% x limits to the specified values
ylim([-la la]);

% end of program