% example 5.3b

clear all; clc; close all

% Input data 
AB = 0.15;  % m 	
AC = 0.30;	% m 	
CD = 0.20;	% m		
DE = 0.14;	% m	
CF = 0.50;  % m
a  = 0.25;	% m
phi0 = 60*pi/180; % rad

% Position of joint A
xA=0; yA=0; rA=[xA yA 0]; 

% Position of joint C
xC=-AC; yC=0; rC=[xC yC 0];

xE=-AC-a;

% allocate/initialize the matrix 
% to have 12 frames
 M = moviein(12);

% at the initial moment phi0 => 
incr = 0 ; 

% the step has to be small for this method
step=pi/15; 
for phi=phi0:step:2*pi+phi0,
% fprintf('phi =  %g deegres \n', phi*180/pi)

xB=AB*cos(phi); yB=AB*sin(phi);
rB=[xB yB 0]; 

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

% select the correct position for D 
%    only for increment == 0
% the selection process is automatic 
%    for all the other steps

if incr == 0
%    
if(phi>=0 && phi<=pi)
if yD1 <= yC xD=xD1; yD=yD1; 
else xD=xD2; yD=yD2; end
 else
if yD1 >= yC xD=xD1; yD=yD1; 
else xD=xD2; yD=yD2; end
end
%   
else
  dist1 = Dist(xD1,yD1,xDold,yDold);
  dist2 = Dist(xD2,yD2,xDold,yDold);
  if dist1 < dist2 xD=xD1; yD=yD1; 
   else xD=xD2; yD=yD2; 
   end
end
xDold=xD;
yDold=yD;

% Position of joint E 
% Distance formula: DE=constant
yEsol = sym('yEsol','real');
eqE1 = (xD-xE)^2+(yD-yEsol)^2-DE^2;
solE = solve(eqE1,'yEsol');
yEpos=eval(solE);
yE1 = yEpos(1); 
yE2 = yEpos(2); 

if incr == 0
   if yE1 > yC 
    yE=yE1;
   else
    yE=yE2;
   end
else
  dist1 = Dist(xE,yE1,xE,yEold);
  dist2 = Dist(xE,yE2,xE,yEold);
  if dist1 < dist2 yE=yE1; 
   else yE=yE2; 
   end
end
yEold=yE;
      
rD = [xD yD 0];
rE = [xE yE 0]; 

% Angles of links with x-axis
phi3 = atan((yC-yB)/(xC-xB));
phi2 = phi3;
phi4 = atan((yE-yD)/(xE-xD));

% Position of joint F 
CF = AB+AC+0.1;
xF = xC+CF*cos(phi2);
yF = yC+CF*sin(phi2);
rF = [xF,yF,0];

% Graphic of the mechanism
%axis manual
% axis equal

la=-xE+0.05;
GR=plot(...
[xA,xB],[yA,yB],'r-o',...
[xD,xF],[yD,yF],'k-o',...
[xC,xC],[yC,yC],'k-o',...
[xD,xE],[yD,yE],'b-o',...
[xE,xE],[-la,la],'r--');
hold on
grid
text(xA,yA,'   A')
text(xB,yB,'   B')
text(xC,yC,'   C')
text(xD,yD,'  D')
text(xE,yE,'  E')
text(xF,yF,'  F')
% x limits 
xlim([-la la]);
% y limits 
ylim([-la la]);

incr=incr+1;
% record the movie
 M(:,incr) = getframe; 
end

 movie2avi(M,'R_RTRRRT.avi');

% end of program