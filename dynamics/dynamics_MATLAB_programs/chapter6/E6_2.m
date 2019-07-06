% example 6.2
% R-RTR-RRT Mechanism
% contour method

% kinematics 
R_RTR_RRT_va
% calculate velocities and accelerations

rC1 = rB/2; 
fprintf('rC1 = [%6.3f,%6.3f,%d] (m)\n', rC1)
rC2 = rB; 
fprintf('rC2 = [%6.3f,%6.3f,%d] (m)\n', rC2)
rC3 = (rD+rC)/2; 
fprintf('aC3 = [%6.3f,%6.3f,%d] (m)\n', rC3)
rC4 = (rC+rE)/2; 
fprintf('rC4 = [%6.3f,%6.3f,%d] (m)\n', rC4)
rC5 = rE; 
fprintf('rC5 = [%6.3f,%6.3f,%d] (m)\n', rC5)

aC1 = aB/2; 
fprintf('aC1 = [%6.3f,%6.3f,%d] (m/s^2)\n', aC1)
aC2 = aB; 
fprintf('aC2 = [%6.3f,%6.3f,%d] (m/s^2)\n', aC2)
aD = 0;
aC3 = (aD+aC)/2; 
fprintf('aC3 = [%6.3f,%6.3f,%d] (m/s^2)\n', aC3)
aC4 = (aC+aE)/2; 
fprintf('aC4 = [%6.3f,%6.3f,%d] (m/s^2)\n', aC4)
aC5 = aE; 
fprintf('aC5 = [%6.3f,%d,%d] (m/s^2)\n', aC5)

fprintf('\n')
fprintf('Dynamic force analysis  \n')
fprintf('Dyad method  \n\n')

h = 0.01; % height of the bar
d = 0.001; % depth of the bar
hSlider = 0.02; % height of the slider 
wSlider = 0.05; % depth of the slider 
rho = 8000; % density of the material
g = 9.807; % gravitational acceleration

fprintf('Inertia forces and inertia moments\n\n')

fprintf('Link 1  \n')
m1 = rho*AB*h*d;
Fin1 = -m1*aC1;
G1 = [0,-m1*g,0];
IC1 = m1*(AB^2+h^2)/12;
alpha1 = [0 0 0];
Min1 = -IC1*alpha1;
fprintf('m1 = %6.3f (kg)\n', m1)
fprintf('Fin1=-m1 aC1=[%6.3f,%6.3f,%d] (N)\n',Fin1)
fprintf('G1 = - m1 g = [%d,%6.3f,%d] (N)\n', G1)
fprintf('IC1 = %g (kg m^2)\n', IC1)
fprintf('Min1=-IC1 alpha1=[%d,%d,%d](N m)\n',Min1)
fprintf('\n')

fprintf('Link 2 \n')
m2 = rho*hSlider*wSlider*d;
Fin2 = -m2*aC2;
G2 = [0,-m2*g,0];
IC2 = m2*(hSlider^2+wSlider^2)/12;
Min2 = -IC2*alpha2;
fprintf('m2 = %6.3f (kg)\n', m2)
fprintf('Fin2=-m2 aC2=[%6.3f,%6.3f,%d] (N)\n',Fin2)
fprintf('G2 = - m2 g = [%d,%6.3f,%d] (N)\n', G2)
fprintf('IC2 = %6.3d (kg m^2)\n', IC2)
fprintf('Min2=-IC2 alpha2=[%d,%d,%6.3d](N m)\n',Min2)
fprintf('\n')

fprintf('Link 3 \n')
m3 = rho*CD*h*d;
Fin3 = -m3*aC3;
G3 = [0,-m3*g,0];
IC3 = m3*(CD^2+h^2)/12;
Min3 = -IC3*alpha3;
fprintf('m3 = %6.3f (kg)\n', m3)
fprintf('Fin3 =-m3 aC3=[%6.3f,%6.3f,%g] (N)\n',Fin3)
fprintf('G3 = - m3 g = [%g, %6.3f, %g] (N)\n', G3)
fprintf('IC3 = %6.3d (kg m^2)\n', IC3)
fprintf('Min3=-IC3 alpha3=[%d,%d,%6.3d](N m)\n',Min3)
fprintf('\n')

fprintf('Link 4 \n')
m4 = rho*CE*h*d;
Fin4 = -m4*aC4;
G4 = [0,-m4*g,0];
IC4 = m2*(CE^2+h^2)/12;
Min4 = -IC4*alpha4;
fprintf('m4 = %6.3f (kg)\n', m4)
fprintf('Fin4=-m4 aC4=[%6.3f,%6.3f,%d] (N)\n',Fin4)
fprintf('G4 = - m4 g = [%d,%6.3f,%d] (N)\n', G4)
fprintf('IC4 = %6.3d (kg m^2)\n', IC4)
fprintf('Min4=-IC4 alpha4=[%d,%d,%6.3d](N m)\n',Min4)
fprintf('\n')

fprintf('Link 5 \n')
m5 = rho*hSlider*wSlider*d;
Fin5 = -m5*aC5;
G5 = [0,-m5*g,0];
IC5 = m5*(hSlider^2+wSlider^2)/12;
Min5 = -IC5*alpha5;
fprintf('m5 = %6.3f (kg)\n', m5)
fprintf('Fin5=-m5 aC5=[%6.3f,%d,%d] (N)\n',Fin5)
fprintf('G5 = - m5 g = [%g,%6.3f,%g] (N)\n',G5)
fprintf('IC5 = %6.3d (kg m^2)\n', IC5)
fprintf('Min5=-IC5 alpha5=[%d,%d,%d](N m)\n',Min5)
fprintf('\n')

Fext = -sign(vE(1))*[500,0,0];
fprintf('Fext = [%d, %6.3f, %d] (N)\n',Fext)
fprintf('\n')

% denote
F1 = Fin1 + G1;
M1 = Min1;
F2 = Fin2 + G2;
M2 = Min2;
F3 = Fin3 + G3;
M3 = Min3;
F4 = Fin4 + G4;
M4 = Min4;
F5 = Fin5 + G5;
M5 = Min5;

% F34 
F34x=sym('F34x','real');
F34y=sym('F34y','real');
F34=[F34x, F34y, 0];
% sum M_E for (4)
% EC x F34 + EC4 x F4 + M4 = 0
eqME4=cross(rC-rE,F34)+cross(rC4-rE,F4)+M4;
% sum F for (4&5) on x-axis
% (F5+Fext+F4+F34).i = 0
eqF45x=dot(F5+Fext+F4+F34,[1 0 0]);
solF34=solve(eqME4(3),eqF45x);
F34xn=eval(solF34.F34x);
F34yn=eval(solF34.F34y);
F34n=[F34xn, F34yn, 0];
fprintf('F34 = [%6.3f,%6.3f,%d] (N)\n',F34n)
fprintf('\n')

% F54 
F54x=sym('F54x','real');
F54y=sym('F54y','real');
F54=[F54x, F54y, 0];
% sum M_C for (4)
% CE x F54 + CC4 x F4 + M4 = 0
eqMC4=cross(rE-rC,F54)+cross(rC4-rC,F4)+M4;
% sum F for (5) on x-axis
% [F5+Fext+(-F54)].i = 0
eqF5x=dot(F5+Fext-F54,[1 0 0]);
solF54=solve(eqMC4(3),eqF5x);
F54xn=eval(solF54.F54x);
F54yn=eval(solF54.F54y);
F54n=[F54xn, F54yn, 0];
fprintf('F54 = [%6.3f,%6.3f,%d] (N)\n',F54n)
fprintf('\n')

% F05 
% sum M_E for (5)
% EP x F05 = 0 => E=P
% F05 acts at E
F05y=sym('F05y','real');
F05=[0, F05y, 0];
% sum M_C for (4&5)
% CE x (F05+F5+Fext) + CC4 x F4 + M4 = 0
eqMC45=cross(rE-rC,F05+F5+Fext)+...
       cross(rC4-rC,F4)+M4;
solF05=solve(eqMC45(3));
F05yn=eval(solF05);
F05n=[0, F05yn, 0];
fprintf('F05 = [%d,%6.3f,%d] (N)\n',F05n)
fprintf('\n')

F43 = -F34n;

% F32 acts at Q(xQ,yQ)
F32x=sym('F32x','real');
F32y=sym('F32y','real');
F32=[F32x, F32y, 0];
xQ = sym('xQ','real');
yQ = sym('yQ','real');
rQ = [xQ yQ 0];
% F32 perpendicular to BD
% F32.BD = 0
eqF32BD = dot(F32,rB-rD); 
% Q is on the line BD
% DB x DQ = 0
eqQ = cross(rB-rD,rQ-rD);
% sum M_B for (2)
% BQ x F32 + M2 = 0
eqMB2=cross(rQ-rB,F32)+M2;
% sum M_D for (3)
% DQ x (-F32) + DC3 x F3 + DC x F43 + M3 = 0
eqMD3=cross(rQ-rD,-F32)+cross(rC3-rD,F3)...
      +cross(rC-rD,F43)+M3;
solF32=solve(eqF32BD,eqQ(3),eqMB2(3),eqMD3(3));
F32xn=eval(solF32.F32x);
F32yn=eval(solF32.F32y);
F32n=[F32xn, F32yn, 0];
fprintf('F32 = [%6.3f,%6.3f,%d] (N)\n',F32n)
xQn=eval(solF32.xQ);
yQn=eval(solF32.yQ);
rQn=[xQn, yQn, 0];
fprintf('rQ = [%6.6f,%6.6f,%d] (m)\n',rQn)
fprintf('\n')

% F03
F03x=sym('F03x','real');
F03y=sym('F03y','real');
F03=[F03x, F03y, 0];
% sum F for (3) on DB
% (F03 + F3 + F43).DB = 0
eqF3DB=dot(F03+F3+F43,rB-rD);
% sum M_B for (3&2)
% BD x F03 + BC x F43 + BC3 x F3 +M3+M2 = 0
eqMB32=cross(rD-rB,F03)+cross(rC-rB,F43)...
      +cross(rC3-rB,F3)+M3+M2;
solF03=solve(eqF3DB,eqMB32(3));
F03xn=eval(solF03.F03x);
F03yn=eval(solF03.F03y);
F03n=[F03xn, F03yn, 0];
fprintf('F03 = [%6.3f,%6.3f,%d] (N)\n',F03n)
fprintf('\n')

% F12
F12x=sym('F12x','real');
F12y=sym('F12y','real');
F12=[F12x, F12y, 0];
% sum F for (2) on DB
% (F12 + F2).DB = 0
eqF2DB=dot(F12+F2,rB-rD);
% sum M_D for (2&3)
% DB x (F12+F2) + DC x F43 + DC3 x F3 +M3+M2 = 0
eqMD23=cross(rB-rD,F12+F2)+cross(rC-rD,F43)...
      +cross(rC3-rD,F3)+M3+M2;
solF12=solve(eqF2DB,eqMD23(3));
F12xn=eval(solF12.F12x);
F12yn=eval(solF12.F12y);
F12n=[F12xn, F12yn, 0];
fprintf('F12 = [%6.3f,%6.3f,%d] (N)\n',F12n)
fprintf('\n')

% F01
F01x=sym('F01x','real');
F01y=sym('F01y','real');
F01=[F01x, F01y, 0];
Mz=sym('Mz','real');
M=[0, 0, Mz];
% sum M_B for (1)
% BA x F01 + BC1 x F1 + M1 + M = 0
eqMB1=cross(-rB,F01)+cross(rC1-rB,F1)+M1+M;
% sum F for (1&2) on DB
% (F01 + F1 + F2).DB = 0
eqF12DB=dot(F01+F1+F2,rB-rD);
% sum M_D for (1&2&3)
% DA x F01 + DC1 x F1 + DB x F2 + DC x F43 
% + DC3 x F3 + M3 + M2 + M1 + M = 0
eqMD123=cross(-rD,F01)+cross(rC1-rD,F1)+...
        cross(rB-rD,F2)+cross(rC-rD,F43)+...
        cross(rC3-rD,F3)+M3+M2+M1+M;
solF01=solve(eqMB1(3),eqF12DB,eqMD123(3));
F01xn=eval(solF01.F01x);
F01yn=eval(solF01.F01y);
Mzn=eval(solF01.Mz);
F01n=[F01xn, F01yn, 0];
Mn=[0, 0, Mzn];
fprintf('F01 = [%6.3f,%6.3f,%d] (N)\n',F01n)
fprintf('\n')
fprintf('Mmot = [%d,%d,%6.3f] (N m)\n',Mn)
fprintf('\n')

% end of program


% F34 = [-500.685,124.998,0] (N)
% 
% F54 = [500.140,-125.037,0] (N)
% 
% F05 = [0,-124.959,0] (N)
% 
% F32 = [625.114,-163.282,0] (N)
% rQ = [0.141421,0.141421,0] (m)
% 
% F03 = [123.647,-38.376,0] (N)
% 
% F12 = [-625.610,162.864,0] (N)
% 
% F01 = [-626.106,162.525,0] (N)
% 
% Mmot = [0,0,111.518] (N m)