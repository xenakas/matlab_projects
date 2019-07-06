% example 6.1
% R-RTR-RRT Mechanism
% dyad method

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
fprintf('Fin3 =- m3 aC3=[%6.3f,%6.3f,%g] (N)\n',Fin3)
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
fprintf('Fin4=- m4 aC4=[%6.3f,%6.3f,%d] (N)\n',Fin4)
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
fprintf('Fin5 =-m5 aC5=[%6.3f,%d,%d] (N)\n', Fin5)
fprintf('G5 = - m5 g = [%g, %6.3f, %g] (N)\n', G5)
fprintf('IC5 = %6.3d (kg m^2)\n', IC5)
fprintf('Min5=-IC5 alpha5=[%d,%d,%d](N m)\n',Min5)
fprintf('\n')

Fext = -sign(vE(1))*[500,0,0];
fprintf('Fext = [%d, %6.3f, %d] (N)\n',Fext)
fprintf('\n')

% Dyad RRT (4 & 5)
F5 = Fin5 + G5 + Fext;
F4 = Fin4 + G4;
M5 = Min5;
M4 = Min4;
rN = [0 yE 0];
FD45=forceRRT(rC,rN,rE,rC4,F4,M4,F5,M5);
F34x = FD45(1);
F34y = FD45(2);
F05x = FD45(3);
F05y = FD45(4);
F54x = FD45(5);
F54y = FD45(6);
F34 = [F34x F34y 0];
F05 = [F05x F05y 0];
F54 = [F54x F54y 0];
fprintf('F05 = [%6.3f,%6.3f,%d] (N)\n',F05)
fprintf('F54 = [%6.3f,%6.3f,%d] (N)\n',F54)
fprintf('F34 = [%6.3f,%6.3f,%d] (N)\n',F34)
fprintf('\n')

% Dyad RTR (2 & 3)
F43 = -F34;
F3 = Fin3 + G3 + F43;
F2 = Fin2 + G2;
M3 = Min3 + cross(rC-rC3,F43);
M2 = Min2;

FD23=forceRTR(rD,rB,rC3,F3,M3,F2,M2);

F03x = FD23(1);
F03y = FD23(2);
F12x = FD23(3);
F12y = FD23(4);
F32x = FD23(5);
F32y = FD23(6);

F03 = [F03x F03y 0];
F12 = [F12x F12y 0];
F32 = [F32x F32y 0];

xQ = FD23(7);
yQ = FD23(8);
rQ = [xQ yQ 0];

fprintf('F03 = [%6.3f,%6.3f,%d] (N)\n',F03)
fprintf('F12 = [%6.3f,%6.3f,%d] (N)\n',F12)
fprintf('F32 = [%6.3f,%6.3f,%d] (N)\n',F32)
fprintf('rQ = [%6.5f,%6.5f,%d] (m)\n',rQ)
fprintf('\n')

% driver
FDR=forceDR(rA,rB,rC1,Fin1+G1,Min1,-F12);
F01 = [FDR(1) FDR(2) 0];
Mmot= [0 0 FDR(3)];
fprintf('F01 = [%6.3f,%6.3f,%d] (N)\n',F01)
fprintf('Mmot = [%d,%d,%6.3f] (N m)\n',Mmot)

% end of program


% 
% 
% Results 
% 
% phi = phi1 = 45 (degrees) 
% rA = [0, 0, 0] (m)
% rD = [ 0.000, -0.400, 0] (m)
% rB = [ 0.141,  0.141, 0] (m)
% rC = [ 0.177,  0.277, 0] (m)
% phi2 = phi3 = 75.361 (degrees) 
% rE = [-0.114,  0.350, 0] (m)
% phi4 = 165.971 (degrees) 
% 
% 
% Velocity and acceleration analysis 
% 
% omega1 = [0, 0,   20.9] (rad/s)
% alpha1 = [0, 0, 0] (rad/s^2)
% 
% vB = vB1 = vB2 = [ -2.96,   2.96, 0] (m/s)
% aB = aB1 = aB2 = [   -62,    -62, 0] (m/s^2)
% 
% omega2=omega3 = [0,0,  6.46](rad/s)
% alpha2=alpha3 = [0,0,  30.4](rad/s^2)
% vC = [-4.374,  1.143, 0] (m/s)
% aC = [-27.947, -22.882, 0] (m/s^2)
% 
% vE = [-4.660,  0.000, 0] (m/s)
% aE = [-17.464,  0.000, 0] (m/s^2)
% 
% omega4 = [0,0,  3.93](rad/s)
% alpha4 = [0,0, -82.5](rad/s^2)
% rC1 = [ 0.071, 0.071,0] (m)
% rC2 = [ 0.141, 0.141,0] (m)
% aC3 = [ 0.088,-0.061,0] (m)
% rC4 = [ 0.031, 0.314,0] (m)
% rC5 = [-0.114, 0.350,0] (m)
% aC1 = [-31.017,-31.017,0] (m/s^2)
% aC2 = [-62.034,-62.034,0] (m/s^2)
% aC3 = [-13.974,-11.441,0] (m/s^2)
% aC4 = [-22.705,-11.441,0] (m/s^2)
% aC5 = [-17.464, 0.000,0] (m/s^2)
% 
% Dynamic force analysis  
% Dyad method  
% 
% Fe = [500,  0.000, 0] (N)
% 
% Inertia forces and inertia moments
% 
% Link 1  
% m1 =  0.016 (kg)
% m1 aC1 = [-0.496,-0.496,0] (N)
% Fin1= - m1 aC1 = [ 0.496, 0.496,0] (N)
% G1 = - m1 g = [0,-0.157,0] (N)
% IC1 = 5.34667e-05 (kg m^2)
% IC1 alpha1=[ 0.000, 0.000,0](N m)
% Min1=-IC1 alpha1=[0, 0, 0](N m)
% 
% Link 2 
% m2 =  0.008 (kg)
% m2 aC2 = [-0.496,-0.496,0] (N)
% Fin2= - m2 aC2 =[ 0.496, 0.496,0] (N)
% G2 = - m2 g = [0,-0.078,0] (N)
% IC2 = 1.933e-06 (kg m^2)
% IC2 alpha2=[0,0,5.871e-05](N m)
% Min2=-IC2 alpha2=[0,0,-5.871e-05](N m)
% 
% Link 3 
% m3 =  0.056 (kg)
% m3 aC3 = [-0.783,-0.641,0] (N)
% Fin3 = - m3 aC3 =[ 0.783, 0.641,-0] (N)
% G3 = - m3 g = [0, -0.549, 0] (N)
% IC3 = 2.287e-03 (kg m^2)
% IC3 alpha3=[0,0,6.945e-02](N m)
% Min3=-IC3 alpha3=[0,0,-6.945e-02](N m)
% 
% Link 4 
% m4 =  0.024 (kg)
% m4 aC4 = [-0.545,-0.275,0] (N)
% Fin4= - m4 aC4 =[ 0.545, 0.275,0] (N)
% G4 = - m4 g = [0,-0.235,0] (N)
% IC4 = 6.007e-05 (kg m^2)
% IC4 alpha4=[0,0,-4.954e-03](N m)
% Min4=-IC4 alpha4=[0,0,4.954e-03](N m)
% 
% Link 5 
% m5 =  0.008 (kg)
% m5 aC5 = [-0.140, 0.000,0] (N)
% Fin5 = - m5 aC5 =[ 0.140,0,0] (N)
% G5 = - m5 g = [0, -0.078, 0] (N)
% IC5 = 1.933e-06 (kg m^2)
% IC5 alpha5=[0, 0, 0](N m)
% Min5=-IC5 alpha5=[0, 0, 0](N m)
% 
% F05 = [ 0.000,-124.959,0] (N)
% F54 = [500.140,-125.037,0] (N)
% F34 = [-500.685,124.998,0] (N)
% 
% F03 = [123.647,-38.376,0] (N)
% F12 = [-625.610,162.864,0] (N)
% F32 = [625.114,-163.282,0] (N)
% rQ = [0.14142,0.14142,0] (m)
% 
% F01 = [-626.106,162.525,0] (N)
% Mmot = [0,0,111.518] (N m)
