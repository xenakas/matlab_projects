% example 5.1

clear all; clc; close all

syms a b c omega x_G x_F y_G y_F z_G z_F

% omega=[omega_x, omega_y omega_z]
theta_x=a/sqrt(a^2+b^2+c^2);
theta_y=b/sqrt(a^2+b^2+c^2);
theta_z=c/sqrt(a^2+b^2+c^2);

omega_x=omega*theta_x;
omega_y=omega*theta_y;
omega_z=omega*theta_z;

omega_v=[omega_x omega_y omega_z];

% coordinates of  O, A, B, C, D, E, F, G
x_O=0; y_O=0; z_O=0;
x_A=a; y_A=0; z_A=0;
x_B=0; y_B=b; z_B=0;
x_C=0; y_C=0; z_C=c;
x_D=a; y_D=b; z_D=0;
x_E=a; y_E=0; z_E=c;
x_F=a; y_F=b; z_F=c;
x_G=0; y_G=b; z_G=c;

r_OA=[x_A y_A z_A];
% fprintf('r_OA =\n')
% pretty(r_OA); fprintf('\n\n')

% velocity of origin O
v_O = [0 0 0];

% velocity of A
v_A=v_O+cross(omega_v,r_OA);
fprintf('v_a = \n')
pretty(v_A); fprintf('\n\n')

% magnitude of v_A
m_v_A=simple(normvec(v_A(1),v_A(2),v_A(3)));
fprintf('|v_A| = \n')
pretty(m_v_A); fprintf('\n\n')

c_x=x_G-x_F;
c_y=y_G-y_F;
c_z=z_G-z_F;

r_FG=[x_G-x_F y_G-y_F z_G-z_F];
fprintf('r_FG =\n')
pretty(r_FG); fprintf('\n\n')

% velocity of F
v_F = [0 0 0];

% velocity of G
v_G=v_F+cross(omega_v,r_FG);
% fprintf('v_G = \n')
% pretty(v_G); fprintf('\n\n')

% magnitude of v_G
m_v_G=simple(normvec(v_G(1),v_G(2),v_G(3)));
fprintf('|v_G| = \n')
pretty(m_v_G); fprintf('\n\n')

delta=m_v_A-m_v_G;
fprintf('|v_A|-|v_G| = %s \n',char(delta))

% show and animate the prism

% numerical values
a=7; b=20; c=5; % (m)
x_O=0; y_O=0; z_O=0;
x_A=a; y_A=0; z_A=0;
x_B=0; y_B=b; z_B=0;
x_C=0; y_C=0; z_C=c;
x_D=a; y_D=b; z_D=0;
x_E=a; y_E=0; z_E=c;
x_F=a; y_F=b; z_F=c;
x_G=0; y_G=b; z_G=c;

hold on
axis([-20 20 -20 20 -20 20])
grid on

t1=text...
(0-.2, 0, 0+1.4,' O','fontsize',12);
t2=text...
(x_A, y_A-1.7, z_A-.2,' A','fontsize',12);
t3=text...
(x_B-.2, y_B, z_B-.1,' B','fontsize',12);
t4=text...
(x_C-1.1, y_C-1.1, z_C+.2,' C','fontsize',12);
t5=text...
(x_D+0.4, y_D+.2, z_D+.2,' D','fontsize',12);
t6=text...
(x_E, y_E-1.5, z_E+.2,' E','fontsize',12);
t7=text...
(x_F-1.1, y_F-1.3, z_F+.5,' F','fontsize',12);
t8=text...
(x_G, y_G, z_G+.7,' G','fontsize',12);

%input the prism vertices
vert = [x_O y_O z_O; x_A y_A z_A;...
        x_E y_E z_E; x_C y_C z_C;... 
        x_D y_D z_D; x_F y_F z_F;...
        x_G y_G z_G; x_B y_B z_B]; 
%input the faces of the prism
fac = [1 2 3 4; ... 
    2 1 8 5; ... 
    3 4 7 6; ... 
    4 1 8 7; ... 
    2 3 6 5; ... 
    5 6 7 8]; 
%draw the prism using patch function
prism=patch('Faces',fac,...
'Vertices',vert,'FaceColor','b');   

view(100,50); 
light('Position',[1 3 2]); 
%alpha sets one of three transparency properties
%depending on what arguments you specify 
alpha(prism,0.3); 

%scale factor 
s_f=0.4;

line([x_O,x_F],[y_O,y_F],[z_O,z_F]);
q=quiver3...
(x_F,y_F,z_F,s_f*x_F,s_f*y_F,s_f*z_F);
set(q,'Color','r');

str1 = {'\omega'};
text(x_F+7,y_F+5,16,str1,'fontsize',16);

xlabel('x'); ylabel('y');zlabel('z');

pause(3);
%delete letters A, B, C, D, E, G,
delete(t2);delete(t3);delete(t4);
delete(t5);delete(t6);delete(t8);
%delete the prism
delete(prism);

%rotate the prism about OF
for t=0:0.01:0.15
%call the rotation matrix 
%given an axis (omega axis)
%and an angle (angle of rotation)
P5 = Rotate(vert',t,a,b,c);
vert = P5'; 
% draw the prism
prism=patch...
('Faces',fac,'Vertices',vert,'FaceColor','b'); 
alpha(prism,0.3);
pause(0.7);
delete(prism);
end

%final position of the prism
prism=patch...
('Faces',fac,'Vertices',vert,'FaceColor','b'); 
alpha(prism,0.3); 
hold off

% end of program