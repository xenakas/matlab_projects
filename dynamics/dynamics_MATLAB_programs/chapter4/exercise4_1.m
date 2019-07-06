% exercise 4.1
clear all; clc; close all

syms theta R r v P Q g t  
h = sym('h(t)');

theta=h/R;
h1=theta*r;
omega=v/R;
v1=omega*r;

T=(1/2)*(P/g)*v^2+(1/2)*(Q/g)*v1^2;

fprintf('kinetic energy of the system \n')
fprintf('T = %s  \n',char(T))

U=P*h-Q*h1;
fprintf('total work \n')
fprintf('U = %s  \n',char(U))

dTU=simplify(T-U);

v=simple(solve(dTU,v));

fprintf('velocity v^2 = \n')
pretty(v(1)^2)

% v = dh/dt
% a = dv/dt
% d(v^2)/dt = 2 v a = 2 (dh/dt) a 

a=simple(diff(v(1)^2,t)/(2*diff(h,t)));

fprintf('acceleration a = \n')
pretty(a)

% MATLAB figure of the system
axis equal
hold on
axis([-70 70 -150 70])
R=30;
x_R=0;
y_R=0;
circle_m(x_R,y_R,R);
r=15;
x_r=0;
y_r=0;
circle_m(x_r,y_r,r);

x_O=0; y_O=0; z_O=0;

disp=1.5;
t1=text(x_O+disp, y_O-2*disp, z_O,'O','fontsize',8);

len_2=60;
len_1=100;
line([R R],[-len_2 0],'Marker','.',...
    'LineStyle','-','Color',[.8 .6 .8])
line([-r -r],[-len_1 0],'Marker','.',...
    'LineStyle','-','Color',[.8 .2 .8])

l=10;
rectangle('Position',[-r-l/2,-len_1-l/2,l,l],...
         'LineWidth',2,'LineStyle','-');
rectangle('Position',[R-l/2,-len_2-l/2,l,l],...
         'LineWidth',2,'LineStyle','-');
Q=20;
quiver(-r,-len_1,0,-Q)
P=30;
quiver(R,-len_2,0,-P)

% end of program