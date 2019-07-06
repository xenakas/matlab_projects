% exercise 4.2

clear all; clc; close 

syms x k xA xB

VA=dsolve('DF = - k*xA','F(0) = 0','xA');
VB=dsolve('DF = - k*xB','F(0) = 0','xB');
fprintf('Potential energy \n')
fprintf('VA = %s  \n',char(VA))
fprintf('VB = %s  \n',char(VB))

U=abs(-(VB-VA));
fprintf('Work done by the particle \n')
fprintf('UAB=-(VB-VA) = %s  \n',char(U))

% numerical application
disp=4;
x_O1=0; % mm
x_O=20; % mm
x_A=40; % mm
x_B=30; % mm

subplot(2,2,1); 
axis manual
axis equal
hold on
grid on
axis([0 70 -10 10])

x_ii = pi:pi/2:x_B;
x_if =  2*sin(x_ii);

plot(x_ii,x_if,'-','LineWidth',2)
title 'Final spring position xB [mm]'
t1=text(x_O1,-disp,'O1','fontsize',8);
t2=text(x_O,-disp, 'O','fontsize',8);
t3=text(x_B,-disp,'B','fontsize',8);

subplot(2,2,3); 
axis manual
axis equal
hold on
grid on
axis([0 70 -10 10])

u_fi = pi:pi/2:x_A;
u_ff =  2*sin(u_fi);
plot(u_fi,u_ff,'-','Color','r','LineWidth',2)
title 'Initial spring position xA [mm]'
t1=text(x_O1,-disp,'O1','fontsize',8);
t2=text(x_O,-disp, 'O','fontsize',8);
t3=text(x_A,-disp,'A','fontsize',8);

subplot(2,2,[2,4]);
kn=1000; % N/m
slist={k,'xA','xB'};
nlist={kn,x_A-x_O,x_B-x_O};
Un = subs(U,slist,nlist);

bar(0,Un*10^(-3),0.5,'r');
title 'work [J]'

