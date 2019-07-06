% example 4.6

clear all; clc; close 

syms x y k g t a v0 real

% equations of motion
xs = dsolve('D2x = - k^2*x');
ys = dsolve('D2y = - k^2*y-g');

fprintf('x = %s  \n',char(xs))
fprintf('y = %s  \n',char(ys))
fprintf('\n')

dxs = diff(xs,t);
dys = diff(ys,t);
fprintf('dx/dt = %s  \n',char(dxs))
fprintf('dy/dt = %s  \n',char(dys))
fprintf('\n\n')

fprintf ...
('t=0,x(0)=a,y(0)=-g/k^2,dx(0)=0,dy(0)=v0=>\n')

eqx = a-subs(xs, t, 0);
eqy = -g/k^2-subs(ys, t, 0);
eqdx = 0-subs(dxs, t, 0);
eqdy = v0-subs(dys, t, 0);

fprintf(' %s = 0 \n',char(eqx))
fprintf(' %s = 0 \n',char(eqy))
fprintf(' %s = 0 \n',char(eqdx))
fprintf(' %s = 0 \n',char(eqdy))
fprintf('\n')

sol=solve(eqx,eqy,eqdx,eqdy,'C2,C3,C5,C6');
C2s = sol.C2;
C3s = sol.C3;
C5s = sol.C5;
C6s = sol.C6;
fprintf('C2 = %s \n',char(C2s))
fprintf('C3 = %s \n',char(C3s))
fprintf('C5 = %s \n',char(C5s))
fprintf('C6 = %s \n',char(C6s))
fprintf('\n\n')

xn = subs(xs,{'C2','C3'},{C2s, C3s});
yn = subs(ys,{'C5','C6'},{C5s, C6s});
fprintf('x = %s \n',char(xn))
fprintf('y = %s \n',char(yn))
fprintf('\n')

axis manual
axis equal
hold on
grid on
axis([-4.5 4.5 -4.5 4.5])

slist={a, g, v0, k};
nlist={2, 9.8, 10, 3};

xnt=subs(xn,slist,nlist);
ynt=subs(yn,slist,nlist);

x0 = subs(xnt,t,0);
y0 = subs(ynt,t,0);

xC = 0;
yC = subs(-g/k^2,slist,nlist);

text(0,0+.2,' O','fontsize',12)
text(x0,y0,'  P_O','fontsize',12)
text(xC,yC+.2,'  C','fontsize',12)
plot(xC,yC,'.','Color','black')

line ...
([0 0],[0 3],'LineWidth',1.2,'Color','red')
line ...
([0 3],[0 0],'LineWidth',1.2,'Color','red')

text(0,3,' y','fontsize',12)
text(3,0,' x','fontsize',12)
plot(0,0,'.','Color','black')
line ...
([xC x0],[yC y0],'LineStyle','--','LineWidth',1.5)
plot(xC,yC,'.','Color','red')


for t_n = 0 : 0.005 : 2.1
   x_t=subs(xnt,t,t_n);
   y_t=subs(ynt,t,t_n);
   hm=plot(x_t,y_t,'.','Color','red');
   pause(0.001)
   ht=plot(x_t,y_t,'.');
   delete(hm);
end
plot(x0,y0,'.','Color','red')

% end of program

% % Plotting parametrically
% close all; clear all;
% t = linspace(0,3*pi,100); % 50 equally spaced points between 0 and pi
% x = t.*cos(t);     % Again, the .* means multiply element by element
% y = t.*sin(t);     % rather than matrix or vector multiplication
% plot(x,y)
% axis equal % Make the scaling the same on both axes

