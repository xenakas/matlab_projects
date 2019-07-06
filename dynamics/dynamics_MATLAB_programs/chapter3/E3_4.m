% example 3.4
clear all; clc; close all

syms t1 t2 t3 V a1 a3 x t real
syms c1 c2 c3 c4 c5  real

list={t1, t2, t3, V};
listn={1, 1, 3, 1}; 
fprintf('segment OA\n'); fprintf('\n')
ddx1 = a1;
fprintf('ddx1 = a1 \n')
fprintf('\n');
dx1 = int(ddx1, t) + c1;
fprintf('dx1 = %s \n',char(dx1))
dx10 = subs(dx1, t, 0);
dx11 = subs(dx1, t, t1)-V;
sol1 = solve(dx10,dx11, 'c1, a1');
c10 = sol1.c1;
a10 = sol1.a1;
fprintf('t=0 ; dx1=0 => c1 = %s \n',char(c10))
fprintf('t=t1; dx1=V => a1 = %s \n',char(a10))
fprintf('\n')
a1s = a10;
a1n = subs(a1s, list, listn);
fprintf('a1 = %s = %s (m/s^2)\n',char(a1s),char(a1n))
fprintf('\n')
v1 = int(a10, t) + c10;
v1n = subs(v1, list, listn);
fprintf('v1 = %s = %s (m/s)\n',char(v1),char(v1n))
fprintf('\n')
x1s = int(v1, t) + c2;
fprintf('x1 = %s \n',char(x1s))
x10 = subs(x1s, t, 0);
c20 = solve(x10,'c2');
fprintf('t=0; x1=0 => c2 = %s \n', char(c20))
x1 = int(v1, t) + c20;
x1n = subs(x1, list, listn);
fprintf('x1 = %s = %s (m)\n',char(x1), char(x1n))
fprintf('\n')
d1 = subs(x1, t, t1);
d1n = subs(x1n, t, 1);
fprintf('d1 = %s = %g (m)\n',char(d1), d1n)
fprintf('\n')
fprintf('segment AB\n'); fprintf('\n')
dx2 = V;
fprintf('dx2 = v2 = V \n')
fprintf('\n')
ddx2 = diff(dx2, t);
fprintf('ddx2 = a2 = %s \n',char(ddx2))
fprintf('\n')
x2s = int(dx2, t) + c3;
fprintf('x2 = %s \n',char(x2s));
x21 = subs(x2s, t, t1)-d1;
c30 = solve(x21,'c3');
fprintf('t=t1; x2=d1 => c3 = %s \n', char(c30))
x2 = int(dx2, t) + c30;
x2n = subs(x2, list, listn);
fprintf('x2 = %s = %s (m)\n',char(x2), char(x2n))
fprintf('\n')
s2 = subs(x2, t, t1+t2);
s2n = double(subs(s2, list, listn));
fprintf('s2 = d1 + d2 = %s = %g (m)\n',char(s2),s2n)
fprintf('\n')
fprintf('segment CD\n'); fprintf('\n')
ddx3 = -a3;
fprintf('ddx3 = -a3 \n');
fprintf('\n');
dx3 = int(ddx3, t) + c4;
fprintf('dx3 = %s \n',char(dx3));
dx32 = subs(dx3, t, t1+t2)-V;
dx33 = subs(dx3, t, t1+t2+t3);
sol3 = solve(dx32,dx33, 'c4, a3');
c40 = sol3.c4;
a30 = sol3.a3;
c4n = double(subs(c40, list, listn)); 
a3n = double(subs(a30, list, listn)); 
fprintf('t=t1+t2   ; dx3=V \n')
fprintf('t=t1+t2+t3; dx3=0 => \n')
fprintf('c4 = %s = %g \n',char(c40),c4n)
fprintf('a3 = %s = %g \n',char(a30),a3n)
fprintf('\n')
a3s = -a30;
v3 = int(-a30, t) + c40;
v3n = subs(v3, list, listn);
fprintf('v3 = %s = %s (m/s)\n',char(v3),char(v3n))
fprintf('\n')
x3s = int(v3, t) + c5;
fprintf('x3 = %s \n',char(x3s));
x12 = subs(x3s, t, t1+t2)-s2;
c50 = solve(x12,'c5');
fprintf('t=t1+t2; x3=s2 =>  \n')
fprintf('c5 = %s \n', char(c50))
c50n = subs(c50, list, listn);
fprintf('c5 = %g \n', double(c50n))
fprintf('\n')
x3 = int(v3, t) + c50;
x3n = subs(x3, list, listn);
fprintf('x3 = %s \n',char(x3))
fprintf('x3 = %s (m)\n',char(x3n))
fprintf('\n')
s3 = subs(x3, t, t1+t2+t3);
s3n = subs(s3, list, listn);
fprintf('s3 = %s \n',char(s3))
fprintf('s3 = %g (m)\n',double(s3n))
fprintf('\n')
% Graphic
y1=0:.01:1;
y2=1:.01:2;
y3=2:.01:5;
y1= 0:.01:1;
y2= 1:.01:2;
y3= 2:.01:5;
px1=subs(x1n,t,y1);
px2=subs(x2n,t,y2);
px3=subs(x3n,t,y3);
pv1=subs(v1n,t,y1);
pv2=double(subs(V,list,listn));
pv3=subs(v3n,t,y3);
pa1=subs(a1n,t,y1);
pa2=0;
pa3=subs(a3n,t,y3);
subplot(3,1,1),...
plot(y1,px1,'b',y2,px2,'k',y3,px3,'r'),...
ylabel('x (m)'), grid,...
subplot(3,1,2),...
plot(y1,pv1,'b',y2,pv2,'k-',y3,pv3,'r'),...
ylabel('v (m/s)'), grid, axis([0 5 0 1.5]) 
subplot(3,1,3),...
plot(y1,pa1,'b',y2,pa2,'k',y3,pa3,'r') 
xlabel('t (s)'), ylabel('a (m/s^2)'), grid,...
axis([0 5 -1 2]) 

% end of program