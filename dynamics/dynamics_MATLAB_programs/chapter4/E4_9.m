% example 4.9

clear all; clc; close all 

syms m g alpha v_0 x t v_B y
% a)
xs=dsolve...
('m*D2x+m*g*sin(alpha)=0','x(0)=0','Dx(0)=v_0');
fprintf('eom IC: t=0=> x(0)=0, Dx(0)=v_0 \n')
fprintf('x = %s \n',char(xs))

vxs=dsolve...
('m*Dv+m*g*sin(alpha)=0','v(0)=v_0');
fprintf('vx = %s \n\n',char(vxs))

ts=solve(subs(x-xs,x,v_0^2/g),t);
ts=simplify(ts);
fprintf('t_AB1 = %s  \n',char(ts(1)))
fprintf('t_AB2 = %s  \n',char(ts(2)))
fprintf('\n')

v_B1=simplify(subs(vxs,t,ts(1)));
v_B2=simplify(subs(vxs,t,ts(2)));

fprintf('v_B1 = %s  \n',char(v_B1));
fprintf('v_B2 = %s  \n',char(v_B2));
fprintf('\n\n');

% numerical data
slist={m, g, v_0, alpha};
nlist={1,9.8,10,pi/10};

tn1=subs(ts(1),slist,nlist);
tn2=subs(ts(2),slist,nlist);
fprintf('t_AB1 = %6.3f (s) \n',tn1)
fprintf('t_AB2 = %6.3f (s) \n',tn2)
fprintf('\n')

v_B1n=subs(v_B1,slist,nlist);
v_B2n=subs(v_B2,slist,nlist);
fprintf('v_B1 = %6.3f (m/s) \n',v_B1n)
fprintf('v_B2 = %6.3f (m/s) \n',v_B2n)
 if v_B1n > 0 
    v_Bn = v_B1n; tn=tn1; vB=v_B1;
 else
    v_Bn = v_B2n; tn=tn2; vB=v_B2;
       end
fprintf('v_B = %6.3f (m/s) \n',v_Bn)
fprintf('t_AB = %6.3f (m/s) \n',tn)
fprintf('\n\n')

l=v_0^2/g;
xB=l*cos(alpha);
yB=l*sin(alpha);
xBn=subs(xB,slist,nlist)
yBn=subs(yB,slist,nlist);

% b)
N=m*g*cos(alpha);
Nn=subs(N,slist,nlist);
fprintf('N = %s = %6.3f (N)\n',char(N),Nn)
fprintf('\n\n')

fprintf('inclined plane \n\n');
% IC: t=0 => x(0)=x_0, y(0)=y_0
% Dx(0)=v_B*cos(alpha), Dy(0)=v_B*sin(alpha)

vx_p=dsolve('m*Dv=0','v(0)=v_B*cos(alpha)');
vy_p=dsolve('m*Dv+m*g=0','v(0)=v_B*sin(alpha)');
fprintf('v_x = %s  \n',char(vx_p))
fprintf('v_y = %s  \n',char(vy_p))
fprintf('\n\n')

x_p=dsolve...
('m*D2x=0','x(0)=0','Dx(0)=v_B*cos(alpha)');
y_p=dsolve...
('m*D2y+m*g=0','y(0)=0','Dy(0)=v_B*sin(alpha)');
fprintf('x = %s  \n',char(x_p))
fprintf('y = %s  \n',char(y_p))
fprintf('\n\n')

x_pn=simplify(subs(x_p,v_B,vB));
y_pn=simplify(subs(y_p,v_B,vB));
fprintf('x = %s  \n',char(x_pn))
fprintf('y = %s  \n',char(y_pn))
fprintf('\n\n')

t_x=solve(x-x_pn,t);
y_eq=subs(y_pn,t,t_x);
fprintf('trajectory of the particle\n')
fprintf('y = %s  \n',char(y_eq));
fprintf('\n\n');

% yC=-v0^2*sin(alpha)/g
eq=subs(y-y_eq,y,-v_0^2*sin(alpha)/g);

x_C=simplify(solve(eq,x));
fprintf('x_C1 = %s \n',char(x_C(1)))
fprintf('x_C2 = %s \n',char(x_C(2)))
fprintf('\n\n');

xCn=subs(x_C,slist,nlist);
 if xCn(1) > 0 
   xC = xCn(1); xCs = x_C(1);
 else
   xC = xCn(1); xCs = x_C(2);
       end
fprintf('xC = %6.3f (m) \n',xC)
fprintf('\n\n')

x_AC=v_0^2*cos(alpha)/g+xCs;
fprintf('x_AC = %s \n',char(x_AC))
fprintf('\n\n')

x_ACval=subs(x_AC,slist,nlist);
fprintf('x_AC = %6.3f (m) \n',x_ACval)
fprintf('\n\n')

axis manual
hold on
grid on
axis([-2 18 -2 5])

a=xBn; b=yBn;
x_O=0; y_O=0;
x_A=a; y_A=0;
x_B=a; y_B=b;

t1=text...
(x_O-0.6, y_O+0.1, 0,' O','fontsize',12);
t2=text(x_A, y_A, 0,' A','fontsize',12);
t3=text(x_B, y_B-0.1, 0,' B','fontsize',12);
t3=text(x_ACval, 0, 0,'  C','fontsize',12);

%input the inclined plane  vertices
vert = [x_O y_O 0; x_A y_A 0; x_B y_B 0;]; 
%input the faces inclined plane
fac = [ 2 1 3]; 
%draw the prism
prism=patch...
('Faces',fac,'Vertices',vert,'FaceColor','y'); 

for t_n = 0 : 0.005 : tn
   slist={m,g,v_0, alpha,'t'};
   nlist={2,9.8,10,pi/10,t_n};
   x_t=subs(xs*cos(alpha),slist,nlist);
   y_t=subs(xs*sin(alpha),slist,nlist);
   hm=plot(x_t,y_t,'.','Color','red');
   ht=plot(x_t,y_t,'.');
   pause(0.001)
   delete(hm);
end

for t_n = 0 : 0.005 : 1.02
   slist1={g,v_0, alpha,'v_B','t'};
   nlist1={9.8,10,pi/10,v_Bn,t_n};
   x_tn = xBn+subs(x_pn,slist1,nlist1);
   y_tn = yBn+subs(y_pn,slist1,nlist1);
   hm=plot(x_tn,y_tn,'.','Color','red');
   ht=plot(x_tn,y_tn,'.');
   pause(0.001)
   delete(hm);
end

% end of program
