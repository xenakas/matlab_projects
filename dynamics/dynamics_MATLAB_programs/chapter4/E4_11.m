% example 4.11

clear all; clc; close 

syms m g v_0 alpha t real 

% IC t=0 =>
% x(0)=0, y(0)=0 
% Dx(0)=v_0*cos(alpha), Dy(0)=v_0*sin(alpha)
xs=dsolve('m*D2x=0','x(0)=0','Dx(0)=v_0*cos(alpha)');
ys=dsolve('m*D2y=-m*g','y(0)=0','Dy(0)=v_0*sin(alpha)');
zs=dsolve('m*D2z=0','z(0)=0','Dz(0)=0');

fprintf('position \n')
fprintf('x = %s  \n',char(xs))
fprintf('y = %s  \n',char(ys))
fprintf('z = %s  \n',char(zs))

Dx=diff(xs,t);
Dy=diff(ys,t);
Dz=diff(zs,t);

fprintf('velocity \n');
fprintf('vx = %s  \n',char(Dx))
fprintf('vy = %s  \n',char(Dy))
fprintf('vz = %s  \n',char(Dz))
fprintf('\n')

syms x y real 
t_x = solve(x-xs,t);
t_y = solve(y-ys,t);

yn=subs(ys,t,t_x);
fprintf('trajectory \n')      
fprintf('y = \n')
pretty(yn)
fprintf('\n\n')

magn_v=subs(sqrt(Dx^2+Dy^2),t,t_y);
magn_v1=simplify(magn_v(1));
fprintf('|v|= %s \n',char(magn_v1))
fprintf('\n')

x_As=solve(yn,x);
x_A = simplify(x_As(2));
fprintf('x_A = %s \n',char(x_A))
fprintf('\n')

x_max=subs(x_A,sin(2*alpha),1);
fprintf('x_max = %s \n',char(x_max))
fprintf('\n')

x_B=simplify(solve(diff(yn,x),x));
fprintf('x_B = %s \n',char(x_B))
fprintf('\n')

y_B=simplify(subs(yn,x,x_B));
fprintf('y_B = %s \n',char(y_B))
fprintf('\n')

y_max=simplify(subs(y_B,sin(alpha),1));
fprintf('y_max = %s \n',char(y_max))
fprintf('\n')

slist={v_0, alpha, g};
nlist={10, pi/3, 9.8};

ynn=subs(yn,slist,nlist);
fprintf('y =  %s',char(ynn))
fprintf('\n')

axis manual
axis equal
hold on
grid on
axis([-1 10 -1 4])

for x_t = 0 : 0.01 : 8.83 
   slist={g,v_0, alpha,x};
   nlist={9.8,10,pi/6,x_t};
   y_t = subs(ynn,slist,nlist);
   hm=plot(x_t,y_t,'k.','Color','red');
   ht=plot(x_t,y_t);
   pause(0.001)
   delete(hm);
end

% end of program
