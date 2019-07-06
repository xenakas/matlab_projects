% example 3.5
clear all; clc; close all

syms R omega t theta l_1 l_2 

theta = omega*t;
x = l_1*cos(theta);
y = l_2*sin(theta);

v_x = diff(x,t);
v_y = diff(y,t);
v = [v_x v_y];
magn_v=sqrt(v(1)^2+v(2)^2);

a_x = diff(x,t,2);
a_y = diff(y,t,2);
a = [a_x a_y];
magn_a = sqrt(a(1)^2+a(2)^2);

rho = (v_x^2+v_y^2)^(3/2)...
    /abs(v_x*a_y-a_x*v_y);

% numerical values
omegan = pi/2;
l_1n = 3;
l_2n = 2;

lists = {l_1, l_2, omega};
listt = {l_1n, l_2n, omegan};

xt = subs(x, lists, listt);
yt = subs(y, lists, listt);
  
vxt = subs(v_x, lists, listt);
vyt = subs(v_y, lists, listt);
vt = [vxt vyt];
magn_vt=sqrt(vt(1)^2+vt(2)^2);

axt = subs(a_x, lists, listt);
ayt = subs(a_y, lists, listt);
at = [axt ayt];
magn_at=sqrt(at(1)^2+at(2)^2);

fprintf('v = [v_x, v_y] = \n')
pretty(v); fprintf('\n\n')
pretty(vt); fprintf('\n\n')
fprintf('|v| = \n')
pretty(magn_v); fprintf('\n\n')
pretty(magn_vt); fprintf('\n\n')

fprintf('a = [a_x, a_y] = \n')
pretty(a); fprintf('\n\n')
pretty(at); fprintf('\n\n')
fprintf('|a| = \n')
pretty(simplify(magn_a)); fprintf('\n\n')
pretty(magn_at); fprintf('\n\n')

fprintf('rho = \n')
pretty(simplify(rho))
rhot = subs(rho, lists, listt);
pretty(simplify(rhot))



scale = 3;
axis manual
axis equal
hold on
grid on
axis([-5 5 -5 5])

syms data
data=[];
i=0;

for tn = 0 : 0.01 : 3*pi/2  
   
   xn = subs(xt, t, tn);
   yn = subs(yt, t, tn);
   
   hm=plot(xn, yn,'k.','Color','red');
   ht=plot(xn, yn);
   title('Trajectory, velocity and acceleration')
   pause(0.001)
   delete(hm);
   
   if tn==0 | tn==0.5 | tn==1 | tn==2 | tn==3 
        pause(0.2)
        vxn = subs(vxt, t, tn);
        vyn = subs(vyt, t, tn);
quiver(xn,yn,vxn/scale,vyn/scale, ...
'color','k','LineWidth',1.3);
        pause(0.7)
        axn = subs(axt, t, tn);
        ayn = subs(ayt, t, tn);
quiver(xn,yn,axn/scale,ayn/scale, ...
'color','r','LineWidth',1.3);

rhon = subs(rhot, t, tn);

        i=i+1;
data{i} = ...
{xn,yn,vxn,vyn,sqrt(vxn^2+vyn^2), ...
 axn,ayn,sqrt(axn^2+ayn^2),rhon};
   end
end

[legend_h, object_h, plot_h,text_strings]=...
 legend('Trajectory','Velocity','Acceleration');
  
set(plot_h(1),'color','b');
set(object_h(1),'color','b');
set(plot_h(2),'color','k');
set(object_h(2),'color','k');
set(plot_h(3),'color','r');
set(object_h(3),'color','r');

% create a table with data
f = figure('Position',[280 300 720 150]);
dat = [data{1};data{2};data{3}];
cnames = ...
{'x','y','v_x','v_y',' v ','a_x','a_y',' a ','rho'};
cformat = {'bank','numeric','numeric','bank',...
    'numeric','bank','numeric','numeric'};
rnames = ...
{'t_0=0','t_1=0.5','t_2=1'};
t = uitable('Parent',f,'Data',dat,'ColumnName',...
cnames,'ColumnFormat',cformat,'RowName',rnames,...
'Position',[20 20 674 100]);
   
% end of program