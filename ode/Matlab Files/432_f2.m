%%  Renormalisation of logistic map f

r=input('value of r:  ');

clf reset

%%  The left-hand plot shows f and f^2

subplot(121)

x0=0:0.025:1;

%% y=f(x0) and z=f^2(x0)

y=r.*x0.*(1-x0);
z=r.*y.*(1-y);
plot(x0,y)
hold on
plot(x0,x0)
plot(x0,z, 'LineWidth', 2)
axis equal

%% Now find where f^2(nxl)=1-(1/r)
%% using the interval halving method

nxl=1-(1/r);
g=inline('(r^2*x*(1-x)*(1-(r*x*(1-x))))-1+(1/r)','x','r');
xl=1-(1/r)+0.00001; xr=1;
d=1;
while d>0.0001
   m=(xl+xr)/2;
   if g(xl,r)*g(m,r)<0
      xr=m;
   else
      xl=m;
   end
   d=(xr-xl);
end
nxr=m;

%% Draw the little `invariant box'

plot([nxr nxr],[nxr nxl])
plot([nxr nxl],[nxl nxl])
plot([nxl nxr],[nxr nxr])
plot([nxl nxl],[nxr nxl])

xlim([0 1])
ylim([0 1])

%%  On the right-hand side draw an expanded version of little box

subplot(122)

x2=linspace(nxl,nxr,100);
y2=r.*x2.*(1-x2); y2=r.*y2.*(1-y2);
plot(x2,y2, 'LineWidth', 2);
hold on
plot(x2,x2);
axis equal
xlim([nxl nxr])
ylim([nxl nxr])
