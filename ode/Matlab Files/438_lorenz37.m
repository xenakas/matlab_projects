%%  Lorenz `presentation'
%%  Requires lorenzde.m

clf
echo

% parameters

sigma=10; r=28; b=8/3;

echo off 
p=input('To continue press [RETURN]'); echo

% non-zero fixed points

x=sqrt(b*(r-1)); y=x; z=r-1;

echo off 
p=input('draw the fixed points (press [RETURN])');

% draw the three fixed points & some axes

plot3([0 0],[0 0],[-40 40],'g')
hold on
plot3([0 0],[-40 40],[0 0])
plot3([-40 40],[0 0],[0 0])
plot3(0,0,0,'x','MarkerSize',15)
plot3(x,y,z,'x','MarkerSize',15)
plot3(-x,-y,-z,'x','MarkerSize',15)
xlabel('x'); ylabel('y'); zlabel('z');

p=input('To continue press [RETURN]');
echo

% eigenvalues and eigenvectors at the origin
% L gives values, V gives vectors

[V L]=eig([-sigma sigma 0; r -1 0; 0 0 -b])

% two stable directions and one unstable direction

echo off
p=input('look close to origin (press [RETURN])');
echo

% draw local direction field near origin
% the picture is nice if you rotate to -141 30

[x0,y0,z0]=meshgrid(-2:1:2,-2:1:2,-2:1:2);
xd=sigma.*(-x0+y0);
yd=r.*x0-y0-x0.*z0;
zd=-b.*z0+x0.*y0;
quiver3(x0,y0,z0,xd,yd,zd)
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
zlim([-2.5 2.5])

a1=3.5*V(:,1);
a2=3.5*V(:,2);
a3=3.5*V(:,3);
plot3([-a1(1) a1(1)],[-a1(2) a1(2)],[-a1(3) a1(3)],'g')
plot3([-a2(1) a2(1)],[-a2(2) a2(2)],[-a2(3) a2(3)],'r')
plot3([-a3(1) a3(1)],[-a3(2) a3(2)],[-a3(3) a3(3)],'w')
plot3([-a3(1) a3(1)],[-a3(2) a3(2)],[-a3(3) a3(3)],'g')

echo off
p=input('To continue press [RETURN]');
echo

% eigenvalues and eigenvectors at the non-zero points

[V L]=eig([-sigma sigma 0; r-z -1 -x; y x -b])

% a stable direction and a `2d unstable focus'

echo off
p=input('look close to one of these fixed points (press [RETURN])');
echo


% draw local direction field near one of these fixed points

[x0,y0,z0]=meshgrid(x-2:1:x+2,y-2:1:y+2,z-2:1:z+2);
xd=sigma.*(-x0+y0);
yd=r.*x0-y0-x0.*z0;
zd=-b.*z0+x0.*y0;
quiver3(x0,y0,z0,xd,yd,zd)
xlim([x-3 x+3]);
ylim([y-3 y+3]);
zlim([z-3 z+3])

% draw some local axes to make picture easier to see

a1=3.5*V(:,1); a2=3.5*real(V(:,2)); a3=3.5*imag(V(:,2));
plot3([x-a1(1) x+a1(1)],[y-a1(2) y+a1(2)],[z-a1(3) z+a1(3)],'g')
plot3([x-a2(1) x+a2(1)],[y-a2(2) y+a2(2)],[z-a2(3) z+a2(3)],'r')
plot3([x-a3(1) x+a3(1)],[y-a3(2) y+a3(2)],[z-a3(3) z+a3(3)],'r')

% rotate to -85 44 to see `rotational part'

echo off

% draw local direction field near the other one of these fixed points

[x0,y0,z0]=meshgrid(-x-2:1:-x+2,-y-2:1:-y+2,z-2:1:z+2);
xd=sigma.*(-x0+y0);
yd=r.*x0-y0-x0.*z0;
zd=-b.*z0+x0.*y0;
quiver3(x0,y0,z0,xd,yd,zd)

% draw some local axes to make picture easier to see

plot3([-x-a1(1) -x+a1(1)],[-y-a1(2) -y+a1(2)],[z-a1(3) z+a1(3)],'g')
plot3([-x-a2(1) -x+a2(1)],[-y-a2(2) -y+a2(2)],[z-a2(3) z+a2(3)],'r')
plot3([-x-a3(1) -x+a3(1)],[-y-a3(2) -y+a3(2)],[z-a3(3) z+a3(3)],'r')

p=input('look at whole picture again (press [RETURN])');

echo

% and "join the dots"?
% show whole picture

xlim([-35 35]); ylim([-35 35]); zlim([-2 35])

% oh dear -- try looking at components of solutions
% as functions of time, and also the `bounding function'
% V(x,y,z)=x^2+y^2+(z-sigma-r)^2

echo off
p=input('To continue press [RETURN]');

figure(2)

tt=1; cla

yzero=[0;1;0];

while yzero~=999

   if tt==1; echo; end

   tspan=[0 50];    % length of time
   yzero=input('initial condition (input 999 to continue):    ');
   if yzero~=999
      [t,y]=ode45('lorenzde',tspan,yzero);
      figure(2)
      subplot 411
      plot(t,y(:,1)); ylabel('x')
      subplot 412
      plot(t,y(:,2)); ylabel('y')
      subplot 413
      plot(t,y(:,3)); ylabel('z')
      V=sqrt(y(:,1).^2+y(:,2).^2+(y(:,3)-sigma-r).^2);
      subplot 414
      plot(t,V); ylabel('V(x,y,z)^{1/2}')
   end

   if tt==1; echo off; end
   tt=tt+1;


end

echo

% this is no good -- look at trajectories moving in R^3

echo off

figure(1)

yzero=[0;1;0];

while yzero~=999

   if tt==1; echo; end

   tspan=[0 50];    % length of time
   yzero=input('initial condition (input 999 to continue):    ');
   if yzero~=999
      [t,y]=ode45('lorenzde',tspan,yzero);
      comet3(y(:,1),y(:,2),y(:,3))
   end

   if tt==1; echo off; end
   tt=tt+1;

end

echo

% still looks complicated... try plotting successive maximum
% values of z against each other...

echo off

p=input('To continue press [RETURN]');


yzero=[0 1 0]
figure
hold on

for k=1:25

tspan=[0 20]; [t,y]=ode45('lorenzde',tspan,yzero);

j=0;
for i=2:length(y)-1
   if y(i,3)>y(i+1,3) & y(i,3)>y(i-1,3)
      j=j+1;
      m(j)=y(i,3);
   end
end

for i=1:j-1
   plot(m(i),m(i+1),'x')
end

yzero=y(940,:);

end