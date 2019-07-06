%% Integrate the Lorenz equations
%% Requires lorenzde.m

tspan=[0 50];
yzero=input('yzero = ');
[t,y]=ode45('lorenzde',tspan,yzero);

%% Draw the trajectory in R^3

subplot 221
plot3(y(:,1),y(:,2),y(:,3))
xlabel('x'); ylabel('y'); zlabel('z');
title('Lorenz attractor')

%% And draw the components too

subplot 222
plot(t,y(:,1))
xlabel('t'); ylabel('x');
subplot 223
plot(t,y(:,2))
xlabel('t'); ylabel('y');
subplot 224
plot(t,y(:,3))
xlabel('t'); ylabel('z');
