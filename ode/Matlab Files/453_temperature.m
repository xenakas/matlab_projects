%% Calculate temperature when given by solution of Exercise 9.4
%% input is k, t1 (initial time), and t2 (final time), and T(t1)

m=3; a=5; f=2; w=pi/12;

% k=0.3640; % this is value from Exercise 9.5(ii)
k=input('k = ');

t1=input('t1 = ');
t2=input('t2 = ');
Tt1=input('T(t1) = ');

u=k^2+w^2;

T=m+a*k*(k*cos(w*(t2-f))+w*sin(w*(t2-f)))/u;
T=T+((Tt1-m-(a*k*(k*cos(w*(t1-f))+w*sin(w*(t1-f)))/u))*exp(-k*(t2-t1)))