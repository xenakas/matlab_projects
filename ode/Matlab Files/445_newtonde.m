function yprime=newtonde(t,y,flag,s,k)

%% Comment out unwanted equations

%% Particles moving in the following potentials:

%%      V(x) = x-x^3/3
 yprime=s.*[y(2); -k*y(2)-1+y(1)^2];

%%      V(x) = x^2/2
%% yprime=s.*[y(2); -k*y(2)-y(1)];

%%      V(x) = x^4/2-x^2
%% yprime=s.*[y(2); -k*y(2)-2*y(1)^3+2*y(1)];

%%      V(x) = x^6/6 - 5x^4/4 + 2x^2
%% yprime=s.*[y(2); -k*y(2)-y(1)^5+5*y(1)^3-4*y(1)];

%% Particles moving on a wire in the following shapes

%%      V(x) = x-x^3/3
%% yprime=s.*[y(2); -k*y(2)-(1-y(1)^2)*(1-2*y(1)*y(2)^2)/(2-2*y(1)^2+y(1)^4)];

%%      V(x) = x^2/2
%% yprime=s.*[y(2); -k*y(2)-y(1)*(1+y(2)^2)/(1+y(1)^2)];

%%      V(x) = x^3/2-x^2
%% n=2*y(1)-2*y(1)^3-(2*y(1)^3-y(1))*(6*y(1)^2-1)*y(2)^2;
%% yprime=s.*[y(2); -k*y(2)+n/(1+(2*y(1)^3-y(1))^2)];
