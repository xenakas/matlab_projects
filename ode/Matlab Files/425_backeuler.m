%% backwards Euler method

T=8;        %% final time

h=0.5;      %% timestep

%% MATLAB does not allow an index 0
%% on a vector, so x_n is x(n+1) here

t(1)=0;     %% initial time
x(1)=.5;    %% initial condition

for n=1:T/h;

t(n+1)=n*h;

%% Use iterative method to find x(n+1) given x(n)

gn=x(n);
g=gn+2*h^3;

while abs(gn-g)>h^3;

%% method is O(h) so approximate x(n+1) to within O(h^3)

g=gn;
gn=x(n)+h*g*(1-g);

end

x(n+1)=gn;

end

%% Plot crosses at numerical values, and join these

plot(t,x,'x','MarkerSize',20)
hold on
plot(t,x)