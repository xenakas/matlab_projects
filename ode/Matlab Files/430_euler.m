%% Euler's method

T=12;        %% final time

h=0.5;       %% timestep

%% MATLAB does not allow an index 0
%% on a vector, so x_n is x(n+1) here

t(1)=0;     %% initial time
x(1)=0;     %% initial condition

for n=1:T/h;

    t(n+1)=n*h;
    x(n+1)=x(n) + h * (t(n) - x(n)^2);

end

[t x]       %% display values

%% Plot crosses at numerical values, and join these

plot(t,x,'x','MarkerSize',20)
hold on
plot(t,x)
