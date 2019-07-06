%% Runge-Kutta scheme

T=12;
h=0.5;

t(1)=0;     %% initial time
x(1)=0;     %% initial condition

for n=1:T/h;

   t(n+1)=n*h;

   f1=t(n)-x(n)^2;
   f2=(t(n)+(h/2))-(x(n)+(h*f1/2))^2;
   f3=(t(n)+(h/2))-(x(n)+(h*f2/2))^2;
   f4=(t(n)+h)-(x(n)+h*f3)^2;
   x(n+1)=x(n) + h * (f1+(2*f2)+(2*f3)+f4)/6;

end

[t x]       %% display values

%% Plot crosses at numerical values, and join these

plot(t,x,'x','MarkerSize',20)
hold on
plot(t,x)
