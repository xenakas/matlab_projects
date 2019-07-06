%% Plot iterates of logsitic map
%% against n

clf
hold on

r=input('r = ');
N=input('N = ');

xlim([0 N]); ylim([0 1]);
xlabel('n','FontSize',20); ylabel('x','FontSize',20)

x0=input('initial condition = ');

while x0>=0 & x0<=1

   x(1)=x0; t(1)=0;

   for n=2:N;

      t(n)=n-1;
      x(n)=r*x(n-1)*(1-x(n-1));

   end

   plot(t,x)
   plot(t,x,'x','MarkerSize',15)

   x0=input('initial condition = ');

end
