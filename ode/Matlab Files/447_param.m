%%  Solution to Exercise 4.1(v)

t=linspace(0,5);

hold on

for A=-3:3;
   for B=-3:3;
      x=B*exp(-t)+A.*t.*exp(-t);
      y=A*exp(-t);
      plot(x,y)
   end
end

hold off
