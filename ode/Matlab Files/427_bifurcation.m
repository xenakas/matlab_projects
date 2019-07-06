%% bifurcation.m
%% draw bifurcation diagram for logistic map

clf
hold on

r1=input('Smallest value of r: ');
r2=input('Largest value of r: ');

%% N is the number of r values used
%% lowering/increasing N will decrease/increase the resolution
%% with a consequent decrease/increase in the time needed

N=400;

for t=0:N

%% for each value of r start at x_0=0.5

   z=0.5;

   r=r1+((r2-r1)*t/N);

%% compute one hundred iterates of f

   for j=1:100
      z=r*z*(1-z);
   end

%% and then plot the next 30

   for j=1:30
      z=r*z*(1-z);
      plot(r,z,'.','MarkerSize',2)
   end

end

xlabel('r','FontSize',20)

xlim([r1 r2])
