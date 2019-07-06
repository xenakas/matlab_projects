%%  solvem.m
%%  solve d^3 copies of the Lorenz equations
%%  starting in a small ball
%%  Requires lorenzde.m

clf
clear

d=input('d = ');

for i=0:1999;
   tspan(i+1)=i*35/1999;
end

for i=0:d;
   for j=0:d;
      for k=0:d;
         m=1+k+(j*(d+1))+(i*(d+1)^2)
         y0(:,m)=[(i-4.5)/100000; 1+(j-4.5)/100000; (k-4.5)/100000];
         [t,y]=ode45('lorenzde',tspan,y0(:,m));
         yi(2000,3,m)=1;
         yi(:,:,m)=y;
      end
   end
end
