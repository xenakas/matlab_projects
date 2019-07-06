%%  lotkaplane.m
%%  draw phase plane diagram of the
%%  Lotka-Volterra model of two species
%%  Requires lotakde.m & darrow.m

hold on

%%  input parameters

A=input('A = ');
a=input('a = ');
b=input('b = ');
B=input('B = ');
c=input('c = ');
d=input('d = ');

X0=0;

%% mX is largest value of x and y in the phase diagram
%% change this number if the stationary points lie
%% outside 0<=x,y<=5

mX=5;

%% X0 is initial condition
%% if X0=999 then program finishes
%% if X0 has non-zero imaginary part then the integration
%%    continues from the end of the previous trajectory

while X0~=999

   X0=input('Initial condition:   ');
   if imag(X0)~=0
      X0=[u(2001,1) u(2001,2)];
   end

%% T is maximum time
%% if T is negative there is no arrow
%% if T has non-zero imaginary part the trajectory is
%%    drawn in bold

   T=input('Maximum time:   ');
   if X0~=999

      RX=X0;

%%  trajectory is calculated both forwards and backwards from X0

      for k = 1:2,

         X0=RX;

         sn=(-1)^k;

         rT=2*abs(real(T)); rN=4001; u(1,1)=mX*2; u(2,2)=mX*2;

%%  make sure that trajectory doesn't get too large

         while (max([max(u(:,1)) max(u(:,2))])>mX)

            rT=rT/2; rN=(rN-1)/2+1;

            tspan=linspace(0,rT,rN);

            [t,u]=ode45('lotkade',tspan,X0,[],sn,A,B,a,b,c,d);

         end

         if imag(T)==0
            plot(u(:,1), u(:,2))
         else
            plot(u(:,1), u(:,2), 'LineWidth',2)
         end


      end

%%  put the arrow on the trajectory (unless T is imaginary)

      D=lotkade(0,X0,[],1,A,B,a,b,c,d);

      if real(T)>0; darrow(X0(1), X0(2), D(1), D(2), 0.05); end

   end
end

hold off