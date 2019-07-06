%%  newtonplane.m
%%  draw phase plane diagrams in whole plane
%%  requires newtonde.m & darrow.m

clf
hold on

dmp=input('damping k = ');

mX=5;

X0=0;

while X0~=999
   X0=input('Initial condition:   ');
   if imag(X0)~=0
      X0=[u(2001,1) u(2001,2)];
   end
   T=input('Maximum time:   ');
   if X0~=999
 %     tspan=linspace(0,T,2001);

      RX=X0;

      for k = 1:2,

         X0=RX;

         sn=(-1)^k;

         rT=2*abs(real(T)); rN=4001; u(1,1)=mX*2; u(2,2)=mX*2;

         while (max([max(u(:,1)) max(u(:,2))])>mX)

            rT=rT/2; rN=(rN-1)/2+1;

            tspan=linspace(0,rT,rN);

            [t,u]=ode45('newtonde',tspan,X0,[],sn,dmp);

         end

            if imag(T)==0
             plot(u(:,1), u(:,2))
          else
             plot(u(:,1), u(:,2), 'linewidth', 2)
          end


      end

      D=newtonde(0,X0,[],1,dmp);

      if real(T)>0; darrow(X0(1), X0(2), D(1), D(2), .1); end


   end
end
