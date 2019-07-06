%%  piphase.m
%%  phase portrait for the nonlinear pendulum
%%  on -pi<=x<=pi

hold off

[x,y]=meshgrid(-pi:.5:pi,-3:.5:3);
DxDt=y; DyDt=-sin(x);
quiver(x,y,DxDt,DyDt)
hold on

plot([-pi -pi],[-4 4],'--')
plot([pi pi],[-4 4],'--')

X0=0;

while X0~=999
   X0=input('Initial condition [999 to end]:   ');
   if imag(X0)~=0
      X0=[u(2001,1) u(2001,2)];
   end
   if X0~=999
      T=input('Maximum time:   ');
      
      tspan=linspace(0,T,2001);

      [t,u]=ode45('pendde',tspan,X0,[],1);

      js=1; je=1;

      while je<2001
         je=je+1;
         if u(je,1)>pi
            comet(u(js:je,1),u(js:je,2))
            for j=je:2001
               u(j,1)=u(j,1)-2*pi;
            end
            js=je;
         end
         if u(je,1)<-pi
            comet(u(js:je,1),u(js:je,2))
            for j=je:2001
               u(j,1)=u(j,1)+2*pi;
            end
            js=je;
         end

      end

      if js~=2001;
         comet(u(js:je,1),u(js:je,2))
      end

   end

end

hold off
