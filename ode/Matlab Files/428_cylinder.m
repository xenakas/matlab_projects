%%  Animate trajectories of the nonlinear pendulum
%%  on the `phase cylinder'

%%  Requires pendde.m

hold off

%% Draw portion of cylinder

[tpx,tpz]=meshgrid(-pi:.1:pi+.1,-3:.1:3);
tpcx=cos(tpx);
tpcy=sin(tpx);
surfl(tpcx,tpcy,tpz), shading flat

hold on

axis off

%% * is pendulum pointing down
%% x is pendulum pointing up

plot3(-1,0,0,'*')
plot3(1,0,0,'x')

X0=0;

%% program ends when initial condition is 999

while X0~=999

   X0=input('Initial condition:   ');
   if imag(X0)~=0
      X0=[u(2001,1) u(2001,2)];
   end
   if X0~=999
   	T=input('Maximum time:   ');
   	if X0~=999

	      tspan=linspace(0,T,2001);

	      [t,u]=ode45('pendde',tspan,X0,[],1);

	      c1(:,1)=cos(pi+u(:,1));
   	   c2(:,2)=sin(pi+u(:,1));
      	c3(:,3)=u(:,2);

			comet3(c1(:,1),c2(:,2),c3(:,3))

		end

   end

end

hold off
