%%   lportrait.m
%%   phase portraits for linear equations
%%   Requires linearde.m & darrow.m

mng=5;  % maximum number of integrations from each IC

a=input('a=');
b=input('b=');
c=input('c=');
d=input('d=');

ew=1;   % width of eigenvectors

%%  calculate eigenvalues and eigenvectors

[V,D]=eig([a b;c d])

[x,y]=meshgrid(-1:.25:1,-1:.25:1);
DxDt=a*x+b*y; DyDt=c*x+d*y;

hold on

%%  if the eigenvalues are real then draw the eigenvectors

if imag(D(1))==0

   x1=linspace(-1,1,101)*V(1,1);
   y1=linspace(-1,1,101)*V(2,1);
   plot(x1,y1,'LineWidth',ew)
   x2=linspace(-1,1,101)*V(1,2);
   y2=linspace(-1,1,101)*V(2,2);
   plot(x2,y2,'LineWidth',ew)

%% draw arrows on eigenvectors

   s=0.05;
   xo=.5*V(1,1); yo=.5*V(2,1); nx=D(1,1)*xo; ny=D(1,1)*yo;
   darrow(xo,yo,nx,ny,s)
   xo=-xo; yo=-yo; nx=-nx; ny=-ny;
   darrow(xo,yo,nx,ny,s)
   xo=.5*V(1,2); yo=.5*V(2,2); nx=D(2,2)*xo; ny=D(2,2)*yo;
   darrow(xo,yo,nx,ny,s)
   xo=-xo; yo=-yo; nx=-nx; ny=-ny;
   darrow(xo,yo,nx,ny,s)
   if D(1,1)*D(2,2)>0
      if abs(D(1,1))>abs(D(2,2))
         xo=.55*V(1,1);yo=.55*V(2,1); nx=D(1,1)*xo; ny=D(1,1)*yo;
         darrow(xo,yo,nx,ny,s)
         xo=-xo; yo=-yo; nx=-nx; ny=-ny;
         darrow(xo,yo,nx,ny,s)
      else
         xo=.55*V(1,2); yo=.55*V(2,2); nx=D(2,2)*xo; ny=D(2,2)*yo;
         darrow(xo,yo,nx,ny,s)
         xo=-xo; yo=-yo; nx=-nx; ny=-ny;
         darrow(xo,yo,nx,ny,s)
      end
   end

end

axis equal

%%   now draw trajectories

tspan=linspace(0,1,101);

X0=0;

%%   enter `999' for initial condition to stop

while X0~=999
   X0=input('Initial condition:   ');

%%   trajectories starting at both X0 and -XO are drawn,
%%   both forwards and backwards in time

   if X0~=999
      for yp=1:2;
          xm=(-1)^yp; ym=(-1)^yp;
          Z0=[xm ym].*X0;

          dir=linearde(0,Z0,[],a,b,c,d);

          xo=Z0(1); yo=Z0(2); nx=dir(1); ny=dir(2); s=.05;
          darrow(xo,yo,nx,ny,s)

          for p=1:2;
              m=(-1)^p;

              j=102; Y0=Z0; gn=0;

              while j==102 & gn<mng

                [t,u]=ode45('linearde',tspan,Y0,[],a*m,b*m,c*m,d*m);
                u(102,:)=[0 0];

                j=2;

                while abs(u(j,1))<=1 & abs(u(j,2))<=1 & ...
                (abs(u(j,1))+abs(u(j,2)))>=0.001 & j<=101
                    plot([u(j-1,1); u(j,1)],[u(j-1,2); u(j,2)])
                    j=j+1;
                end

                Y0=u(101,:); gn=gn+1;

              end

           end

       end

   end
end

axis off
axis equal
xlim([-1 1])
ylim([-1 1])
