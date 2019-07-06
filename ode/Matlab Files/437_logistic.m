%%  logistic.m
%%  draw cobweb diagrams for logistic map

hold off
r=input('r=:  ');
x=0:0.01:1;
y=r.*x.*(1-x);
plot(x,x)
hold on
plot(x,y, 'LineWidth', 2)
z=input('x_0=:  ');
n=1; while n~=0
   n=input( ...
'# iterations (0=stop, .1=clear, .2=change IC, .3=change r):  ');
   if n==0.1
      hold off; plot(x,x); hold on; plot(x,y, 'LineWidth', 2);
   end
   if n==0.2
      z=input('new initial condition:  ');
   end
   if n==0.3
      r=input('new value of r:   ');
      hold off
      y=r.*x.*(1-x);
      plot(x,x)
      hold on
      plot(x,y)
      z=input('z_0=:  ');
   end
   if n~=0
      for j=1:abs(n)
         nz=r*z*(1-z);
         if n>0
            plot([z z],[z nz])
            ms=min(8,25*abs(nz-z));
            if z<nz plot(z,(z+nz)/2,'^','markersize',ms)
            else
               plot(z,(z+nz)/2,'v','markersize',ms)
            end
         plot([z nz],[nz nz])
            if z<nz plot((z+nz)/2,nz,'>','markersize',ms)
            else
               plot((z+nz)/2,nz,'<','markersize',ms)
            end
      	end
      z=nz;
   	end
	end
end
