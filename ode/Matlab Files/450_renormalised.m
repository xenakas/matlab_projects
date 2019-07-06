%% Draw cobweb diagrams for renormalised version of logistic map

hold off
r=input('r=:  ');
plot([0 1],[0 1])

nxl=1-(1/r);
g=inline('(r^2*x*(1-x).*(1-(r*x*(1-x))))-1+(1/r)','x','r');
xl=1-(1/r)+0.00001; xr=1;
d=1;
while d>0.0001
   m=(xl+xr)/2;
   if g(xl,r)*g(m,r)<0
      xr=m;
   else
      xl=m;
   end
   d=(xr-xl);
end
nxr=m;

xp=linspace(0,1,100);
xv=nxl+xp*(nxr-nxl);
yv=r.*xv.*(1-xv); yv=r.*yv.*(1-yv);
yp=(yv-nxl)/(nxr-nxl);

hold on

plot(xp,yp,'linewidth',2)

axis equal; xlim([0 1]); ylim([0 1]);


z=input('z_0=:  ');
n=1;
while n~=0
   n=input('number of iterations (0 to stop, .1 to clear, .2 to change IC):  ');
   if n==0.1
      hold off; plot(x,x); hold on; plot(xp,yp, 'LineWidth', 2); axis equal;
      xlim([0 1]); ylim([0 1])
   end
   if n==0.2
      z=input('new initial condition:  ');
   end
   if n~=0
      for j=1:abs(n)
         za=nxl+z*(nxr-nxl);
         nza=r*za*(1-za); nza=r*nza*(1-nza);
         nz=(nza-nxl)/(nxr-nxl);
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
