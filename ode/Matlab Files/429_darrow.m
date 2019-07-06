function darrow(xo,yo,nx,ny,s)

%% Draw an arrow at (xo,yo) in the (nx,ny)
%% direction, of size s

if s==0;   s=0.4;   end

%% theta is (half) the angle the sides of the arrow make
%% to each other

theta=pi/9;

l=sqrt(nx^2+ny^2); if l~=0; nx=nx/l; ny=ny/l; end

if nx~=0
   ppx=-(ny/nx); ppy=1;
   nl=sqrt(ppx^2+ppy^2);
   ppx=ppx/nl; ppy=ppy/nl;
else
   ppx=1; ppy=0;
end

plot([xo-s*nx+s*tan(theta)*ppx xo],[yo-s*ny+s*tan(theta)*ppy yo])
plot([xo-s*nx-s*tan(theta)*ppx xo],[yo-s*ny-s*tan(theta)*ppy yo])
