function yprime=lorenzde(t,y)
yprime=[10*(y(2)-y(1)); 28*y(1)-y(2)-y(1)*y(3); y(1)*y(2)-8*y(3)/3];
