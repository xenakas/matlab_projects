function dx = eomE6_8(t,x)
dx = zeros(4,1);
dx(1) = x(2);
dx(2) = (2*(300*sin(3*t) - (x(2)*x(3)*x(4))/1000))/(x(3)^2/1000 + 40);
dx(3) = x(4);
dx(4) = (x(2)^2*x(3))/2 - (981*2^(1/2))/200;