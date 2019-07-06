function H=circle_m(x_c,y_c,radius)
grid on
nr_points=1000;
theta = linspace(0,2*pi,nr_points);
x = x_c + radius*cos(theta);
y = y_c + radius*sin(theta);
plot(x,y);