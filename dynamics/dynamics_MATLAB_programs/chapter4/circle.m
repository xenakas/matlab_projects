function H=circle(x_c,y_c,radius)


axis equal
hold on
grid on

for theta = 0 : 0.005 : 2*pi 
x=x_c+radius*cos(theta);
y=y_c+radius*sin(theta);
plot(x,y);
end

