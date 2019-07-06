% driver.m
% driver link 
% position, velocity, and acceleration

function out = driver(xA, yA, AB, phi, omega, alpha);
xB = xA + AB*cos(phi);
yB = yA + AB*sin(phi);

vBx = - AB*omega*sin(phi);
vBy =   AB*omega*cos(phi);

aBx = - AB*omega^2*cos(phi) - AB*alpha*sin(phi);
aBy = - AB*omega^2*sin(phi) + AB*alpha*cos(phi);

out = [xB, yB, vBx, vBy, aBx, aBy];

end
