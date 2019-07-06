% linva.m
% linear velocity and acceleration

function out=...
 linva(xM,yM,xN,yN,vMx,vMy,aMx,aMy,dtheta,ddtheta);

vNx = vMx - dtheta*(yN - yM);
vNy = vMy + dtheta*(xN - xM);

aNx = aMx - ddtheta*(yN - yM) - dtheta^2*(xN - xM);
aNy = aMy + ddtheta*(xN - xM) - dtheta^2*(yN - yM);

out = [vNx vNy aNx aNy];
end