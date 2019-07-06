% vaRRT.m
% velocity and acceleration RRT dyad 
function out=vaRRT(xM,yM,xN,yN,xP,yP,...
 vMx,vMy,vNx,vNy,aMx,aMy,aNx,aNy,dtheta,ddtheta);

theta = atan((yP-yN)/(xP-xN));

% velocity
syms vPxSol vPySol
eqRRT1v =...
(xM-xP)*(vMx-vPxSol)+(yM-yP)*(vMy-vPySol);
eqRRT2v =...
sin(theta)*(vPxSol-vNx)+cos(theta)*dtheta*(xP-xN)...
-cos(theta)*(vPySol-vNy)+sin(theta)*dtheta*(yP-yN)';

solRRTv = solve(eqRRT1v, eqRRT2v,'vPxSol','vPySol');
vPx = eval(solRRTv.vPxSol);
vPy = eval(solRRTv.vPySol);

% acceleration
syms aPxSol aPySol
eqRRT1a =...
(xM-xP)*(aMx-aPxSol)+(vMx-vPx)^2+...
(yM-yP)*(aMy-aPySol)+(vMy-vPy)^2;
eqRRT2a =...
sin(theta)*(aPxSol-aNx)-cos(theta)*(aPySol-aNy)...
+(2*cos(theta)*(vPx-vNx)-sin(theta)*dtheta*(xP-xN)...
+2*sin(theta)*(vPy-vNy)...
+cos(theta)*dtheta*(yP-yN))*dtheta...
+(cos(theta)*(xP-xN)+sin(theta)*(yP-yN))*ddtheta';
solRRTa = solve(eqRRT1a, eqRRT2a,'aPxSol','aPySol');

aPx = eval(solRRTa.aPxSol);
aPy = eval(solRRTa.aPySol);

out = [vPx vPy aPx aPy];
end