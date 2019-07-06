% vaRRR.m
% velocity and acceleration RRR dyad 
function out=vaRRR...
(xM,yM,xN,yN,xP,yP,vMx,vMy,vNx,vNy,aMx,aMy,aNx,aNy);

% velocity
vPxs=sym('vPxs','real');
vPys=sym('vPys','real');
eqRRR1v =...
(xM-xP)*(vMx-vPxs)+(yM-yP)*(vMy-vPys);
eqRRR2v =...
(xN-xP)*(vNx-vPxs)+(yN-yP)*(vNy-vPys);

solRRRv=solve(eqRRR1v, eqRRR2v,'vPxs','vPys');
vPx = eval(solRRRv.vPxs);
vPy = eval(solRRRv.vPys);

% acceleration
aPxs=sym('aPxs','real');
aPys=sym('aPys','real');
eqRRR1a =...
(xM-xP)*(aMx-aPxs)+(vMx-vPx)^2+...
(yM-yP)*(aMy-aPys)+(vMy-vPy)^2;
eqRRR2a =...
(xN-xP)*(aNx-aPxs)+(vNx-vPx)^2+...
(yN-yP)*(aNy-aPys)+(vNy-vPy)^2;
solRRRa=solve(eqRRR1a, eqRRR2a,'aPxs','aPys');
aPx = eval(solRRRa.aPxs);
aPy = eval(solRRRa.aPys);

out = [vPx vPy aPx aPy];
end