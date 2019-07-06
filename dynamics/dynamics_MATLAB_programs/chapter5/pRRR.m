% pRRR.m
% position RRR dyad 

function out = pRRR(xM, yM, xN, yN, MP, PN);
xP=sym('xP','real');
yP=sym('yP','real');

eqRRR1 = (xM-xP)^2+(yM-yP)^2-MP^2;
eqRRR2 = (xN-xP)^2+(yN-yP)^2-PN^2;

solRRR = solve(eqRRR1, eqRRR2);
xPpos = eval(solRRR.xP);
yPpos = eval(solRRR.yP);
xP1 = xPpos(1); xP2 = xPpos(2); 
yP1 = yPpos(1); yP2 = yPpos(2); 

out = [xP1 yP1 xP2 yP2];
end