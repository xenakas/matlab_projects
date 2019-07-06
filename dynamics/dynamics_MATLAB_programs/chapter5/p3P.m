% p3P.m
% position 3 points on a line 

function out = p3P(xM, yM, xN, yN, MP);
xP=sym('xP','real');
yP=sym('yP','real');
 eqP1 = (xM-xP)^2 + (yM-yP)^2 - MP^2;
 eqP2 = (yP-yM)/(xP-xM)-(yM-yN)/(xM-xN);
 solP = solve(eqP1, eqP2);
 xPpos = eval(solP.xP);
 yPpos = eval(solP.yP);
 xP1 = xPpos(1); xP2 = xPpos(2); 
 yP1 = yPpos(1); yP2 = yPpos(2); 
out = [xP1 yP1 xP2 yP2];
end