% pRRT.m
% position RRT dyad 

function out = pRRT(xM, yM, xN, yN, MP, Theta);
xP=sym('xP','real');
yP=sym('yP','real');
if Theta==pi/2 || Theta==3*pi/2 
   xP1 = xN; xP2 = xN;
   eqRRT = (xM-xN)^2+(yM-yP)^2-MP^2;
   solRRT = solve(eqRRT);
   yP1 = eval(solRRT(1));
   yP2 = eval(solRRT(2));
elseif Theta==0 || Theta==pi
    yP1 = yN; yP2 = yN;
    eqRRT = (xM-xP)^2+(yM-yN)^2-MP^2;
    solRRT = solve(eqRRT);
    xP1 = eval(solRRT(1));
    xP2 = eval(solRRT(2));
else
    eqRRT1 = (xM-xP)^2 + (yM-yP)^2 - MP^2';
    eqRRT2 = tan(Theta) - (yP-yN)/(xP-xN)';
    solRRT = solve(eqRRT1, eqRRT2);
    xPpos = eval(solRRT.xP);
    yPpos = eval(solRRT.yP);
    xP1 = xPpos(1); xP2 = xPpos(2); 
    yP1 = yPpos(1); yP2 = yPpos(2);
  end  
out = [xP1 yP1 xP2 yP2];
end
