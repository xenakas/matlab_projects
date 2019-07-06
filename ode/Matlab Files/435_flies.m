%% Swarm of flies from the Lorenz equations
%% Run M-file solvem.m first

hold on

xlim([-20 20])
ylim([-30 30])
zlim([0 50])

for m=1:(d+1)^3;
   hend(m)=plot3(y0(1,m),y0(2,m),y0(3,m),'Erase','background');
   hstart(m)=plot3(y0(1,m),y0(2,m),y0(3,m),'Erase','none');
end

pause=input('Rotate axes now, then press [RETURN]');

for i=floor(length(yi)/3):length(yi)-1,

   for m=1:(d+1)^3
      set(hstart(m),'xdata',yi(i,1,m),'ydata',yi(i,2,m),'zdata',yi(i,3,m));
      drawnow
      set(hend(m),'xdata',yi(i,1,m),'ydata',yi(i,2,m),'zdata',yi(i,3,m));
      drawnow
   end

end
