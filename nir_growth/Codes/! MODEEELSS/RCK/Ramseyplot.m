%Plot of the results of the Ramsey Cass Koopmans model

%Calculation of "exact" analytical solution
zstart=start(1)^(1-alpha);
koefa=(delta+nPop+gx);
z=1/(theta*koefa)+(zstart-1/(theta*koefa)).*exp(-(1-alpha)*koefa*t);
kan=z.^(1/(1-alpha));
can=kan.^alpha.*(1-1/theta);

%calculation of difference between analytical solution and relaxation
%solution
relerrork=(k-kan)./kan;
relerrorc=(c-can)./can;
disp(['max error c: ',num2str(max(abs(relerrorc)))]);
disp(['max error k: ',num2str(max(abs(relerrork)))]);
disp(['mean error: ',num2str((sqrt(relerrork*relerrork'+relerrorc*relerrorc'))/(N*M))]);


plot(abs(relerrork))
title('Ramsey-Cass-Koopman model:  relative error of capital')
ylabel('error');
xlabel('grid point');
pause
close
plot(abs(relerrorc))
title('Ramsey-Cass-Koopman model:  relative error of consumption')
ylabel('error');
xlabel('grid point');
pause
close
plot(k,c,'bx-')
hold on
plot(kan,can,'g-')
hold on
plot(kss,css,'r+')
title('Ramsey-Cass-Koopman model: phase diagram (k,c)')
xlabel('k');
ylabel('c');
pause
close
plot(t,c,'bx-')
hold on
plot(t,can,'g-')
title('Ramsey-Cass-Koopman model: consumption')
xlabel('t');
ylabel('c');
axis([0 100 -inf inf])
pause
close
plot(t,k,'bx-')
hold on
plot(t,kan,'g-')
title('Ramsey-Cass-Koopman: capital')
axis([0 100 -inf inf])
xlabel('t');
ylabel('k');
pause
close

