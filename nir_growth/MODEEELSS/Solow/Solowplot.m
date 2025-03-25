%Plot of the results of Solow model


plot(t,k,'b-')
hold on
plot(0,k(1),'bx')
xlabel('t');
ylabel('k');
axis([0 100 -inf inf])
pause
close
plot(t,c,'b-')
hold on
plot(0,c(1),'bx')
xlabel('t');
ylabel('c');
axis([0 100 -inf inf])
pause
close