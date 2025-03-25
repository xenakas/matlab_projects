%Plot of the results of Jones (1995)



subplot(3,2,1)
plot(t,c,'b-')
hold on
plot(0,c(1),'bx')
plot(0,start(3),'rx')
title('i')
xlabel('                                              t');
ylabel('c');
axis([0 150 -inf inf])

subplot(3,2,2)
plot(t,theta,'b-')
hold on
plot(0,theta(1),'bx')
plot(0,start(5),'rx')
title('ii')
xlabel('                                              t');
ylabel('\phi');
axis([0 150 -inf inf])

subplot(3,2,3)
plot(t,v,'b-')
hold on
plot(0,v(1),'bx')
plot(0,start(4),'rx')
title('iii')
xlabel('                                              t');
ylabel('v_a');
axis([0 150 -inf inf])

subplot(3,2,4)
plot(t,k,'b-')
hold on
plot(0,k(1),'bx')
plot(0,start(1),'rx')
title('iv')
xlabel('                                              t');
ylabel('k');
axis([0 150 -inf inf])

subplot(3,2,5)
plot(t,a,'b-')
hold on
plot(0,a(1),'bx')
plot(0,start(2),'rx')
title('v')
xlabel('                                              t');
ylabel('a');
axis([0 150 -inf inf])

subplot(3,2,6)
plot(k,a,'b-')
hold on
plot(k(1),a(1),'rx')
plot(k(end),a(end),'gx')
title('vi')
xlabel('                                              k');
ylabel('a');
% print -depsc 'Figure1.eps';
pause
close