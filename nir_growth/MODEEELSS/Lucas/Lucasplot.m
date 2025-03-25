%Plot of the results of Lucas-Uzawa

for i=1:max(size(t))
    if t(i)<40
        Ende=i;
    end
end

subplot(3,2,1)
plot(t,c,'b-')
hold on
plot(0,c(1),'bx')
title('i')
xlabel('                                              t');
ylabel('c');
axis([0 40 -inf inf])

subplot(3,2,3)
plot(t,h,'b-')
hold on
plot(0,h(1),'bx')
title('iii')
xlabel('                                              t');
ylabel('h');
axis([0 40 -inf inf])

subplot(3,2,4)
plot(t,k,'b-')
hold on
plot(0,k(1),'bx')
title('iv')
xlabel('                                              t');
ylabel('k');
axis([0 40 -inf inf])

subplot(3,2,5)
plot(t,u,'b-')
hold on
plot(0,u(1),'bx')
title('v')
xlabel('                                              t');
ylabel('u');
axis([0 40 -inf inf])

ssk=0;
ssh=0;
for iii=0:max(k)*1.4/1000:max(k)*1.4
    kss=iii;
    hss=(A*beta*kss^(beta-1)*uss^(1-beta)/(rho+gr2*sigma))^(1/(-1+beta-gamma));
    ssk=[ssk kss];
    ssh=[ssh hss];  
end

subplot(3,2,6)
plot(ssk,ssh,'k-')
hold on
plot(k,h,'b-')
hold on
plot(k(1),h(1),'rx')
plot(k(end),h(end),'gx')
title('vi')
xlabel('                                              k');
ylabel('h');
% print -depsc 'Figure1.eps';
pause
close