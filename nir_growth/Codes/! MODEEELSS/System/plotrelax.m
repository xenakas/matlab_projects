function plotrelax(t, x, n1, time)
% plots the result of relax.m up to 16 variables
% 
% COMMAND:
% plotrelax(t, x, n1, time)
% 
% INPUT ARGUMENTS:
% t: time vector
% x: matrix of variables. Each row i desribes the value of variable i at time vector t
% n1: (optinal, default value is 2) number of state variables
% time: (optional, default value is 100) The variables are plotted from t0
% to time.
%     
% Copyright by Trimborn, Koch, Steger, 2008

if nargin < 4, time=100; end 
if nargin < 3, n1=2; end 
[n,len]=size(x);
p=1;
q=1;
if n==2
    p=2;
elseif n==3
    p=3;
elseif n==4
    p=2;
    q=2;
elseif n==5 | n==6
    p=3;
    q=2;
elseif n>6 & n<10
    p=3;
    q=3;
elseif n>9 & n<13
    p=4;
    q=3;
elseif n>12
    p=4;
    q=4;
end

for i=1:min(16,n)
    subplot(p,q,i);  
    plot(t,x(i,:))
    title(num2str(i));
    xlabel('time');
    ylabel(['variable ',num2str(i)]);
    hold on   
    axis([0 time -inf inf])
end
pause
close

if n1==1 | n1==2 | n==2
    plot(x(1,:),x(2,:))
    hold on
    plot(x(1,end),x(2,end),'rx')
    title('phase diagram');
    xlabel('variable 1');
    ylabel('variable 2');
    pause
    close
elseif n1>2 | n==3
    plot3(x(1,:),x(2,:),x(3,:))
    hold on
    plot3(x(1,end),x(2,end),x(3,end),'rx')
    title('phase diagram');
    xlabel('variable 1');
    ylabel('variable 2');
    zlabel('variable 3');
    pause
    close
end