% example 6.3(b)
% dynamics of a pendulum
% using the function: R(t,x) 
% defined in an external program R.m

clear all; clc; close all

tfinal=5;
time=[0 tfinal];
% x(1)(0)=0  x(2)(0)=0 
x0=[0 0]; 

% solve a second order ODE 
% using MATLAB ODE45 function

[t,x]=ode45(@R,time,x0);

subplot(3,1,1),...
plot(t,x(:,1),'r'),...
xlabel('t'),ylabel('\theta'),grid,...
subplot(3,1,2),...
plot(t,x(:,2)),...
xlabel('t'),ylabel('\omega'),grid,...
subplot(3,1,3),...
plot(x(:,1),x(:,2)),...
xlabel('\theta'), ylabel('\omega'),grid

% end of program