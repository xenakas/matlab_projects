% Example 4.1-1
% Program for Kalman Filter performance study
clear all; close all;
k=[0:150];
p0=1; r=10; q=0.01; s=p0;
p=p0./(1+k*(p0/r)); % Predicted error Covariance p
K=p./r; % Kalman Gain K
for i=2:length(k)
    s(i)=s(i-1)+q;
    s(i)=(1-K(i)).^2*s(i)+K(i).^2*r; % Actual Error covariance s
end
plot(k,p,k,s); grid on; legend('p_k','s_k'); xlabel('k'); % Plot p

