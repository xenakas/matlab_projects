% Example 2.4-5
% Define variables
close all; clear all;
T=1;
sigpsq=150; sigvsq=150; sigmsq=30;
A=[1 T ; 0 1];
Q=[sigpsq 0 ; 0 sigvsq];
H=[1 0];

% Initial values 
p=[30 0 ; 0 30];
alpha(1)=0.5;
beta(1)=0;

% Calculate error covariances and filter gains
for i=1:10
    pm=A*p*A'+Q; % apriori error covariance
    k=pm*H'*inv(H*pm*H'+sigmsq); % kalman gain
    p=(eye(2)-k*H)*pm; % aposteriori error covariance
    alpha(i+1)=k(1); 
    beta(i+1)=k(2);
end

% Plot the Filter Gains
figure(1); plot([0:10],alpha,'k-','linewidth',2); hold on; plot([0:10],beta,'k--','linewidth',2); grid on;
xlabel('k'); ylabel('Filter Gains');
legend('\alpha_k','\beta_k_/_T','location',[0.72 0.7 0.16 0.15])