% Example 2.3-1
% Program to plot the apriori and aposteriori error covariance
x=[0;10]; p=[2 0 ; 0 3]; % Initialization
A=[1 1 ; 0 1]; G=[0;1]; Q=1; R=2; H=[1 0];
p22(1,1)=p(2,2);
for i=1:10
    pm=A*p*A'+G*Q*G';
    pm22(i+1,:)=pm(2,2); % apriori error covariance p(2,2) for velocity s
    p=inv(inv(pm)+H'*inv(R)*H);
    p22(i+1,:)=p(2,2); % apostiriori error covariance p(2,2) for velocity s
end
plot([0:10],p22,'k.',[1:10],pm22([2:11]),'k*'); % Plotting
grid on; legend('aposteriori error covariance','apriori error covariance');
xlabel('k'); ylabel('Velocity error variances');