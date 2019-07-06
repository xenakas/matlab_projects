% Example 2.4-6
% Multi-rate Kalman Filter for radar tracking
% Time update every 50msec (T) and measurement update every 1 sec (Td)

% Define/Initialize variables
clear all; close all
i=1; T=0.05; Td=1;
sigpsq=150*T; sigvsq=150*T; sigmsq=30/Td; % Qs=Q*T, Rs=Q/Td 
Q=[sigpsq 0 ; 0 sigvsq];
H=[1 0];
p=[30 0 ; 0 30]; p1(1)=p(1,1); p4(1)=p(2,2);

% Loop for calculating the error variances
for t=0:T:10-T
    A=[1 T ; 0 1]; i=i+1;
    p=A*p*A'+Q;
    k=p*H'*inv(H*p*H'+sigmsq);
    % Take Measurement at T=1,2...10 secs
    if (mod(t,1)==0 && t~=0)        
        p=(eye(2)-k*H)*p;
    end
    p1(i)=p(1,1); p4(i)=p(2,2);
end

% Plotting
figure(1); plot([0:0.05:10],p1,'k'); hold on; plot([0:0.05:10],p4,'k-','linewidth',2);
legend('p_1','p_4'); grid on;
xlabel('Time (Secs.)'); ylabel('Error Variances');