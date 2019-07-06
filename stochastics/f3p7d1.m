% Example 3.7-1
% Multi-rate Continuous Discrete KF implementation for alpha-beta tracker
% A=2x2, G=2x1, X0=2x1, P0=2x2, Q=1x1, ttu=ix1, H=1x2, z=mx1, R=1x1,
% tmu=mx1
% Example: f3p7d1([0 1;0 0],[0;1],[0;0],[10 10 10 10],1,[0:0.1:10],[1 0],[5 ; 10 ; 15 ; 20 ; 25 ; 30 ; 35 ; 40 ; 45 ; 50]+randn(10,1),1,[1:10])
function f3p7d1(A,G,X0,P0,Q,ttu,H,Z,R,tmu)
     global P1;
    X=X0'; P=P0; j=1;
    for i=1:length(ttu)-1
        % Time update for states
        X(i+1,:)=(A*X(i,:)')'; 
        % Time update for error covariance
        P1=vectortomatrix(P(i,:)); P1=A*P1+P1*A'+G*Q*G'; P1=matrixtovector(P1);
        [t1 P1] = ode45(@tu_ecov,[ttu(i) ttu(i+1)],P(i,:));
        P(i+1,:)=P1(end,:);
        % Measurement Update
        if (ttu(i+1)>=tmu(j)) % if true then do measurement update
            [X(i+1,:),P(i+1,:)]=measurement_update(X(i+1,:)',P(i+1,:),Z(j,1),H,R);
            j=j+1;
        end
    end
    % Plotting
    plot(ttu,P(:,1),'k-',ttu,P(:,4),'k--');
    legend('p_1','p_4'); grid on;
return

% Function used by ode45() for finding P- 
function Pd=tu_ecov(t,P)
    global P1;
    Pd=P1';
return

% Measurement Update
function [X,P]=measurement_update(X,P,Z,H,R)
    P=vectortomatrix(P);
    K=P*H'*inv(H*P*H'+R);
    P=(eye(length(P))-K*H)*P;
    X=X+K*(Z-H*X);
    P=matrixtovector(P);
    X=X';
return

function P1=matrixtovector(P)
    P1=[];
    for i=1:length(P)
        P1=[P1 P(i,:)];
    end
return

function P1=vectortomatrix(P)
    lng=sqrt(length(P)); P1=[];
    for k=1:lng:length(P)
        P1=[P1 ; P(k:k+lng-1)];
    end
return
