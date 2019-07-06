% Example 5.3-1
% FM Demodulation

function f5p3d1()

    % px(1:3)=Error covariances, px(4:5)=State estimates
    % Sampling period T=0.1 sec, Measurement every 1 sec.
    px=[10 10 10 10 10];
    i=0; wc=1; r=1; del=0.1; del_mu=1;
    global sig; sig=1;
    rnd=randn(1000,1);
    for t=0:del:10  
        
        % Time Update
        t1=t+del;
        [td px] = ode45(@diffeqns,[t t1],px);
        px=px(end,:);
        
        % Measurement Update
        if (mod(t1,del_mu)==0)
            i=i+1;
            h=[0 sqrt(2)*sig*cos(wc*t1+px(5))];
            d=2*sig^2*px(3)*cos(wc*t1+px(5))+r/del;
            % Kalman Gain
            k=[px(2)*h(2)/d ; px(3)*h(2)/d]; 
            % Error covariances
            px(1)=px(1)-k(1)*h(2)*px(2);
            px(2)=(1-k(2)*h(2))*px(2);
            px(3)=(1-k(2)*h(2))*px(3);
            % Available measurement zk
            z=sqrt(2)*sig*sin(wc*t1+px(5))+rnd(i);
            % State estimate after measurement update
            px(4)=px(4)+k(1)*(z-sqrt(2)*sig*sin(wc*t1+px(5)));
            px(5)=px(5)+k(2)*(z-sqrt(2)*sig*sin(wc*t1+px(5)));
        end
        
        % Plotting
        subplot(3,1,1); plot(t1,px(1),'k.'); hold on; ylabel('p_1'); grid on;
        subplot(3,1,2); plot(t1,px(2),'k.'); hold on; ylabel('p_2'); grid on;
        subplot(3,1,3); plot(t1,px(3),'k.'); hold on; ylabel('p_4'); grid on;
        xlabel('Time (Secs.)');
        
    end  
return

function dpx=diffeqns(t,px)
    global sig;
    a=1;
    % Differential equations for Error covariance
    dpx(1)=-2*a*px(1)+2*a*sig^2;
    dpx(2)=-a*px(2)+px(1);
    dpx(3)=2*px(2);
    % Differential equations for states
    dpx(4)=-a*px(4);
    dpx(5)=px(4);
    dpx=dpx';
return
