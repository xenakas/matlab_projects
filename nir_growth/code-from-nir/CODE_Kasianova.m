
randn('state',100)           
T = 10; N = 1000; dt = T/N;
dW = zeros(1,N);             
W = zeros(1,N);              

dW(1) = sqrt(dt)*randn;      
W(1) = dW(1);                
for j = 2:N
   dW(j) = sqrt(dt)*randn;   
   W(j) = W(j-1) + dW(j); 
end

figure(1)
plot([0:dt:T],[0,W],'r-')    
xlabel('t','FontSize',16) 
ylabel('W(t)','FontSize',16,'Rotation',0)


dt = zeros(1,N);             
t = zeros(1,N);              
dt(1) = T/N;
t(1) = dt(1);                
for j = 2:N
    dt(j) = T/N
    t(j) = t(j-1) + dt(j); 
end

figure(2)
plot([0:dt:T],[0,t],'r-')    
xlabel('t','FontSize',16) 
ylabel('t','FontSize',16,'Rotation',0)




m = 0.8;
s = 0.6;
dx = zeros(1,N);             
x = zeros(1,N);              
x(1) = 10;                
for j = 2:N
    x(j) = x(1) * exp(m*t(j)) ; 
   % dx(j) = m * x(j) * dt;
end


figure(3)
plot([0:dt:T],[0,x],'r-')    
xlabel('t','FontSize',16) 
ylabel('X(t)','FontSize',16,'Rotation',0)


dy = zeros(1,N);             
y = zeros(1,N);              
y(1) = 10;                
for j = 2:N
    y(j) = x(1) * exp((m - s*s/2) *t(j) + s*W(j)) ; 
 %dy(j) = m * y(j) * dt + s * x(j) * dW(j);
   
end


figure(4)
plot([0:dt:T],[0,y],'r-')    % plot W against t
xlabel('t','FontSize',16) 
ylabel('X(t)','FontSize',16,'Rotation',0)

