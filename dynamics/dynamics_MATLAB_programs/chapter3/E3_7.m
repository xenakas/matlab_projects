% example 3.7
clear all; clc; close all

syms a b omega t 

% parametric equations
theta = omega*t;
x = a*cos(theta);
y = a*sin(theta);
z = b*t;

% velocity components on x,y,z
v_x=diff(x,t);
v_y=diff(y,t);
v_z=diff(z,t);
vP=[v_x v_y v_z];
% magnitude of the velocity
magn_v=sqrt(vP(1)^2+vP(2)^2+vP(3)^2);

fprintf('velocity components v_x, v_y and v_z')
fprintf('\n')
pretty(vP); fprintf('\n\n')
fprintf('magnitude of the velocity is')
fprintf('\n')
pretty(simplify(magn_v)); fprintf('\n\n')

tau = vP/magn_v;
dtau = diff(tau,t);
magn_dtau = sqrt(dtau(1)^2+dtau(2)^2+dtau(3)^2);
nu = dtau/magn_dtau;
beta = cross(tau,nu);

% using the TrFrenet function
% [Tangent,Normal,Binormal]=TrFrenet(x,y,z,t);
% tau = Tangent;
% nu = Normal;
% beta = Binormal;

tau = simplify(tau);
nu = simplify(nu);
beta = simplify(beta);

fprintf('tau_x = ')
pretty(tau(1)); 
fprintf('\n');
fprintf('tau_y = ')
pretty(tau(2)); 
fprintf('\n');
fprintf('tau_z = ')
pretty(tau(3)); 
fprintf('\n\n');

fprintf('nu_x = ')
pretty(nu(1)); 
fprintf('\n');
fprintf('nu_y = ')
pretty(nu(2)); 
fprintf('\n');
fprintf('nu_z = ')
pretty(nu(3)); 
fprintf('\n\n');

fprintf('beta_x = ')
pretty(beta(1)); 
fprintf('\n');
fprintf('beta_y = ')
pretty(beta(2)); 
fprintf('\n');
fprintf('beta_z = ')
pretty(beta(3)); 
fprintf('\n\n');

fprintf('velocity components on tau, nu, beta')
fprintf('\n');
v_tau = simplify(vP*tau.');
v_nu = simplify(vP*nu.');
v_beta = simplify(vP*beta.');

fprintf('v_tau = ')
pretty(v_tau); 
fprintf('\n')
fprintf('v_nu = ')
pretty(v_nu); 
fprintf('\n')
fprintf('v_beta = ')
pretty(v_beta); 
fprintf('\n\n')

% acceleration components
a_x=diff(x,t,2);
a_y=diff(y,t,2);
a_z=diff(z,t,2);
aP=[a_x a_y a_z];
% magnitude of the acceleration
magn_a=sqrt(aP(1)^2+aP(2)^2+aP(3)^2);

fprintf('acceleration components a_x, a_y and a_z')
fprintf('\n')
pretty(aP); fprintf('\n\n')
fprintf('magnitude of the acceleration is')
fprintf('\n')
pretty(simplify(magn_a)); fprintf('\n\n')

fprintf('acceleration components on tau, nu, beta')
fprintf('\n')
a_tau = simplify(aP*tau.');
a_nu = simplify(aP*nu.');
a_beta = simplify(aP*beta.');

fprintf('a_tau = ')
pretty(a_tau); 
fprintf('\n')
fprintf('a_nu = ')
pretty(a_nu); 
fprintf('\n');
fprintf('a_beta = ')
pretty(a_beta); 
fprintf('\n\n')

omegan=1.0;
an= 1;
bn = omegan/(2*pi);

axis manual
axis equal
axis([-2 2 -2 2 -2 2])

grid on

az = 37.5;
el = 30;
view(az, el);
hold on

scale_factor=1;

start_value=0;
end_value=pi;
step=pi/50;

i=0;
for tn = start_value : step : end_value
   i=i+1;
   slist1={a,b,omega,t};
   nlist1={an,bn,omegan,tn};
   
   xn = subs(x,slist1,nlist1);
   yn = subs(y,slist1,nlist1);
   zn = subs(z,slist1,nlist1);
   
   taun = subs(tau,slist1,nlist1)/scale_factor;
   nun = subs(nu,slist1,nlist1)/scale_factor;
   betan = subs(beta,slist1,nlist1)/scale_factor;
  
   ht = plot3(xn,yn,zn,'k.','Color','r');
   hm = plot3(xn,yn,zn,'k.','Color','b');
   
   title('trajectory of the particle')
   % Frenet frame represented along the trajectory
   ptau=quiver3(xn,yn,zn,taun(1),taun(2),taun(3),...
       'Color','b','LineWidth',1);
   pnu=quiver3(xn,yn,zn,nun(1),nun(2),nun(3),...
       'Color','r','LineWidth',1);
pbeta=quiver3(xn,yn,zn,betan(1),betan(2),betan(3),...
       'Color','k','LineWidth',1);
   
   pause(0.05)

   delete(ht); 
   delete(ptau); 
   delete(pnu);
   delete(pbeta);
       
   x_nn(i) = xn;
   y_nn(i) = yn;
   z_nn(i) = zn;
   for j=1:3
   taunn(i,j)=taun(j);
   nunn(i,j)=nun(j);
   betann(i,j)=betan(j);
   end
end

i=0;
% plot Frenet frame   
for tn = start_value : step : end_value
i=i+1;
    if (i>5) && (mod(i,5)==0)     
    quiver3(x_nn(i),y_nn(i),z_nn(i),...
        taunn(i,1),...
        taunn(i,2),...
        taunn(i,3),'color','b');
    quiver3(x_nn(i),y_nn(i),z_nn(i),...
        nunn(i,1),...
        nunn(i,2),...
        nunn(i,3),'color','r');
    quiver3(x_nn(i),y_nn(i),z_nn(i),...
        betann(i,1),...
        betann(i,2),...
        betann(i,3),'color','k');
    end
    title('Frenet frames')
end

% end of program