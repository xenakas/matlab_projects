%% find time of death for Exercise 9.5(iv)

%% value of k from findk.m

k=0.3640;

%% initial guess

t0=2;

w=pi/12;

%% change N for more accurate value

N=20;

for i=1:N;

   dn=17-(5*k*(k*cos(5*w)+w*sin(5*w)))/(k^2+w^2);
   nn=34-(5*k*(k*cos(w*(t0-2))+w*sin(w*(t0-2))))/(k^2+w^2);

   x=dn/nn;

   nt=7+(log(x)/k);

   t0=nt

end
