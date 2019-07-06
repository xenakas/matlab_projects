%% find k for Exercise 9.5(ii)
%% first guess

k=0;

w=pi/12;

%% change N for more iterations and so
%% a more accurate value for k

N=20;

for i=1:N;

    dn=12-(5*k*(k*cos(6*w)+w*sin(6*w))/(k^2+w^2));
    nn=17-(5*k*(k*cos(5*w)+w*sin(5*w))/(k^2+w^2));

    x=dn/nn;

    nk=-log(x);

    k=nk

end
