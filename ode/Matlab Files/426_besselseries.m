%% Bessel function of order nu from its series expansion
%% Sum up to power x^(2N)

nu=input('nu = ');
N=input('N= ');

a(1)=-1/(4*(1+nu));

%% a(i) is coefficient of x^(2i)

for i=2:N;
   a(i)=-a(i-1)/(4*i*(i+nu));
end

x=linspace(0,10,1000);

y=1+0.*x;

hold on

for k=1:N;
   y=y+a(k)*x.^(2*k);
   ay=y.*x.^nu;
   if (k/2)~=floor(k/2) & k>1;
      plot(x,ay);
   end
end

%% Change the normalisation

z=besselj(nu,x)*2^nu*gamma(1+nu);

plot(x,z,'linewidth',2)

ylim([1.5*min(z),1.5*max(z)])
