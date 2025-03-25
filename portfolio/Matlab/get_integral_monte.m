function [I_monte] = get_integral_monte(a, b)

n = 100000000;
X = (b-a) * rand(n,1) + a;
X = get_g(X);
I_monte = (b-a) * mean(X);

function [g] = get_g(x)
g = exp(sin(x));