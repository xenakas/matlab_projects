% Требуется вычислить интеграл $\int_{-1}^{1}\sqrt{1 - x^2}\,dx = \pi / 2$.

I_theo = pi/2;

n = 1000000;
X = 2 * rand(n,1) - 1;
X = sqrt(1 - X.^2);
I_monte = 2 * mean(X);
