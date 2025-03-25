function [w_opt_port, D_opt_port] = get_opt_port(E, V, E_port)
n = length(E);

% Начальная точка в процессе оптимизации.
w0 = (1/n) * ones(n,1);
%__________________________________________________________________________
% Параметры для функции fmincon.
fun = @(w) get_port_var(w, V);
x0 = w0;
A = [];
b = [];
Aeq = [ones(1,n); 
       E'];
beq = [1; E_port];
lb = zeros(n,1);
ub = ones(n,1);
nonlcon = [];
options = optimset('Algorithm', 'interior-point', 'Display', 'off');

%__________________________________________________________________________
% Минимизация дисперсии портфеля при помощи функции
% fmincon.
[w_opt_port, D_opt_port] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%__________________________________________________________________________



function [D] = get_port_var(w, V)
D = w' * V * w;