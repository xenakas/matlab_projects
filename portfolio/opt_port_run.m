clear; clc;

% P --- матрица цен
% R --- матрица доходностей
% E_port --- ожидаемая доходность портфеля

P = xlsread('data.xls','prices');
R = price2ret(P);
E_port = 0.0031;


E = mean(R)'
V = cov(R)
[w_opt_port, D_opt_port] = get_opt_port(E, V, E_port)