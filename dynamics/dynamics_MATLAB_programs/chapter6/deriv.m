function fout = deriv(f, g)
% deriv differentiates f with respect to g=g(t)
% the variable g=g(t) is a function of time
syms t xx dxx
lg = {diff(g, t), g}; 
lx = {dxx, xx};
f1 = subs(f, lg, lx);
f2 = diff(f1, xx);
fout = subs(f2, lx, lg);