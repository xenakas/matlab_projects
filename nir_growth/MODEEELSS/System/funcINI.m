function funcINI = funcINI(t,x)
% System file of the relaxation algorithm by Trimborn, Koch, and Steger.
% Copyright by Trimborn, Koch, Steger, 2008

globalpar; 
%extraction of variables from x
varex;
%include shocks
shock;
%Evaluation of residuals
initbound;