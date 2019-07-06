% example 3.6
clear all; clc; close all 

R0n = 4.0;  % (m)
R1n = 14.0; % (m)
an =0.0889; % (m/s^2)

syms t a dr R0 R1

v=a*t;
r=int(v);

lpe=int(v,t, 0,t);
rpe=int(1,dr, R0,R1);
eq=lpe-rpe;
soleq=solve(eq,t);
tf=simplify(soleq); 

slist={R0,R1,a};
nlist={R0n,R1n,an};

tfn=abs(subs(tf,slist,nlist))

tn = 0:0.1:tfn;
rn = subs(R0n+r,{a,t},{an,tn});

polar(tn,rn)

% end of program