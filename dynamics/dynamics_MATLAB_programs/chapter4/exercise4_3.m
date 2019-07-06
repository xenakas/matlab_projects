% exercise 4.3

clear all; clc; close 

syms x y z G m M xA yA zA xB yB zB

r_v = [x y z];
fprintf('r = [%s %s %s]\n',...
    char(r_v(1)),char(r_v(2)),char(r_v(3)))

r = sqrt(x^2+y^2+z^2);

F = -G*m*M*r_v/r^3;
fprintf('Fx = %s  \n',char(F(1)))
fprintf('Fy = %s  \n',char(F(2)))
fprintf('Fz = %s  \n',char(F(3)))

% potential energy
V = -F*r_v.'; % .? array transpose
V = simplify(V);
fprintf('V = %s \n',char(V))

slist={x,y,z};
nlistA={xA,yA,zA};
nlistB={xB,yB,zB};

VA=subs(V,slist,nlistA);
VB=subs(V,slist,nlistB);

VAB = simplify(VB-VA);
fprintf('VAB = VB-VA = %s \n',char(VAB ))

% work
U = simplify(F*r_v.');
UA=subs(U,slist,nlistA);
UB=subs(U,slist,nlistB);
UAB = simplify(UB-UA);
fprintf('UAB = UB-UA = %s \n',char(UAB))

if VAB==-UAB; fprintf('UAB = -VAB \n');
else fprintf('UAB not equal to -VAB \n'); 
end;

% end of program