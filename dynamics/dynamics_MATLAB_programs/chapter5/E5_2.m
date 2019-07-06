% example 5.2

clear all; clc; close all

syms theta r t omega alpha  v_B21 a_B21 

r_AB=[r 0 0];
v_A=[0 0 0];
omega_v=[0 0 omega];

v_B1=v_A+cross(omega_v,r_AB);

fprintf('v_B1 = ')
pretty(v_B1); fprintf('\n\n')

v_B2B1=[v_B21 0 0];

v_B2=v_B1+v_B2B1;
fprintf('v_B2=')
pretty(v_B2); fprintf('\n\n')

% a_B2
a_A=diff(v_A,t);
alpha_v=[0 0 alpha];

a_B1=a_A+cross(alpha_v,r_AB)+...
    cross(omega_v,cross(omega_v,r_AB));
fprintf('a_B1=')
pretty(a_B1); fprintf('\n\n')

a_B2B1=[a_B21 0 0];

a_B21cor=2*cross(omega_v,v_B2B1);
fprintf('Coriolis acceleration: a_B21cor =')
pretty(a_B21cor); fprintf('\n\n')

a_B2=a_A+a_B2B1+...
    cross(alpha_v,r_AB)+...
    2*cross(omega_v,v_B2B1)+...
    cross(omega_v,cross(omega_v,r_AB));
fprintf('a_B2=')
pretty(a_B2); fprintf('\n\n')

% x0 y0 z0 reference frame

RF0=...
[  cos(theta) sin(theta) 0 ;
  -sin(theta) cos(theta) 0 ;
   0 0 1];

r_AB0=r_AB*RF0;
fprintf('r_AB0 = ')
pretty(r_AB0); fprintf('\n\n')

v_B10=v_A+cross(omega_v,r_AB0);

fprintf('v_B10=')
pretty(v_B10); fprintf('\n\n')

m_vB1=normvec(v_B1(1),v_B1(2),v_B1(3));
m_vB10=simplify...
    (normvec(v_B10(1),v_B10(2),v_B10(3)));

fprintf('|v_B1| = ')
pretty(m_vB1); fprintf('\n')
fprintf('|v_B10| = ')
pretty(m_vB10); fprintf('\n\n')

v_B2B10=v_B2B1*RF0;
v_B20=v_B10+v_B2B10;

fprintf('v_B20x=')
pretty(simplify(v_B20(1)))
fprintf('v_B20y=')
pretty(simplify(v_B20(2)))
fprintf('v_B20z=')
pretty(simplify(v_B20(3)))
fprintf('\n\n')

a_B21cor0=2*cross(omega_v,v_B2B10);
fprintf('Coriolis acceleration RF0: a_B21cor0 =')
pretty(a_B21cor0); fprintf('\n\n')

a_B10=a_A+cross(alpha_v,r_AB0)+...
    cross(omega_v,cross(omega_v,r_AB0));

fprintf('a_B10x=')
pretty(simplify(a_B10(1)))
fprintf('a_B10y=')
pretty(simplify(a_B10(2)))
fprintf('a_B10z=')
pretty(simplify(a_B10(3)))
fprintf('\n\n')

a_B2B10=a_B2B1*RF0;
a_B21cor0=2*cross(omega_v,v_B2B10);

a_B20=a_A+a_B2B10+...
    cross(alpha_v,r_AB0)+...
    2*cross(omega_v,v_B2B10)+...
    cross(omega_v,cross(omega_v,r_AB0));

fprintf('a_B20x=')
pretty(simplify(a_B20(1)))
fprintf('a_B20y=')
pretty(simplify(a_B20(2)))
fprintf('a_B20z=')
pretty(simplify(a_B20(3)))
fprintf('\n\n')

% numerical data
slist={theta,r,omega,alpha,v_B21,a_B21};
nlist={pi/4,1,-1,-2,3,2};

v_B1n=subs(v_B1,slist,nlist);
v_B10n=subs(v_B10,slist,nlist);

v_B2n=subs(v_B2,slist,nlist);
v_B20n=subs(v_B20,slist,nlist);

a_B1n=subs(a_B1,slist,nlist);
a_B10n=subs(a_B10,slist,nlist);

a_B21corn=subs(a_B21cor,slist,nlist);
a_B21cor0n=subs(a_B21cor0,slist,nlist);

a_B2n=subs(a_B2,slist,nlist);
a_B20n=subs(a_B20,slist,nlist);

fprintf...
('vB1=[%6.3f %6.3f %g] (m/s)\n',v_B1n)
fprintf...
('RF0: vB10=[%6.3f %6.3f %g] (m/s)\n',v_B10n)
fprintf...
('|vB1|=|vB10|=%6.3f (m/s)\n',norm(v_B10n))
fprintf('\n')

fprintf...
('vB2=[%6.3f %6.3f %g] (m/s)\n',v_B2n)
fprintf...
('RF0: vB20=[%6.3f %6.3f %g] (m/s)\n',v_B20n)
fprintf...
('|vB2|=|vB20|=%6.3f (m/s)\n',norm(v_B20n))
fprintf('\n')

fprintf...
('aB1=[%6.3f %6.3f %g] (m/s^2)\n',a_B1n)
fprintf...
('RF0: aB10=[%6.3f %6.3f %g] (m/s^2)\n',a_B10n)
fprintf...
('|aB1|=|aB10|=%6.3f (m/s^2)\n',norm(a_B10n))
fprintf('\n')

fprintf...
('aB21cor=[%6.3f %6.3f %g] (m/s^2)\n',a_B21corn)
fprintf...
('RF0: aB21cor0=[%6.3f %6.3f %g] (m/s^2)\n',...
a_B21cor0n)
fprintf...
('|aB21cor|=|aB21cor0|=%6.3f (m/s^2)\n',...
norm(a_B21cor0n))
fprintf('\n')

fprintf...
('aB2=[%6.3f %6.3f %g] (m/s^2)\n',a_B2n)
fprintf...
('RF0: aB20=[%6.3f %6.3f %g] (m/s^2)\n',a_B20n)
fprintf...
('|aB2|=|aB20|=%6.3f (m/s^2)\n',norm(a_B20n))
fprintf('\n')

% end of program

% 
% v_B1 = 
%   +-             -+
%   | 0, omega r, 0 |
%   +-             -+
% 
% 
% v_B2=
%   +-                 -+
%   | v_B21, omega r, 0 |
%   +-                 -+
% 
% 
% a_B1=
%   +-       2              -+
%   | - omega  r, alpha r, 0 |
%   +-                      -+
% 
% 
% Coriolis acceleration: a_B21cor =
%   +-                   -+
%   | 0, 2 omega v_B21, 0 |
%   +-                   -+
% 
% 
% a_B2=
%   +-             2                              -+
%   | a_B21 - omega  r, alpha r + 2 omega v_B21, 0 |
%   +-                                            -+
% 
% 
% r_AB0 = 
%   +-                             -+
%   | r cos(theta), r sin(theta), 0 |
%   +-                             -+
% 
% 
% v_B10=
%   +-                                          -+
%   | -omega r sin(theta), omega r cos(theta), 0 |
%   +-                                          -+
% 
% 
% |v_B1| = 
%         2  2 1/2
%   (omega  r )
% 
% |v_B10| = 
%         2  2 1/2
%   (omega  r )
% 
% 
% v_B20x=
%   v_B21 cos(theta) - omega r sin(theta)
% v_B20y=
%   v_B21 sin(theta) + omega r cos(theta)
% v_B20z=
%   0
% 
% 
% Coriolis acceleration RF0: a_B21cor0 =
%   +-                                                        -+
%   | (-2) omega v_B21 sin(theta), 2 omega v_B21 cos(theta), 0 |
%   +-                                                        -+
% 
% 
% a_B10x=
%                        2
%   - r (cos(theta) omega  + alpha sin(theta))
% a_B10y=
%             2
%   - r (omega  sin(theta) - alpha cos(theta))
% a_B10z=
%   0
% 
% 
% a_B20x=
%                     2
% - r cos(theta) omega  - 2 v_B21 sin(theta) omega +
% 
%    a_B21 cos(theta) - alpha r sin(theta)
% a_B20y=
%                     2
% - r sin(theta) omega  + 2 v_B21 cos(theta) omega +
% 
%    a_B21 sin(theta) + alpha r cos(theta)
% a_B20z=
%   0
% 
% 
% vB1=[ 0.000 -1.000 0] (m/s)
% RF0: vB10=[ 0.707 -0.707 0] (m/s)
% |vB1|=|vB10|= 1.000 (m/s)
% 
% vB2=[ 3.000 -1.000 0] (m/s)
% RF0: vB20=[ 2.828  1.414 0] (m/s)
% |vB2|=|vB20|= 3.162 (m/s)
% 
% aB1=[-1.000 -2.000 0] (m/s^2)
% RF0: aB10=[ 0.707 -2.121 0] (m/s^2)
% |aB1|=|aB10|= 2.236 (m/s^2)
% 
% aB21cor=[ 0.000 -6.000 0] (m/s^2)
% RF0: aB21cor0=[ 4.243 -4.243 0] (m/s^2)
% |aB21cor|=|aB21cor0|= 6.000 (m/s^2)
% 
% aB2=[ 1.000 -8.000 0] (m/s^2)
% RF0: aB20=[ 6.364 -4.950 0] (m/s^2)
% |aB2|=|aB20|= 8.062 (m/s^2)
