%% makematrix.m
%% construct matrix with specified
%% eigenvalues & eigenvectors

l1=input('First eigenvalue = ');
v1=input('First eigenvector = '); v1=conj(v1');

if imag(l1)~=0
   v2=imag(v1); v1=real(v1);
   M=[real(l1) imag(l1); -imag(l1) real(l1)];
else
   l2=input('Second eigenvalue = ');
   v2=input('Second eigenvector = '); v2=v2';
   M=[l1 0; 0 l2];
   if l1==l2 & norm(v1-v2)~=0;
      M=[l1 1; 0 l1];
   end

end

P=[v1, v2];

A=P*M*inv(P)

[V D]=eig(A)
