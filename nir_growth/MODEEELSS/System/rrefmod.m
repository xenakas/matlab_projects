function [A,b] = rrefmod(A,b)
%RREFMOD   Reduced row echelon form.
%   form: [B , x] = RREFMOD(A,b)
%   The modified algorithm transforms the rows of the full rank (m x n) matrix A (m<=n) such that B has the
%   following shape: B = IR with Identity matrix I and a (m x n-m) matrix R. In addition the
%   algorithm performs the same operations with the vector b. Output is the
%   transformed matrix B and corresponding vector x. For transforming the
%   rows the algorithm includes a column pivot search.
%
%   See RREF.

[m,n] = size(A);
q=length(b);
if q~=m
    disp('Warning: Wrong dimension of vector!');
end

% Compute the default tolerance if none was provided.
 tol = max(m,n)*eps*norm(A,'inf'); 

% Loop over the entire matrix.
i = 1;
j = 1;
while (i <= m) & (j <= n)
   % Find value and index of largest element in the remainder of column j.
   [p,k] = max(abs(A(i:m,j))); k = k+i-1;
   if (p <= tol)
      % The column is negligible, zero it out.
      A(i:m,j) = zeros(m-i+1,1);
      j = j + 1;
   else
      % Swap i-th and k-th rows.
      A([i k],j:n) = A([k i],j:n);
      b([i k])=b([k i]);
      % Divide the pivot row by the pivot element.
      b(i)=b(i)/A(i,j);
      A(i,j:n) = A(i,j:n)/A(i,j);
      % Subtract multiples of the pivot row from all the other rows.
      for k = [1:i-1 i+1:m]
         b(k)=b(k)-A(k,j)*b(i);
         A(k,j:n) = A(k,j:n) - A(k,j)*A(i,j:n);
      end
      i = i + 1;
      j = j + 1;
   end
end

% B=inv(A(:,1:m));
% b=B*b;
% if n>m
%     A(:,m+1:n)=B*A(:,m+1:n);
% end
% A(:,1:m)=eye(m);