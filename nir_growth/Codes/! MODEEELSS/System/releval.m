function releval=releval(tt,t,x);
%Version 3.1
%evaluates variable x at knots tt by linear interpolation, while
%the values of x are only known at knots t.
%Input t and x as maching row vectors or x as maching matrix
%
%example: releval([1 5],t,k) Evaluates capital stock at time 1 and 5
%
%Copyright by Trimborn, Koch, Steger, 2008
[rows, cols]=size(x);
releval=zeros(rows,length(tt)); 
if length(t)~=cols
    disp(['t and x do not match!']);
end  
index=1;
for i=1:length(tt)     
    for ii=1:length(t)
        if tt(i)>t(ii)
            index=ii;
        end
    end
    a=(x(:,index)-x(:,index+1))/(t(index)-t(index+1));
    b=x(:,index)-a*t(index);
    if tt(i)<inf
        releval(:,i)=a*tt(i)+b;
    else
        releval(:,i)=x(:,end);
    end
end
