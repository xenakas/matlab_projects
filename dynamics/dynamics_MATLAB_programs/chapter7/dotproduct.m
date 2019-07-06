function val = dotproduct(a,b)

% symbolic dot product function
% a.b=a_x*b_x+a_x*b_x+a_x*b_x
% function accepts sym as the input argument
val = sum(a.*b);
