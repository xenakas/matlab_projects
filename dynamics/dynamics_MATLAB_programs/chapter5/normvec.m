function val = normvec(v_x,v_y,v_z)
%norm The symbolic norm function
%     of a vector v=v_x*i+v_y*j+v_z*k. This function
%     accepts a sym as the input argument.
val=sqrt(factor(v_x^2+v_y^2+v_z^2));