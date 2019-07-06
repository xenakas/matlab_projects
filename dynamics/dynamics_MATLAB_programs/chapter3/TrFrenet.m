function[vTangent,vNormal,vBinormal]=TrFrenet(x,y,z,t) 

[v_r,a_r,magn_v,magn_a] = part_der(x,y,z,t);
% Tangent
vTangent = v_r./magn_v;
dvTangent=diff(vTangent,t);
magn_dvTangent = ...
sqrt(dvTangent(1)^2+dvTangent(2)^2+dvTangent(3)^2);
% Normal
vNormal = dvTangent/magn_dvTangent;
% Binormal
vBinormal = cross(vTangent,vNormal);
