function [value,isterminal,direction] = eventE6_8(t,x)

value = x(3);

isterminal = 1;

direction = 0;

% If isterminal vector is set to 1, integration will halt 
% when a zero-crossing is detected.

% The elements of the direction vector are -1, 1, or 0, 
% specifying that the corresponding event must be decreasing, 
% increasing, or that any crossing is to be detected