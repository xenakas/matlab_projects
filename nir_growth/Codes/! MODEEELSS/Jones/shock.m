%shock.m
if t==-1        % At time -1 the parameters should have the values for the old steady state
    alphaF=1;
    alphaJ=1;
else
    alphaF=1.3;     % These are the parameter values for the new steady state
    alphaJ=0.9;
end
