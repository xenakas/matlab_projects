function GeneralEqmConditions=HeteroAgentStationaryEqm_Case1_subfn(p, V0Kron, n_d, n_a, n_s, n_p, pi_s, d_grid, a_grid, s_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Parameters, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions)

N_d=prod(n_d);
N_a=prod(n_a);
N_s=prod(n_s);
N_p=prod(n_p);

l_p=length(n_p);

%% 
for ii=1:l_p
    Parameters.(GEPriceParamNames{ii})=p(ii);
end

%     ReturnFnParams(IndexesForPricesInReturnFnParams)=p;
[~,Policy]=ValueFnIter_Case1(V0Kron, n_d,n_a,n_s,d_grid,a_grid,s_grid, pi_s, ReturnFn, Parameters, DiscountFactorParamNames,ReturnFnParamNames,vfoptions);

%Step 2: Calculate the Steady-state distn (given this price) and use it to assess market clearance
StationaryDistKron=StationaryDist_Case1(Policy,n_d,n_a,n_s,pi_s,simoptions);
SSvalues_AggVars=SSvalues_AggVars_Case1(StationaryDistKron, Policy, SSvaluesFn, Parameters, SSvalueParamNames, n_d, n_a, n_s, d_grid, a_grid, s_grid,2); % The 2 is for Parallel (use GPU)

% The following line is often a useful double-check if something is going wrong.
%    SSvalues_AggVars

% use of real() is a hack that could disguise errors, but I couldn't
% find why matlab was treating output as complex
GeneralEqmConditionsVec=real(GeneralEqmConditions_Case1(SSvalues_AggVars,p, GeneralEqmEqns, Parameters,GeneralEqmEqnParamNames));

if heteroagentoptions.multiGEcriterion==0 %only used when there is only one price 
    GeneralEqmConditions=GeneralEqmConditionsVec;
elseif heteroagentoptions.multiGEcriterion==1 %the measure of market clearance is to take the sum of squares of clearance in each market 
    GeneralEqmConditions=sqrt(sum(GeneralEqmConditionsVec.^2));                                                                                                         
end

GeneralEqmConditions=gather(GeneralEqmConditions);

if heteroagentoptions.verbose==1
    [p,GeneralEqmConditions]
end

end