function [AIC, LL, parshat] = test_collusion_ga (D)

nvars = 1;

LowerBound = 0;
   
UpperBound = 1001;


[x,fval] = ga(@(pars) loglike_collusion(pars,D), nvars,[],[],[],[], LowerBound, UpperBound,[],gaoptimset('PopulationSize',100,'Generations',15,'TolFun', 10^(-16)));

parshat = x;
LL = -fval;
AIC = 2*nvars-2*LL;

end