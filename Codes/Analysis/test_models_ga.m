function [AIC, LL, parshat] = test_models_ga (fitmodel,D)

nvars = [3,2,2,2,3,3,1,4];


LowerBounds = [...    
    0,0,0,nan;
    0,-1,nan,nan;
    0,0,nan,nan;
    0,0,nan,nan;
    0,-1,-1,nan;
    0,0,-2,nan;
    0,nan,nan,nan;
    0,0,0,0];
   
UpperBounds = [...
    1001,1,1,nan;
    1001,1,nan,nan;
    1001,1,nan,nan;
    1001,1,nan,nan;
    1001,1,1,nan;
    1001,1,2,nan;
    1001,nan,nan,nan;
    1001,1,1,1];


[x,fval] = ga(@(pars) loglike(pars,D,fitmodel), nvars(fitmodel),[],[],[],[], LowerBounds(fitmodel, 1:nvars(fitmodel)), UpperBounds(fitmodel, 1:nvars(fitmodel)),[],gaoptimset('PopulationSize',100,'Generations',15,'TolFun', 10^(-16)));

parshat = x;
LL = -fval;
AIC = 2*nvars(fitmodel)-2*LL;

end


%[x,fval] = ga(@(pars) loglike(pars,D,fitmodel), nvars(fitmodel),[],[],[],[], LowerBounds(fitmodel, 1:nvars(fitmodel)), UpperBounds(fitmodel, 1:nvars(fitmodel)));
% ,[],gaoptimset('PopulationSize',50,'GenerationsSize',50)
% ,'Display', 'iter')

        