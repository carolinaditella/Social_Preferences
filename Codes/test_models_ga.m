function [AIC, LL, parshat] = test_models_ga (fitmodel,D)

nvars = [3,2,2,2,3,3,1,4];


LowerBounds = [...    
    0,0,0,nan;      % FEHR-SCHMIDT
    0,nan,nan,nan;  % LINEAR
    0,0,nan,nan;    % RALWSIAN
    0,0,nan,nan;    % COBB-DOUGLAS
    0,-1,-1,nan;    % CHARNESS-RABIN
    0,0,-2,nan;     % CES
    0,nan,nan,nan;  % LEXICOGRAPHIC
    0,0,0,0];       % COX-CES
   
UpperBounds = [...
    1001,1,1,nan;
    1001,nan,nan,nan;
    1001,1,nan,nan;
    1001,1,nan,nan;
    1001,1,1,nan;
    1001,1,2,nan;
    1001,nan,nan,nan;
    1001,1,1,1];


[x,fval] = ga(@(pars) loglike(pars,D,fitmodel), nvars(fitmodel),[],[],[],[], LowerBounds(fitmodel, 1:nvars(fitmodel)), UpperBounds(fitmodel, 1:nvars(fitmodel)),[],gaoptimset('PopulationSize',100,'Generations',15,'TolFun', 10^(-16)));

parshat = x;
LL = -fval;
AIC = -fval - nvars(fitmodel);

end




        