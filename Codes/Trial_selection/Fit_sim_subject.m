% FIT SIMULATED SUBJECT

function [AIC winning_model winning_params] = Fit_sim_subject (D, initialVars)

[fitmodels,  nfmodels, maxnumpars, nparshat_vec] = destructure_initial_variables(initialVars);


LL = nan(1,nfmodels);
AIC = nan(1,nfmodels);
parshatM = nan(nfmodels, maxnumpars);


for fitmodelidx = 1:nfmodels % fit all models to each subject
    f = fitmodels(fitmodelidx)
    
    nparshat = nparshat_vec(f);
    
    [AICtemp, LLtemp, parshattemp] = test_models_ga(f,D); 
    
    AIC(1,fitmodelidx) = AICtemp;
    LL(1,fitmodelidx) = LLtemp;
    parshatM(fitmodelidx,1:nparshat) = parshattemp;
    
end

% find winning model:
idx = find(AIC(:) == max(AIC(:)));
location_wm = idx(randi(length(idx)));
winning_model = fitmodels(location_wm);
winning_params = parshatM(location_wm, :);

end
