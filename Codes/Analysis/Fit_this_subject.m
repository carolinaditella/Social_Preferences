function [AIC LL parshatM winning_model winning_params D] = Fit_this_subject (Data, initialVars, analysisVars, partner)


[fitmodels, nfmodels, maxnumpars, nparshat_vec,] = destructure_initial_variables(initialVars);

%[partners,n_partners,s_n_list,n_subjects] = destructure_analysis_variables(analysisVars);


if strcmp(partner,'friend') == 1
    D(:,1:5) = Data.choice_friend(:,1:5);
    D(:,10:12) = Data.choice_friend(:,6:8);
elseif strcmp(partner, 'stranger') == 1
    D(:,1:5) = Data.choice_stranger(:,1:5);
    D(:,10:12) = Data.choice_stranger(:,6:8);
end

D(:,6) = (D(:,1)-D(:,2))>=0; % these are the charness-rabin indicator variables. 
D(:,7) = (D(:,2)-D(:,1))>0;
D(:,8) = (D(:,3)-D(:,4))>=0;
D(:,9) = (D(:,4)-D(:,3))>0;

% exclude nan trials:
exc_trials = find(isnan(D(:,5)));
exc_trials_ii = find(D(:,3)>101 | D(:,4)>101);
D(exc_trials,:) = [];
D(exc_trials_ii,:) = [];


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
idx = find(AIC(:) == min(AIC(:)));
location_wm = idx(randi(length(idx)));
winning_model = fitmodels(location_wm);
winning_params = parshatM(location_wm, :);

end
