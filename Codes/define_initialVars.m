%% define_initialVars 

% For the list of models we worked with:
models_names = {'FEHR-SCHMIDT','UTILITARIAN','RAWLSIAN','COBB-DOUGLAS','CHARNES-RABIN','CES','LEXICOGRAPHIC','COX_CES'};
% models' equatios are defined in loglike.m and models' parameter ranges in
% test_models_ga.m. 

% we decided to test the Charness-Rabin, Cobb-Douglas and Rawlsian models

fitmodels = [5 4 3]; % these later changed to 1,2 and 3 
nfmodels = length(fitmodels);
nparshat_vec = [3 2 2 2 3 3 1 4]; 
maxnumpars = max(nparshat_vec);

initialVars = structure_initial_variables (fitmodels, nfmodels, maxnumpars, nparshat_vec); 
