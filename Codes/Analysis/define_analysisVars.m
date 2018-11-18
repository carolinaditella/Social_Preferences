%% define_analysisVars 

% For the list of models we worked with:
models_names = {'FEHR-SCHMIDT','UTILITARIAN','RAWLSIAN','COBB-DOUGLAS','CHARNES-RABIN','CES','LEXICOGRAPHIC','COX_CES'};
% models' equatios are defined in loglike.m and models' parameter ranges in
% test_models_ga.m. 

% we decided to test the Charness-Rabin, Cobb-Douglas and Rawlsian models



partners = {'friend','stranger'};
n_partners = length(partners);
s_n_list = 301:330; % this might better be outsiede of the structure. Do I neeed to structure really here?  
n_subjects = length(s_n_list);


analysisVars = structure_analysis_variables (partners,n_partners,s_n_list, n_subjects); 
